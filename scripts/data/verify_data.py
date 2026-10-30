#!/usr/bin/env python3
"""Compare every manifest row with disk: raw, truth and FASTA.  Stdlib only.

Sources: data/manifests/{tier_s,dev_sealed,d9,d4_ref}.tsv (path rule as in fetch_tier_s.dest) and fasta.tsv.
Expected sha256: tier_s.sha256, dev_sealed.sha256, sealed.sha256, fasta.sha256 (+ fasta.tsv column).
Status per row:  OK | MISSING (no file, no .part) | PARTIAL (.part present or size < expected) | MISMATCH (size > expected,
sha differs, or a ZIP-member/other check fails).  Column `sha` = verified | cached | none (no sha recorded) | sealed-stat-only.
Sealed raw files (role=sealed, not -truth) are chmod 000 by design: checked with os.stat only, never chmod'ed or opened
(only bench/sealed_run.py may unlock them); their sha256 in sealed.sha256 is not re-verified here.  A sealed file whose
mode is not 000 is flagged in `note` (that means it was left unsealed).
sha cache: data/manifests/.sha_cache.json, key = path -> {key: "size:int(mtime)", sha}; a file is re-hashed only when
size or mtime changed.  Extra files on disk under data/{raw,sealed,truth,fasta} not in any manifest are listed as EXTRA.
Usage: verify_data.py [--out results/<dir>/verify.tsv] [--no-hash] [--jobs N]   exit code 0 iff no MISSING/PARTIAL/MISMATCH.
"""
import argparse, csv, hashlib, json, os, stat, sys
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DATA = os.path.join(ROOT, "data")
MAN = os.path.join(DATA, "manifests")
CACHE_P = os.path.join(MAN, ".sha_cache.json")


def rows_of(name):
    with open(os.path.join(MAN, name + ".tsv"), newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def dest(r):  # identical to fetch_tier_s.dest
    truth = r["id"].endswith("-truth")
    if r.get("role") == "sealed":
        base = os.path.join(DATA, "truth", "sealed") if truth else os.path.join(DATA, "sealed")
    else:
        base = os.path.join(DATA, "truth" if truth else "raw")
    return os.path.join(base, r["accession"], r["filename"])


def read_sha(name):
    p = os.path.join(MAN, name)
    out = {}
    if os.path.exists(p):
        for l in open(p).read().splitlines():
            if l.strip():
                s, f = l.split("  ", 1); out[f] = s
    return out


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""): h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out"); ap.add_argument("--no-hash", action="store_true"); ap.add_argument("--jobs", type=int, default=16)
    a = ap.parse_args()
    cache = json.load(open(CACHE_P)) if os.path.exists(CACHE_P) else {}
    exp = {}
    for n in ("tier_s.sha256", "dev_sealed.sha256", "sealed.sha256", "fasta.sha256", "derived_truth.sha256", "d9.sha256",
              "d4_ref.sha256"): exp.update(read_sha(n))

    items = []  # dict(kind, id, set, role, rel, size, sha, sealed)
    seen = set()
    for m in ("tier_s", "dev_sealed", "d9", "d4_ref"):   # d9 / d4_ref: M6 (2026-10-27)
        for r in rows_of(m):
            p = dest(r); rel = os.path.relpath(p, DATA)
            sealed = r["role"] == "sealed" and not r["id"].endswith("-truth")
            kind = "truth" if r["id"].endswith("-truth") else ("sealed_raw" if sealed else ("raw" if r["filename"].lower().endswith((".raw", ".zip")) else "meta"))
            items.append(dict(src=m, kind=kind, id=r["id"], set=r["set_id"], role=r["role"], rel=rel, size=int(r["expected_size"]) if r["expected_size"] else None,
                              sha=exp.get(rel), sealed=sealed))
            seen.add(rel)
    for rel, sh in sorted(read_sha("derived_truth.sha256").items()):  # outputs of make_pool_truth.py (deterministic)
        items.append(dict(src="derived_truth", kind="truth", id="derived", set="D2/T2", role="sealed" if "/sealed/" in rel else "dev", rel=rel, size=None, sha=sh, sealed=False))
        seen.add(rel)
    for r in rows_of("fasta"):
        rel = f"fasta/{r['name']}.fasta"
        items.append(dict(src="fasta", kind="fasta", id=r["name"], set="-", role="ref", rel=rel, size=int(r["bytes"]), sha=r["sha256"] or exp.get(rel), sealed=False))
        seen.add(rel)

    def check(it):
        p = os.path.join(DATA, it["rel"]); note = []
        part = p + ".part"
        try: st = os.stat(p)
        except FileNotFoundError: st = None
        if st is None:
            if os.path.exists(part):
                return it, "PARTIAL", os.path.getsize(part), "none", f".part {os.path.getsize(part)}/{it['size']}"
            return it, "MISSING", 0, "none", ""
        size = st.st_size
        if os.path.exists(part): note.append("stale .part present")
        if it["size"] is not None:
            if size < it["size"]: return it, "PARTIAL", size, "none", f"{size}/{it['size']}"
            if size > it["size"]: return it, "MISMATCH", size, "none", f"size {size} > expected {it['size']}"
        if it["sealed"]:
            if stat.S_IMODE(st.st_mode) != 0: note.append(f"NOT SEALED: mode {oct(stat.S_IMODE(st.st_mode))}")
            note.append("sha in sealed.sha256 not re-verified (stat only)" if it["sha"] else "no sha recorded")
            return it, "OK", size, "sealed-stat-only", "; ".join(note)
        if a.no_hash: return it, "OK", size, "skipped", "; ".join(note)
        key = f"{size}:{int(st.st_mtime)}"; c = cache.get(it["rel"])
        if c and c.get("key") == key: s, how = c["sha"], "cached"
        else: s, how = sha256(p), "verified"; cache[it["rel"]] = {"key": key, "sha": s}
        if it["sha"] is None: return it, "OK", size, "none", "; ".join(note + ["no sha256 recorded in any manifest"])
        if s != it["sha"]: return it, "MISMATCH", size, how, "; ".join(note + [f"sha256 {s[:12]} != expected {it['sha'][:12]}"])
        return it, "OK", size, how, "; ".join(note)

    with ThreadPoolExecutor(a.jobs) as ex: res = list(ex.map(check, items))
    tmp = CACHE_P + f".{os.getpid()}"; json.dump(cache, open(tmp, "w")); os.replace(tmp, CACHE_P)

    # extras: files on disk not in a manifest
    extra = []
    for top in ("raw", "sealed", "truth", "fasta"):
        for d, _, fs in os.walk(os.path.join(DATA, top)):
            for f in fs:
                rel = os.path.relpath(os.path.join(d, f), DATA)
                if rel not in seen and not f.endswith(".part"): extra.append(rel)
    parts = [os.path.relpath(os.path.join(d, f), DATA) for d, _, fs in os.walk(DATA) for f in fs if f.endswith(".part")]

    hdr = ["status", "kind", "id", "set", "role", "path", "size", "expected_size", "sha", "note"]
    lines = ["\t".join(hdr)]
    for it, st, size, how, note in res:
        lines.append("\t".join(map(str, [st, it["kind"], it["id"], it["set"], it["role"], it["rel"], size, it["size"] if it["size"] is not None else "", how, note])))
    for rel in sorted(extra): lines.append("\t".join(["EXTRA", "", "", "", "", rel, str(os.stat(os.path.join(DATA, rel)).st_size) if os.access(os.path.join(DATA, rel), os.R_OK) else "", "", "", "not in any manifest"]))
    if a.out:
        os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
        open(a.out, "w").write("\n".join(lines) + "\n")

    # summary table
    cnt = {}
    for it, st, *_ in res: cnt.setdefault(it["kind"], {}).setdefault(st, 0); cnt[it["kind"]][st] += 1
    sts = ["OK", "MISSING", "PARTIAL", "MISMATCH"]
    print(f"{'kind':<12}" + "".join(f"{s:>10}" for s in sts) + f"{'total':>8}")
    tot = {s: 0 for s in sts}
    for k in sorted(cnt):
        print(f"{k:<12}" + "".join(f"{cnt[k].get(s,0):>10}" for s in sts) + f"{sum(cnt[k].values()):>8}")
        for s in sts: tot[s] += cnt[k].get(s, 0)
    print(f"{'ALL':<12}" + "".join(f"{tot[s]:>10}" for s in sts) + f"{sum(tot.values()):>8}")
    for it, st, size, how, note in res:
        if st != "OK" or "NOT SEALED" in note or "no sha" in note: print(f"  {st:<8} {it['rel']}  {note}")
    if extra: print(f"EXTRA files not in a manifest: {len(extra)}"); [print("  ", e) for e in sorted(extra)]
    if parts: print("stray .part files:", *parts)
    gb = sum(s for _, st, s, *_ in res if st == "OK") / 1e9
    print(f"bytes on disk for OK rows: {gb:.2f} GB")
    return 0 if all(st == "OK" for _, st, *_ in res) else 1


if __name__ == "__main__":
    sys.exit(main())
