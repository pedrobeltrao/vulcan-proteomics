#!/usr/bin/env python3
"""A6.3 HN1 selection: frame + exclusions (metadata only), the seeded walk with cluster slots, per-project listing /
file guard / file order, and the label-free stopping rule on capped counts.  Admission (A6.4) is an injected callback
(scripts/a6/admit_hn1.py provides it); this module never reads a label.

Phases:
  --exclusions   (metadata only; allowed before the tag) build config/a6_excluded_accessions.tsv (A6.3.2 items 1-6) and
                 check the pool size against reviews/checks/a6_power.json (1,415 at drafting).
  --walk         (agent S, step 4; ONLY after the prereg-v0.4.5 tag) the walk.  Refuses if the tag is absent: the walk
                 seed S = int(first 8 hex of the tag commit sha, 16) cannot be computed before the tag, and no PRIDE
                 listing is made before it.

Reuse: reviews/checks/a6_power.pool (frame, exclusions, person/lab keys), scripts/e9/select_e9 (INSTR, BAD_KW,
BAD_FILE, pride_files, ftp_url), scripts/e9/stream.download (size + PRIDE SHA-1 verified), cap_factors (a4_power,
verbatim; or scripts/a6/a6_stats if it exports it).

File guard (A6.3.4(b)), stricter reading recorded (A6.1.4): a RAW is unusable if its PRIDE SHA-1 equals any SHA-1 in a
committed manifest (pride_sha1 columns), its (file name, size) equals any manifest row, or its file name equals the
basename of any file in a committed *.sha256 file (those carry no size, so the name alone is matched), or it matches a
RAW already selected in HN1; after download its sha256 must not occur in any committed *.sha256.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import glob
import hashlib
import json
import os
import random
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
for p in (ROOT, HERE, os.path.join(ROOT, "reviews", "checks"), os.path.join(ROOT, "scripts", "e9"),
          os.path.join(ROOT, "scripts", "c9")):
    if p not in sys.path:
        sys.path.insert(0, p)

import numpy as np  # noqa: E402

TAG = "prereg-v0.4.5"
CET = dt.timezone(dt.timedelta(hours=1))
DEADLINE = dt.datetime(2026, 10, 28, 16, 30, tzinfo=CET)
EXCL_TSV = os.path.join(ROOT, "config", "a6_excluded_accessions.tsv")
OUT_DIR = os.path.join(ROOT, "results", "2026-10-28_A6_hn1")
SEL_LOG = os.path.join(OUT_DIR, "selection_log.tsv")
HN1_TSV = os.path.join(ROOT, "data", "manifests", "hn1.tsv")
RUNS_PER_PROJECT, PER_CLUSTER, HARD_CAP = 3, 5, 400
SIZE_MIN, SIZE_MAX = 0.1e9, 2.0e9
CAP_SHARE, CAP_HASH_SEED = 0.10, 47
STOP = {"J_min": 20, "Jeff_min": 10, "proxy_min": 1200, "A_cand_floor": 120_000, "g1_mult": 10}
G1_NMIN_FLOOR = 11_100
LOG_COLS = ["walk_pos", "accession", "event", "cluster_size_before", "n_listed_raw", "n_usable_raw", "n_downloaded",
            "n_runs_admitted", "A", "A_cand", "proxy", "cum_runs", "J", "capped_proxy", "capped_A_cand",
            "Jeff_proxy", "Jeff_A_cand", "wall_s", "time_cet", "note"]


# ------------------------------------------------------------------------------------------------ cap factors
try:  # prefer the A6 statistics module if it exports the (verbatim) copy
    from a6_stats import cap_factors  # type: ignore  # noqa: E402
except Exception:  # noqa: BLE001
    def cap_factors(c, share):
        """Copied verbatim from reviews/checks/a4_power.cap_factors (A4.3.9 cluster cap)."""
        c = np.asarray(c, float)
        if c.sum() <= 0 or c.size * share < 1:
            return np.ones_like(c)
        x = c.copy()
        for _ in range(200):
            L = share * x.sum()
            if (x <= L * (1 + 1e-9)).all():
                break
            x = np.minimum(x, L)
        return np.where(c > 0, x / np.maximum(c, 1e-12), 1.0)


def jeff(x):
    x = np.asarray(x, float)
    return float(x.sum() ** 2 / (x ** 2).sum()) if x.sum() > 0 else 0.0


def capped3(cp, ca, cA, share=CAP_SHARE):
    """A6.3.5: per-cluster factors = min of the proxy, A_cand and |A| factors; None if J * share < 1."""
    cp, ca, cA = (np.asarray(v, float) for v in (cp, ca, cA))
    if len(cp) * share < 1:
        return None
    f = np.minimum.reduce([cap_factors(cp, share), cap_factors(ca, share), cap_factors(cA, share)])
    return cp * f, ca * f, cA * f, f


def v1_minimum(cp, ca, cA, g1_nmin=G1_NMIN_FLOOR) -> tuple[bool, dict]:
    """A6.3.6 stop 1 on CAPPED counts (label-free)."""
    r = capped3(cp, ca, cA)
    a_min = max(STOP["A_cand_floor"], STOP["g1_mult"] * int(g1_nmin))
    if r is None:
        return False, {"J": len(cp), "capped_proxy": 0.0, "capped_A_cand": 0.0, "Jeff_proxy": 0.0,
                       "Jeff_A_cand": 0.0, "A_cand_min": a_min}
    xp, xa, _, f = r
    d = {"J": len(cp), "capped_proxy": float(xp.sum()), "capped_A_cand": float(xa.sum()), "Jeff_proxy": jeff(xp),
         "Jeff_A_cand": jeff(xa), "A_cand_min": a_min, "factors": [float(x) for x in f]}
    ok = (d["J"] >= STOP["J_min"] and d["Jeff_proxy"] >= STOP["Jeff_min"] and d["Jeff_A_cand"] >= STOP["Jeff_min"]
          and d["capped_proxy"] >= STOP["proxy_min"] and d["capped_A_cand"] >= a_min)
    return ok, d


def cap_keep(run: str, store_index: int, f_j: float) -> bool:
    """A6.3.5: cluster j's A spectrum is kept iff u(sha256(run | store index | 47)) < f_j (label-free)."""
    h = hashlib.sha256(f"{run}|{int(store_index)}|{CAP_HASH_SEED}".encode()).digest()
    return int.from_bytes(h[:8], "big") / 2.0 ** 64 < f_j


# ------------------------------------------------------------------------------------------------ tag / seed
def tag_sha(tag: str = TAG) -> str | None:
    r = subprocess.run(["git", "-C", ROOT, "rev-parse", "--verify", "--quiet", tag + "^{commit}"],
                       capture_output=True, text=True)
    s = r.stdout.strip()
    return s if r.returncode == 0 and re.fullmatch(r"[0-9a-f]{40}", s) else None


def walk_seed(sha: str | None) -> int:
    """S = int(first 8 hex digits of the prereg-v0.4.5 tag commit sha, 16).  Refuses without a sha."""
    if not sha or not re.fullmatch(r"[0-9a-f]{40}", sha):
        raise SystemExit(f"select_hn1: no {TAG} tag commit: the walk seed cannot be derived before the tag (A6.3.3)")
    return int(sha[:8], 16)


def walk_order(pool_accs, S: int) -> list[str]:
    order = sorted(pool_accs)
    random.Random(S).shuffle(order)
    return order


def acc_digits(acc: str) -> int:
    return int(re.sub(r"\D", "", acc))


def file_order(files: list[dict], S: int, acc: str) -> list[dict]:
    """A6.3.4(c): usable RAWs ordered by random.Random(S + accession digits) (shuffle of the name-sorted list)."""
    out = sorted(files, key=lambda f: f["fileName"])
    random.Random(S + acc_digits(acc)).shuffle(out)
    return out


# ------------------------------------------------------------------------------------------------ pool / exclusions
def build_pool():
    """-> (pool accessions, person/lab keysets, info) exactly as reviews/checks/a6_power.pool (A6.3.1-A6.3.2)."""
    import a6_power as A6P
    keep, keys, _usable, info = A6P.pool()
    return keep, keys, info


def tracked_pxd_a6():
    """A6.3.2 item 4 = reviews/checks/a6_power.tracked_pxd, skipping this module's own outputs (the exclusions file
    and the HN1 manifest/log only hold accessions that are already excluded or walked), so that the exclusions file is
    a fixed point: re-running after it is committed reproduces it byte for byte."""
    own = {"t4_frame.tsv", os.path.basename(EXCL_TSV), os.path.basename(EXCL_TSV) + ".sha256"}
    files = subprocess.run(["git", "ls-files", "data/manifests", "results", "config"], cwd=ROOT, capture_output=True,
                           text=True, check=True).stdout.split()
    acc = set()
    for f in files:
        if os.path.basename(f) in own:
            continue
        try:
            acc |= set(re.findall(r"PXD\d{6}", open(os.path.join(ROOT, f), errors="ignore").read()))
        except (OSError, UnicodeDecodeError):
            pass
    return acc


def exclusions_detail():
    """Every excluded accession with its reasons (A6.3.2 items 1-6; metadata only)."""
    import io
    import contextlib
    import pandas as pd
    import a4_r3_walk as R3
    import a6_power as A6P
    import c9_power as CP
    import pool_count as PCNT
    PCNT.OUT = os.devnull
    reasons: dict[str, set] = {}

    def add(accs, why):
        for a in accs:
            if isinstance(a, str) and re.fullmatch(r"PXD\d{6}", a):
                reasons.setdefault(a, set()).add(why)
    for f in sorted(glob.glob(os.path.join(ROOT, "data", "manifests", "*.tsv"))):
        if os.path.basename(f) == "t4_frame.tsv":
            continue
        try:
            for r in csv.DictReader(open(f), delimiter="\t"):
                for k in ("accession", "Accession", "project"):
                    add([r.get(k) or ""], "1:manifest:" + os.path.basename(f))
        except Exception:  # noqa: BLE001
            pass
    for f in glob.glob(os.path.join(ROOT, "results", "**", "t4_selector", "**", "t4_selection.tsv"), recursive=True):
        add([r["accession"] for r in csv.DictReader(open(f), delimiter="\t")], "2:t4_draw")
    add([t[0] for t in json.load(open(os.path.join(ROOT, "data", "manifests", "e9_d8.log.json")))["tried"]],
        "3:e9_tried")
    add(["PXD001468"], "3:D0")
    add(tracked_pxd_a6(), "4:tracked_file")
    add([re.match(r"files_(PXD\d{6})", os.path.basename(f)).group(1)
         for f in glob.glob(os.path.join(ROOT, "data/cache/pride_meta/files_PXD*"))], "5:cached_listing")
    # item 6: same-lab rule, evaluated over the C9 frame pool (as a6_power.pool)
    with contextlib.redirect_stdout(io.StringIO()):
        projects_c9, _, _ = CP.pool()
    crawl = pd.read_csv(os.path.join(ROOT, "data/cache/pride_meta/all_projects_classified.tsv"), sep="\t",
                        dtype=str).fillna("").set_index("accession")
    prior = R3.prior_null_accessions()
    kp = R3.keysets(sorted(prior), crawl)
    prior_keys = set().union(*kp.values()) if kp else set()
    cand = sorted(a for a, _ in projects_c9)
    kc = R3.keysets(cand, crawl)
    add([a for a in cand if kc[a] & prior_keys], "6:same_lab")
    return reasons


def write_exclusions(path: str = EXCL_TSV) -> dict:
    reasons = exclusions_detail()
    with open(path, "w") as fh:
        fh.write("accession\treasons\n")
        for a in sorted(reasons):
            fh.write(f"{a}\t{';'.join(sorted(reasons[a]))}\n")
    sha = hashlib.sha256(open(path, "rb").read()).hexdigest()
    with open(path + ".sha256", "w") as fh:
        fh.write(f"{sha}  config/{os.path.basename(path)}\n")
    keep, _keys, info = build_pool()
    excl = set(reasons)
    assert not (set(keep) & excl), "pool and exclusions overlap"
    exp = json.load(open(os.path.join(ROOT, "reviews", "checks", "a6_power.json")))["pool"]["a6_pool"]
    return {"n_excluded": len(excl), "sha256": sha, "pool": len(keep), "pool_at_drafting": exp,
            "pool_matches_drafting": len(keep) == exp, "pool_info": info}


# ------------------------------------------------------------------------------------------------ file guard
def known_files():
    """Committed manifests: PRIDE SHA-1s, (name, size) pairs, *.sha256 basenames and sha256 values."""
    files = subprocess.run(["git", "-C", ROOT, "ls-files", "*.sha256", "*/cache_raw.sha256", "data/manifests/*.tsv"],
                           capture_output=True, text=True, check=True).stdout.split()
    sha1, name_size, names, sha256 = set(), set(), set(), set()
    for f in files:
        p = os.path.join(ROOT, f)
        if f.endswith(".sha256"):
            for line in open(p, errors="ignore"):
                parts = line.split()
                if len(parts) >= 2:
                    sha256.add(parts[0].lower())
                    names.add(os.path.basename(parts[-1]))
        else:
            try:
                for r in csv.DictReader(open(p, errors="ignore"), delimiter="\t"):
                    if r.get("pride_sha1"):
                        sha1.add(r["pride_sha1"].lower())
                    fn = r.get("filename") or r.get("file_name") or r.get("fileName")
                    sz = r.get("expected_size") or r.get("size") or r.get("fileSizeBytes")
                    if fn and sz and str(sz).isdigit():
                        name_size.add((os.path.basename(fn), int(sz)))
            except Exception:  # noqa: BLE001
                pass
    return {"sha1": sha1, "name_size": name_size, "names": names, "sha256": sha256}


def pride_sha1(f: dict) -> str:
    for c in f.get("checksums") or []:
        if isinstance(c, dict) and "sha" in str(c.get("type", c.get("name", ""))).lower():
            return str(c.get("value", "")).lower()
    return str(f.get("checksum", "") or "").lower()


def usable_raws(files: list[dict], guard: dict, selected: set) -> tuple[list[dict], list[dict]]:
    """A6.3.4(a)-(b): .raw, 0.1-2.0 GB, not BAD_FILE, has an https URL, passes the file guard."""
    import select_e9 as SE
    ok, rej = [], []
    for f in files:
        name, size = f.get("fileName", ""), int(f.get("fileSizeBytes") or 0)
        why = None
        if f.get("fileCategory", {}).get("value") != "RAW" or not name.lower().endswith(".raw"):
            why = "not raw"
        elif not SIZE_MIN <= size <= SIZE_MAX:
            why = "size"
        elif SE.BAD_FILE.search(name):
            why = "BAD_FILE"
        elif not SE.ftp_url(f):
            why = "no url"
        else:
            s1 = pride_sha1(f)
            if (s1 and s1 in guard["sha1"]) or (name, size) in guard["name_size"] or name in guard["names"] \
                    or (s1 and s1 in selected) or (name, size) in selected:
                why = "file guard"
        (rej if why else ok).append(f if not why else {"fileName": name, "why": why})
    return ok, rej


def list_project(acc: str, sha: str | None):
    """A6.3.4(a) listing (PRIDE API).  Only after the tag: refuses otherwise.  NEVER called before the tag."""
    walk_seed(sha)
    import select_e9 as SE
    return SE.pride_files(acc)


def download_raw(row: dict, work: str) -> str | None:
    """scripts/e9/stream.download (size + PRIDE SHA-1 verified) with the sha256 log redirected to HN1's own file
    (results/2026-10-28_A6_hn1/hn1_raw.sha256) so no E9 manifest is touched.  row: id, filename, url, expected_size,
    pride_sha1.  -> local path or None."""
    import queue
    import stream as ST
    os.makedirs(os.path.join(work, "raw"), exist_ok=True)
    ST.SHA = os.path.join(OUT_DIR, "hn1_raw.sha256")
    ready, failed = queue.Queue(), []
    ST.download(row, work, ST.Budget(1 << 62), ready, failed)
    if failed or ready.empty():
        return None
    return ready.get()[1]


# ------------------------------------------------------------------------------------------------ the walk
class Clusters:
    """Union-find over ADMITTED projects linked by any shared person/lab key (A6.3.3)."""

    def __init__(self):
        self.par, self.tot, self.key_owner = [], {}, {}

    def find(self, i):
        while self.par[i] != i:
            self.par[i] = self.par[self.par[i]]
            i = self.par[i]
        return i

    def linked(self, ks):
        return {self.find(self.key_owner[k]) for k in ks if k in self.key_owner}

    def size_if_added(self, ks):
        return sum(self.tot[r][3] for r in self.linked(ks)) + 1

    def admit(self, ks, proxy, a_cand, a):
        linked = self.linked(ks)
        i = len(self.par)
        self.par.append(i)
        t = [proxy, a_cand, a, 1]
        for r in linked:
            self.par[r] = i
            o = self.tot.pop(r)
            t = [t[j] + o[j] for j in range(4)]
        self.tot[i] = t
        for k in ks:
            self.key_owner[k] = i

    def arrays(self):
        v = np.array(list(self.tot.values()), float).reshape(-1, 4)
        return v[:, 0], v[:, 1], v[:, 2]


def _now():
    return dt.datetime.now(CET)


def run_walk(order, keys, admit, list_fn, *, S, log_path=SEL_LOG, deadline=DEADLINE, hard_cap=HARD_CAP,
             g1_nmin=G1_NMIN_FLOOR, now=_now, guard=None):
    """The walk (A6.3.3-A6.3.6).  admit(acc, ordered_usable_files, max_runs) -> {"runs": [{"run", "A", "A_cand",
    "proxy"}...], "n_downloaded": int, "note": str}; a project is admitted iff >= 1 run passes.  list_fn(acc) ->
    PRIDE file records.  Returns the walk summary; every step is written to the selection log."""
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    guard = guard if guard is not None else known_files()
    cl, selected, listed = Clusters(), set(), []
    st = {"runs": 0, "admitted": [], "end": "pool exhausted", "walked": 0, "cluster_full": 0, "projects": {}}
    fh = open(log_path, "w")
    fh.write("\t".join(LOG_COLS) + "\n")

    def log(**kw):
        kw.setdefault("time_cet", now().isoformat(timespec="seconds"))
        fh.write("\t".join(str(kw.get(c, "")) for c in LOG_COLS) + "\n")
        fh.flush()
    for pos, acc in enumerate(order):
        if now() >= deadline:
            st["end"] = "walk deadline"
            break
        if st["runs"] >= hard_cap:
            st["end"] = "hard cap"
            break
        st["walked"] += 1
        ks = keys.get(acc) or {"?" + acc}
        size_before = cl.size_if_added(ks) - 1
        if size_before + 1 > PER_CLUSTER:
            st["cluster_full"] += 1
            log(walk_pos=pos, accession=acc, event="cluster full", cluster_size_before=size_before)
            continue
        t0 = time.time()
        listed.append(acc)
        files = list_fn(acc)
        ok, _rej = usable_raws(files, guard, selected)
        ordered = file_order(ok, S, acc)
        res = admit(acc, ordered, min(RUNS_PER_PROJECT, hard_cap - st["runs"]))
        runs = res.get("runs", [])
        if now() >= deadline:                                  # A6.3.6 stop 3: last project "reached, not used"
            st["end"] = "walk deadline"
            st["projects"][acc] = "reached, not used"
            log(walk_pos=pos, accession=acc, event="reached, not used", n_listed_raw=len(files),
                n_usable_raw=len(ok), n_downloaded=res.get("n_downloaded", ""), wall_s=round(time.time() - t0, 1))
            break
        for f in ordered[:res.get("n_downloaded", 0)]:
            s1 = pride_sha1(f)
            if s1:
                selected.add(s1)
            selected.add((f["fileName"], int(f.get("fileSizeBytes") or 0)))
        if not runs:
            st["projects"][acc] = "rejected"
            log(walk_pos=pos, accession=acc, event="rejected", cluster_size_before=size_before,
                n_listed_raw=len(files), n_usable_raw=len(ok), n_downloaded=res.get("n_downloaded", 0),
                n_runs_admitted=0, cum_runs=st["runs"], wall_s=round(time.time() - t0, 1), note=res.get("note", ""))
            continue
        A = sum(r["A"] for r in runs)
        ac = sum(r["A_cand"] for r in runs)
        pr = sum(r["proxy"] for r in runs)
        cl.admit(ks, pr, ac, A)
        st["runs"] += len(runs)
        st["admitted"].append(acc)
        st["projects"][acc] = "admitted"
        cp, ca, cA = cl.arrays()
        ok_min, d = v1_minimum(cp, ca, cA, g1_nmin)
        log(walk_pos=pos, accession=acc, event="admitted", cluster_size_before=size_before, n_listed_raw=len(files),
            n_usable_raw=len(ok), n_downloaded=res.get("n_downloaded", 0), n_runs_admitted=len(runs), A=A,
            A_cand=ac, proxy=pr, cum_runs=st["runs"], J=d["J"], capped_proxy=round(d["capped_proxy"], 1),
            capped_A_cand=round(d["capped_A_cand"], 1), Jeff_proxy=round(d["Jeff_proxy"], 2),
            Jeff_A_cand=round(d["Jeff_A_cand"], 2), wall_s=round(time.time() - t0, 1), note=res.get("note", ""))
        if ok_min:
            st["end"] = "v1 minimum"
            break
    else:
        st["end"] = "pool exhausted"
    if st["end"] == "pool exhausted" and st["runs"] >= hard_cap:
        st["end"] = "hard cap"
    fh.close()
    cp, ca, cA = cl.arrays()
    st["v1_minimum_reached"] = st["end"] == "v1 minimum"
    st["final_counts"] = v1_minimum(cp, ca, cA, g1_nmin)[1]
    st["listed"] = listed
    st["clusters"] = cl
    return st


def write_hn1_manifest(st: dict, path: str = HN1_TSV):
    """A6.3.7: every LISTED project (admitted, rejected after listing, reached-not-used) gets role dev-null."""
    with open(path, "w") as fh:
        fh.write("accession\trole\tstatus\tamendment\n")
        for a in st["listed"]:
            fh.write(f"{a}\tdev-null\t{st['projects'].get(a, 'listed')}\tA6\n")


def main():  # pragma: no cover
    ap = argparse.ArgumentParser()
    ap.add_argument("--exclusions", action="store_true")
    ap.add_argument("--walk", action="store_true")
    a = ap.parse_args()
    if a.exclusions:
        print(json.dumps({k: v for k, v in write_exclusions().items() if k != "pool_info"}, indent=1))
    if a.walk:
        S = walk_seed(tag_sha())                       # refuses before the tag
        raise SystemExit("select_hn1 --walk: run through scripts/a6/admit_hn1.py (agent S), which injects admission; "
                         f"walk seed {S}")


if __name__ == "__main__":
    main()
