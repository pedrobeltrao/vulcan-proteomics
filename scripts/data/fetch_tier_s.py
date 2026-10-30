#!/usr/bin/env python3
"""Fetch Tier S benchmark data, driven by data/manifests/tier_s.tsv.  Stdlib only.

Manifest columns: id, accession, filename, url, expected_size, purpose, role, set_id, member
  * role is dev | sealed (DESIGN.md 6.6); set_id is the D*/T* id.  member (optional) turns the row into a
    single-member extraction from a remote .zip (HTTP Range on the central directory; expected_size is the
    uncompressed size), used to pull small truth tables out of multi-GB PRIDE SEARCH archives.
  * sealed raw rows are written to data/sealed/<acc>/ and are chmod 000 right after their sha256 is recorded;
    sealed truth rows go to data/truth/sealed/<acc>/ and stay readable.  This script NEVER reads a sealed file
    after the checksum step (size via stat only); the sha256 of sealed files lives in data/manifests/sealed.sha256.
  * url / expected_size may be empty for PRIDE rows: they are resolved from the PRIDE
    archive API v3 (/projects/<ACC>/files, paginated).  `--resolve` writes the resolved
    values back into the manifest so later runs need no API call.
  * dev rows whose id ends in "-truth" are written to data/truth/<accession>/, other dev rows to
    data/raw/<accession>/ (sealed rows: see above).
Behaviour: sequential, resumable (HTTP Range on <file>.part), size-verified, idempotent
(complete files with the right size are skipped), writes data/manifests/tier_s.sha256
(`sha256  path-relative-to-data/`, sha256sum -c compatible from within data/).
Usage: fetch_tier_s.py [--manifest tier_s|dev_sealed] [--roles dev|sealed] [--resolve] [--only S1,S2] [--no-hash] [--check]
Wrappers: fetch_dev.py (dev_sealed.tsv, dev rows), fetch_sealed.py (dev_sealed.tsv, sealed rows).
"""
import argparse, csv, fcntl, hashlib, json, os, struct, sys, threading, time, urllib.request, urllib.error, urllib.parse, zlib
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DATA = os.path.join(ROOT, "data")
MANIFEST = os.path.join(DATA, "manifests", "tier_s.tsv")
SHA_OUT = os.path.join(DATA, "manifests", "tier_s.sha256")
SEALED_SHA = os.path.join(DATA, "manifests", "sealed.sha256")
COLS = ["id", "accession", "filename", "url", "expected_size", "purpose", "role", "set_id", "member"]
PRIDE = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/%s/files"
UA = {"User-Agent": "tier-s-fetch/1.0 (" + os.environ.get("VULCAN_CONTACT_EMAIL", "noreply@example.com") + " research)"}


def http(url, headers=None, method="GET", timeout=120):
    req = urllib.request.Request(url, headers={**UA, **(headers or {})}, method=method)
    return urllib.request.urlopen(req, timeout=timeout)


def retry(fn, n=5):
    for i in range(n):
        try:
            return fn()
        except (urllib.error.URLError, ConnectionError, TimeoutError, OSError) as e:
            if isinstance(e, urllib.error.HTTPError) and e.code in (403, 404):
                raise
            print(f"  retry {i+1}/{n}: {e}", file=sys.stderr)
            time.sleep(3 * (i + 1))
    raise RuntimeError("giving up")


def pride_files(acc, cache={}):
    if acc not in cache:
        out, page = [], 0
        while True:
            r = json.load(retry(lambda: http(PRIDE % acc + f"?pageSize=100&page={page}")))
            out += r
            if len(r) < 100:
                break
            page += 1
        cache[acc] = {f["fileName"]: f for f in out}
    return cache[acc]


def zip_dir(url):
    """Return {member: (method, comp_size, uncomp_size, local_header_offset)} from a remote zip (non-zip64 only)."""
    n = int(retry(lambda: http(url, method="HEAD")).headers["Content-Length"])
    rng = lambda a, b: retry(lambda: http(url, {"Range": f"bytes={a}-{b}"})).read()
    t = rng(max(0, n - 70000), n - 1)
    i = t.rfind(b"PK\x05\x06")
    cdsz, cdoff = struct.unpack("<LL", t[i + 12:i + 20])
    cd, i, ents = rng(cdoff, cdoff + cdsz - 1), 0, {}
    while cd[i:i + 4] == b"PK\x01\x02":
        (_, _, _, _, meth, _, _, _, cs, us, nl, el, cl, _, _, _, off) = struct.unpack("<4sHHHHHHLLLHHHHHLL", cd[i:i + 46])
        ents[cd[i + 46:i + 46 + nl].decode()] = (meth, cs, us, off)
        i += 46 + nl + el + cl
    return ents


def zip_member(row):
    ents = zip_dir(row["url"])
    if row["member"] not in ents:
        raise SystemExit(f"{row['id']}: member {row['member']} not in {row['url']}")
    return ents[row["member"]]


def read_manifest():
    with open(MANIFEST, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_manifest(rows):
    import io
    buf = io.StringIO(newline="")
    w = csv.DictWriter(buf, COLS, delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(rows)
    if os.path.exists(MANIFEST) and open(MANIFEST, newline="").read() == buf.getvalue():
        return  # unchanged: do not touch (parallel runs / git status stay quiet)
    tmp = MANIFEST + f".tmp{os.getpid()}"
    with open(tmp, "w", newline="") as fh: fh.write(buf.getvalue())
    os.replace(tmp, MANIFEST)


def resolve(row):
    if row["url"] and row["expected_size"]:
        return
    acc = row["accession"]
    if row.get("member"):
        if not row["url"]:
            raise SystemExit(f"{row['id']}: member rows need an explicit zip url")
        row["expected_size"] = str(zip_member(row)[2]); return
    if not acc.startswith("PXD"):
        if not row["url"]:
            raise SystemExit(f"{row['id']}: non-PRIDE rows need an explicit url")
        if not row["expected_size"]:
            row["expected_size"] = http(row["url"], method="HEAD").headers["Content-Length"]
        return
    f = pride_files(acc).get(row["filename"])
    if f is None:
        raise SystemExit(f"{row['id']}: {row['filename']} not found in {acc} (file list changed?)")
    if not row["url"]:
        ftp = [l["value"] for l in f["publicFileLocations"] if "FTP" in l["name"]][0]
        row["url"] = ftp.replace("ftp://ftp.pride.ebi.ac.uk/", "https://ftp.pride.ebi.ac.uk/")
    row["expected_size"] = str(f["fileSizeBytes"])


def is_sealed_raw(row):
    return row.get("role") == "sealed" and not row["id"].endswith("-truth")


def dest(row):
    truth = row["id"].endswith("-truth")
    if row.get("role") == "sealed":
        base = os.path.join(DATA, "truth", "sealed") if truth else os.path.join(DATA, "sealed")
    else:
        base = os.path.join(DATA, "truth" if truth else "raw")
    return os.path.join(base, row["accession"], row["filename"])


def download(row):
    path, size = dest(row), int(row["expected_size"])
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if os.path.exists(path) and os.path.getsize(path) == size:
        print(f"[skip] {row['id']} {row['filename']} ({size} B)"); finalize(row, path); return
    if row.get("member"):
        meth, cs, us, off = zip_member(row)
        h = retry(lambda: http(row["url"], {"Range": f"bytes={off}-{off + 29}"})).read()
        nl, el = struct.unpack("<HH", h[26:30])
        raw = retry(lambda: http(row["url"], {"Range": f"bytes={off + 30 + nl + el}-{off + 30 + nl + el + cs - 1}"})).read()
        data = zlib.decompress(raw, -15) if meth == 8 else raw
        if len(data) != size:
            raise SystemExit(f"FAILED size check: {row['filename']}")
        with open(path + ".part", "wb") as out: out.write(data)
        os.rename(path + ".part", path); print(f"[done] {row['filename']} (zip member)"); return
    if os.path.exists(path):  # wrong size: do not trust the name; keep the bytes only as a resumable prefix
        os.rename(path, path + ".part")
    if os.path.exists(path + ".part") and not os.access(path + ".part", os.W_OK):
        os.chmod(path + ".part", 0o644)  # an incomplete sealed download is not yet sealed data
    part = path + ".part"
    for attempt in range(8):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if have > size:
            os.remove(part); have = 0
        if have < size:
            hdr = {"Range": f"bytes={have}-"} if have else {}
            print(f"[get ] {row['id']} {row['filename']} from {have}/{size}", flush=True)
            try:
                r = http(row["url"], hdr)
                if have and r.status != 206:
                    have = 0  # server ignored Range
                with open(part, "ab" if have else "wb") as out:
                    while True:
                        b = r.read(1 << 20)
                        if not b: break
                        out.write(b)
            except (urllib.error.URLError, ConnectionError, TimeoutError, OSError) as e:
                print(f"  interrupted: {e}", file=sys.stderr); time.sleep(5 * (attempt + 1)); continue
        if os.path.getsize(part) == size:
            if not pride_sha1_ok(row, part):
                print(f"  PRIDE sha1 mismatch for {row['filename']}: discarding and re-fetching from 0", file=sys.stderr)
                os.remove(part); continue
            os.rename(part, path); print(f"[done] {row['filename']}"); finalize(row, path); return
    raise SystemExit(f"FAILED size check: {row['filename']}")


def pride_sha1_ok(row, path):
    """Verify a completed download against the SHA-1 the PRIDE archive publishes (absent for old projects: then size only)."""
    if not row["accession"].startswith("PXD") or row.get("member"):
        return True
    try:
        exp = pride_files(row["accession"]).get(row["filename"], {}).get("checksum")
    except Exception as e:  # API trouble must not block a size-verified file
        print(f"  (no PRIDE checksum: {e})", file=sys.stderr); return True
    if not exp:
        return True
    h = hashlib.sha1()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""): h.update(b)
    return h.hexdigest() == exp


CACHE_P = os.path.join(DATA, "manifests", ".sha_cache.json")
_LOCK = threading.Lock()


def sha_cached(path):
    """sha256 via cache keyed on size:mtime.  Unreadable (sealed, mode 000) files must already be cached."""
    st = os.stat(path); key = f"{st.st_size}:{int(st.st_mtime)}"; rel = os.path.relpath(path, DATA)
    load = lambda: json.load(open(CACHE_P)) if os.path.exists(CACHE_P) else {}
    c = load().get(rel)
    if c and c.get("key") == key:
        return c["sha"]
    if not os.access(path, os.R_OK):
        raise SystemExit(f"{rel}: unreadable and not in sha cache; cannot checksum without unsealing")
    print(f"[hash] {rel}", flush=True); sha = sha256(path)
    with _LOCK, open(CACHE_P + ".lock", "w") as lk:  # serialise read-modify-write across threads and processes
        fcntl.flock(lk, fcntl.LOCK_EX)
        cache = load(); cache[rel] = {"key": key, "sha": sha}
        tmp = CACHE_P + f".{os.getpid()}.{threading.get_ident()}"
        json.dump(cache, open(tmp, "w")); os.replace(tmp, CACHE_P)
    return sha


def finalize(row, path):
    """Sealed raw files: checksum (while still readable), then chmod 000."""
    if is_sealed_raw(row):
        sha_cached(path); os.chmod(path, 0o000); print(f"[seal] {os.path.relpath(path, DATA)} chmod 000")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(1 << 22)
            if not b: break
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--resolve", action="store_true", help="fill url/expected_size in manifest and exit")
    ap.add_argument("--only", help="comma-separated ids (prefix match on id before '-', e.g. S1,S5)")
    ap.add_argument("--no-hash", action="store_true")
    ap.add_argument("--manifest", default="tier_s", help="manifest name under data/manifests/ (default tier_s)")
    ap.add_argument("--roles", help="comma-separated roles to fetch (dev,sealed); default all")
    ap.add_argument("--jobs", type=int, default=1, help="parallel downloads (threads); PRIDE FTP gives 1-4 MB/s per connection")
    ap.add_argument("--check", action="store_true", help="verify existing files vs sizes only")
    a = ap.parse_args()
    global MANIFEST, SHA_OUT
    MANIFEST = os.path.join(DATA, "manifests", a.manifest + ".tsv")
    SHA_OUT = os.path.join(DATA, "manifests", a.manifest + ".sha256")
    rows = read_manifest()
    for r in rows: resolve(r)
    write_manifest(rows)
    if a.resolve: return
    sel = set(a.only.split(",")) if a.only else None
    roles = set(a.roles.split(",")) if a.roles else None
    rows = [r for r in rows if (not sel or r["id"].split("-")[0] in sel) and (not roles or r["role"] in roles)]
    if a.jobs > 1 and not a.check:
        def safe(r):
            try: download(r); return None
            except BaseException as e: return f"{r['id']} {r['filename']}: {e}"
        with ThreadPoolExecutor(a.jobs) as ex: errs = [e for e in ex.map(safe, rows) if e]
        if errs: raise SystemExit("FAILED:\n  " + "\n  ".join(errs))
        rows = [] if not a.check else rows
    for r in rows:
        if a.check:
            p = dest(r); ok = os.path.exists(p) and os.path.getsize(p) == int(r["expected_size"])
            print(("OK  " if ok else "MISS"), r["id"], r["filename"]); continue
        download(r)
    if a.no_hash or a.check: return
    lines, sealed = [], []
    for r in read_manifest():
        p = dest(r)
        if not os.path.exists(p): continue
        line = f"{sha_cached(p)}  {os.path.relpath(p, DATA)}"
        (sealed if is_sealed_raw(r) else lines).append(line)
    srt = lambda L: "\n".join(sorted(L, key=lambda l: l.split("  ")[1])) + "\n"
    tmp = SHA_OUT + f".tmp{os.getpid()}"
    with open(tmp, "w") as fh: fh.write(srt(lines))
    os.replace(tmp, SHA_OUT)
    print(f"wrote {SHA_OUT} ({len(lines)} files)")
    if sealed:
        # merge with other manifests' sealed entries
        with _LOCK, open(SEALED_SHA + ".lock", "w") as lk:
            fcntl.flock(lk, fcntl.LOCK_EX)
            old = {l.split("  ")[1]: l for l in open(SEALED_SHA).read().split("\n") if l} if os.path.exists(SEALED_SHA) else {}
            old.update({l.split("  ")[1]: l for l in sealed})
            tmp = SEALED_SHA + f".tmp{os.getpid()}"
            with open(tmp, "w") as fh: fh.write(srt(list(old.values())))
            os.replace(tmp, SEALED_SHA)
        print(f"wrote {SEALED_SHA} ({len(old)} sealed files; unreadable by design: verify only via bench/sealed_run.py)")

if __name__ == "__main__":
    main()
