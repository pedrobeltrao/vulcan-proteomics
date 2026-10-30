#!/usr/bin/env python3
"""Cross-check readable (non-sealed) PRIDE-sourced files against the SHA-1 published by the PRIDE archive API v3.
Independent of our own sha256 (catches a resumed download whose prefix came from a corrupt earlier run).
Sealed raw files are skipped (mode 000; never read).  Zip-member truth rows are skipped.  Output: TSV on stdout."""
import csv, hashlib, json, os, sys, urllib.request
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(__file__))
from verify_data import rows_of, dest, DATA
UA = {"User-Agent": "pride-sha1/1.0 (" + os.environ.get("VULCAN_CONTACT_EMAIL", "noreply@example.com") + " research)"}
def api(acc):
    out, page = [], 0
    while True:
        r = json.load(urllib.request.urlopen(urllib.request.Request(f"https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{acc}/files?pageSize=100&page={page}", headers=UA), timeout=120))
        out += r
        if len(r) < 100: break
        page += 1
    return {f["fileName"]: f["checksum"] for f in out}
def sha1(p):
    h = hashlib.sha1()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""): h.update(b)
    return h.hexdigest()
rows = [r for m in ("tier_s", "dev_sealed") for r in rows_of(m) if r["accession"].startswith("PXD") and not r["member"]
        and not (r["role"] == "sealed" and not r["id"].endswith("-truth"))]
chk = {a: api(a) for a in sorted({r["accession"] for r in rows})}
def one(r):
    p = dest(r); exp = chk[r["accession"]].get(r["filename"])
    if not os.path.exists(p): return (r["accession"], r["filename"], "MISSING", "", exp or "")
    s = sha1(p)
    return (r["accession"], r["filename"], "OK" if s == exp else ("NO_PRIDE_CHECKSUM" if not exp else "MISMATCH"), s, exp or "")
with ThreadPoolExecutor(16) as ex: res = list(ex.map(one, rows))
print("accession\tfilename\tstatus\tsha1_local\tsha1_pride")
for r in res: print("\t".join(r))
sys.exit(0 if all(r[2] == "OK" for r in res) else 1)
