#!/usr/bin/env python3
"""Operator wrapper for A2.9.4 steps 2-3: runs the UNMODIFIED tagged selector (scripts/a2/select_dn2.py main()) with
an I/O class that differs from select_dn2.PrideIO only in (results/2026-10-27_A2_dn2/DEVIATIONS.md items 1-3, 7, 9):
  * RAWs go to local NVMe ($VULCAN_A2_WORK/raw/dn2, default /tmp/vulcan_a2/raw/dn2);
  * each download is verified by PRIDE SHA-1 (API `checksum`, else the project's checksum.txt) plus size; sha256 and
    SHA-1 are computed in one pass; mismatch -> delete and re-fetch (max 3), then abort the walk;
  * every walk decision / download is appended to $W/walk_events.jsonl (label-free content only);
  * RAWs and stores of marker-screen-excluded projects are deleted after the screen.
Selection logic, order, seeds, picks, admission and outputs are exactly the tagged code's.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "a2"))
import select_dn2 as S  # noqa: E402

W = os.environ.get("VULCAN_A2_WORK", "/tmp/vulcan_a2")
RAW = os.path.join(W, "raw", "dn2")
EVENTS = os.path.join(W, "walk_events.jsonl")
RUN_NO = os.environ.get("A2_WALK_RUN", "1")


def event(**kw):
    kw = {"t": time.strftime("%Y-%m-%dT%H:%M:%S"), "run": RUN_NO, **kw}
    with open(EVENTS, "a") as fh:
        fh.write(json.dumps(kw, default=str) + "\n")
    print(json.dumps(kw, default=str), flush=True)


# DEVIATIONS item 17: admit.gpu_exhaustive crashes on a run with ZERO certified spectra (engine.search_exhaustive
# returns k2=None for an empty candidate set).  Guard: such a run returns empty arrays without searching, so every
# spectrum keeps code 0 (not certified), exactly what admit_project computes for n = 0.  Non-empty runs: unchanged call.
import numpy as _np  # noqa: E402
import admit as ADM  # noqa: E402

_orig_gpu_exhaustive = ADM.gpu_exhaustive


def _gpu_exhaustive_guarded(idx_dir, stores):
    f = _orig_gpu_exhaustive(idx_dir, stores)

    def run_one(run, cert):
        if len(cert["spec"]) == 0:
            event(ev="empty_certified_run", run_name=run)
            return {"max_hs": _np.zeros(0), "top_seq": [], "partner_seq": []}
        return f(run, cert)
    return run_one


ADM.gpu_exhaustive = _gpu_exhaustive_guarded


def _earlier_sha256(acc, fname):
    if not os.path.exists(EVENTS):
        return None
    for line in open(EVENTS):
        e = json.loads(line)
        if e.get("ev") == "downloaded" and e.get("accession") == acc and e.get("filename") == fname:
            return e["sha256"]
    return None


def hashes(p):
    h256, h1 = hashlib.sha256(), hashlib.sha1()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h256.update(b)
            h1.update(b)
    return h256.hexdigest(), h1.hexdigest()


class OpsIO(S.PrideIO):
    def __init__(self, cfg, raw_dir):  # raw_dir of the tagged main() (data/raw/dn2) is replaced by local NVMe
        super().__init__(cfg, RAW)
        self._cksum_txt = {}

    def project_record(self, acc):
        event(ev="candidate", accession=acc)
        return super().project_record(acc)

    def list_files(self, acc):
        out = [{"fileName": f["fileName"], "fileSizeBytes": f["fileSizeBytes"], "url": self._url(f),
                "pride_sha1": (f.get("checksum") or "").strip().lower()}
               for f in self._files(acc) if f.get("fileCategory", {}).get("value") == "RAW"]
        event(ev="listed", accession=acc, n_raw_records=len(out))
        return out

    def proteome(self, row):
        p = super().proteome(row)
        event(ev="proteome", accession=row["accession"], proteome=p)
        return p

    def _checksum_txt(self, acc, url):
        if acc not in self._cksum_txt:
            m = {}
            try:
                txt = self._get(url.rsplit("/", 1)[0] + "/checksum.txt").decode(errors="replace")
                for line in txt.splitlines():
                    parts = line.split()
                    if len(parts) >= 2 and len(parts[-1]) == 40:
                        m[parts[0]] = parts[-1].lower()
                    elif len(parts) >= 2 and len(parts[0]) == 40:
                        m[parts[-1]] = parts[0].lower()
            except Exception:  # noqa: BLE001  (absent checksum.txt is allowed: recorded as no_pride_checksum)
                pass
            self._cksum_txt[acc] = m
        return self._cksum_txt[acc]

    def download(self, acc, files):
        out = []
        d = os.path.join(self.raw_dir, acc)
        os.makedirs(d, exist_ok=True)
        for f in files:
            p = os.path.join(d, f["fileName"])
            exp1 = f.get("pride_sha1") or self._checksum_txt(acc, f["url"]).get(f["fileName"], "")
            src = "api" if f.get("pride_sha1") else ("checksum.txt" if exp1 else "none")
            prev256 = _earlier_sha256(acc, f["fileName"])   # item 19: same bytes as any earlier run of this walk
            for attempt in range(6):
                if not (os.path.exists(p) and os.path.getsize(p) == f["fileSizeBytes"]):
                    if attempt > 0 and os.path.exists(p):
                        os.remove(p)
                    for k in range(10):    # DEVIATIONS item 18: transient HTTP errors (403) are retried with backoff
                        if subprocess.run(["curl", "-sSf", "--retry", "5", "-C", "-", "-o", p, f["url"]]).returncode == 0:
                            break
                        event(ev="curl_retry", accession=acc, file=f["fileName"], attempt=k)
                        time.sleep(30 * (k + 1))
                    else:
                        raise RuntimeError(f"curl failed 10 times: {f['url']}")
                if os.path.getsize(p) != f["fileSizeBytes"]:
                    event(ev="size_mismatch", accession=acc, file=f["fileName"], attempt=attempt,
                          got=os.path.getsize(p), expected=f["fileSizeBytes"])
                    os.remove(p)
                    continue
                s256, s1 = hashes(p)
                if (not exp1 or s1 == exp1) and (prev256 is None or s256 == prev256):
                    break
                event(ev="sha1_mismatch", accession=acc, file=f["fileName"], attempt=attempt)
                os.remove(p)
            else:
                event(ev="ABORT_download_unverifiable", accession=acc, file=f["fileName"])
                raise SystemExit(f"download of {acc}/{f['fileName']} could not be verified; walk aborted")
            row = {"filename": f["fileName"], "url": f["url"], "expected_size": f["fileSizeBytes"], "sha256": s256,
                   "path": p, "sha1": s1, "pride_sha1": exp1, "pride_sha1_source": src,
                   "sha1_status": "OK" if exp1 else "no_pride_checksum"}
            event(ev="downloaded", accession=acc, **{k: v for k, v in row.items() if k != "path"})
            out.append(row)
        return out

    def marker_and_admit(self, acc, files, prot):
        t0 = time.time()
        res = super().marker_and_admit(acc, files, prot)
        mk = res["marker"]
        event(ev="marker_admit", accession=acc, upid=prot["upid"], exclude=mk["exclude"],
              detected_members=mk["detected_members"], trigger_members=mk["trigger_members"],
              counts=res["counts"], seconds=round(time.time() - t0))
        if mk["exclude"] and os.environ.get("A2_DELETE_EXCLUDED", "0") == "1":   # item 18: kept until walk end
            for f in files:
                for q in (f["path"], os.path.join(W, "store_dn2", f["filename"][:-4] + ".vspec.arrow"),
                          os.path.join(W, "store_dn2", f["filename"][:-4] + ".vspec.arrow.json")):
                    if os.path.exists(q):
                        os.remove(q)
        return res


D16_FIX = "d6d199c"   # the committed D16 fix of scripts/a2/select_dn2.py (DEVIATIONS item 14)
A2_PATHS = ("scripts/a2", "vulcan/gates", "config/a2.yaml", "config/a2_panels.tsv", "config/a2_excluded_accessions.tsv")


def require_tag_d16(tag, paths=A2_PATHS):
    """Stricter replacement of select_dn2.require_tag after D16: every A2 path identical to the tag except
    scripts/a2/select_dn2.py, which must equal the D16 fix commit; tag and fix are ancestors of HEAD; clean tree."""
    def git(*a):
        return subprocess.run(["git", "-C", ROOT, *a], capture_output=True, text=True)
    for ref in (tag, D16_FIX):
        if git("merge-base", "--is-ancestor", ref, "HEAD").returncode != 0:
            raise SystemExit(f"run_dn2_walk: HEAD does not contain {ref}")
    if git("diff", "--quiet", tag, "--", *paths, ":!scripts/a2/select_dn2.py").returncode != 0:
        raise SystemExit("run_dn2_walk: A2 code/config differ from the tag outside the D16 fix")
    if git("diff", "--quiet", D16_FIX, "--", "scripts/a2/select_dn2.py").returncode != 0:
        raise SystemExit("run_dn2_walk: select_dn2.py differs from the D16 fix commit")
    if git("status", "--porcelain", "--", *paths).stdout.strip():
        raise SystemExit("run_dn2_walk: uncommitted changes under A2 paths")


def main():
    os.makedirs(W, exist_ok=True)
    if not os.path.exists(os.path.join(W, "host_threshold.json")):
        raise SystemExit("run scripts/a2/host_threshold.py first (DEVIATIONS item 5)")
    cfg = S.a2_config.load()
    excl = S.load_exclusions(cfg)
    event(ev="start", n_exclusions=len(excl),
          exclusions_sha256=hashlib.sha256("\n".join(sorted(excl)).encode()).hexdigest(),
          git_head=subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True,
                                  text=True).stdout.strip())
    S.PrideIO = OpsIO
    S.require_tag = require_tag_d16
    S.main()
    event(ev="end")


def swap(acc, out):
    """P1 swap test (admit.swap_test) with the item-17 guard active."""
    S.require_tag = require_tag_d16
    r = ADM.swap_test(acc, S.a2_config.load())
    json.dump(r, open(out, "w"), indent=1)
    event(ev="swap_test", accession=acc, result=r)


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "--swap-test":
        swap(sys.argv[2], sys.argv[3])
    else:
        main()
