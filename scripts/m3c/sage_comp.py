#!/usr/bin/env python3
"""Sage 0.14.7 (+LDA, predict_rt on) and Sage+mokapot 0.10.0 comparators on the IT-CID dev sets (D4 all 12 fractions,
D10cf all 29 fractions), M1e set-up: same databases (prod = protein FASTA with Sage decoys; P2001-P2003 / C =
the M1e peptide-level entrapment DBs with rev_ decoys), fragment +-0.5 Da, precursor +-10 ppm, isotope 0..2,
TMT K + N-term static on D4.  Spectra: MGFs written from OUR stores (identical centroids to the Vulcan search;
E2/M4/M5 practice).  mokapot: --seed 11 --max_workers 1 --keep_decoys --decoy_prefix rev_ on Sage's PIN (M1e).

    python scripts/m3c/sage_comp.py mgf
    python scripts/m3c/sage_comp.py run D4|D10cf [--dbs prod P2001 ...] [--threads 16] [--lock-prod]
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

sys.path.insert(0, os.path.join(C.ROOT, "scripts", "baselines"))
import run_searches as RS  # noqa: E402

SAGE = os.path.join(C.ROOT, "tools", "bin", "sage")
MOKAPOT = os.path.join(os.environ.get("VULCAN_MOKAPOT_ENV", "/tmp/vulcan_m1e/envs/mokapot"), "bin", "mokapot")
TIMED = [sys.executable, os.path.join(C.ROOT, "scripts", "baselines", "timed.py")]
CACHE = os.path.join(C.WORK, "comp")
DBS = ["prod", "P2001", "P2002", "P2003", "C"]
RS.DATASETS["D10cf"] = {"accession": "PXD015406", "family": "human_contam", "files": [],
                        "instrument": "Orbitrap Fusion, label-free, MS2 CID ion trap", "prec_ppm": 10.0,
                        "frag_da": 0.5, "bucket_size": 32768}


def mgf(run):
    return os.path.join(C.WORK, "mgf", run + ".mgf")


def runs(ds):
    return [os.path.basename(s)[:-len(".vspec.arrow")] for s in C.DEV[ds]["stores"]()]


def write_mgfs():
    from vulcan.gpu import sage_bridge as SB
    from vulcan.gpu.engine import load_store
    os.makedirs(os.path.join(C.WORK, "mgf"), exist_ok=True)
    for ds in C.DEV:
        for s in C.DEV[ds]["stores"]():
            run = os.path.basename(s)[:-len(".vspec.arrow")]
            if not os.path.exists(mgf(run)):
                SB.write_mgf(load_store(s), mgf(run))
                print("mgf", run, flush=True)


def run(ds, dbs, threads, lock_prod):
    os.makedirs(os.path.join(C.RES, "configs"), exist_ok=True)
    for db in dbs:
        out = os.path.join(CACHE, "sage", ds, db)
        os.makedirs(out, exist_ok=True)
        lock = lock_prod and db == "prod"
        if not os.path.exists(os.path.join(out, "results.sage.tsv")):
            cfg = RS.sage_config(ds, db, out)
            cfg["mzml_paths"] = [mgf(r) for r in runs(ds)]
            cpath = os.path.join(C.RES, "configs", f"sage_{ds}_{db}.json")
            json.dump(cfg, open(cpath, "w"), indent=1)
            cmd = ["env", "RUST_LOG=info", f"RAYON_NUM_THREADS={threads}", SAGE, cpath, "--write-pin", "--batch-size",
                   str(max(1, min(len(cfg["mzml_paths"]), threads))), "--disable-telemetry-i-dont-want-to-improve-sage"]
            subprocess.run(TIMED + ["--label", f"3c_sage:{ds}:{db}@t{threads}", "--threads", str(threads), "--log",
                                    os.path.join(C.RES, "timings.jsonl"), "--stdout", os.path.join(out, "sage.log"),
                                    "--meta", json.dumps({"stage": "search", "tool": "sage-0.14.7", "dataset": ds,
                                                          "db": db, "n_files": len(cfg["mzml_paths"])})]
                           + (["--lock"] if lock else []) + ["--"] + cmd, check=True, cwd=C.ROOT)
        mo = os.path.join(CACHE, "sage_mokapot", ds, db)
        if not os.path.exists(os.path.join(mo, "mokapot.peptides.txt")):
            os.makedirs(mo, exist_ok=True)
            cmd = [MOKAPOT, "--dest_dir", mo, "--seed", "11", "--max_workers", "1", "--keep_decoys", "--decoy_prefix",
                   "rev_", os.path.join(out, "results.sage.pin")]
            subprocess.run(TIMED + ["--label", f"3c_mokapot:sage:{ds}:{db}", "--threads", "1", "--log",
                                    os.path.join(C.RES, "timings.jsonl"), "--stdout", os.path.join(mo, "mokapot.log"),
                                    "--meta", json.dumps({"stage": "rescoring", "tool": "sage+mokapot", "dataset": ds,
                                                          "db": db})]
                           + (["--lock"] if lock else []) + ["--"] + cmd, check=True, cwd=C.ROOT)
        print(ds, db, "done", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd")
    ap.add_argument("ds", nargs="?")
    ap.add_argument("--dbs", nargs="+", default=DBS)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--lock-prod", action="store_true")
    a = ap.parse_args()
    if a.cmd == "mgf":
        write_mgfs()
    else:
        run(a.ds, a.dbs, a.threads, a.lock_prod)
