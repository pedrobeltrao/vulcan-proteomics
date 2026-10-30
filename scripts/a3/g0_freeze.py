#!/usr/bin/env python3
"""A3.1.4: G0 (label blindness, DESIGN §6.4) for the FULL v2_sep pipeline at the A3 freeze commit.

As scripts/e7/gates_v2.py g0 (D1 + E2001: the E2 index vs the E2 role-swapped, row-shuffled index
/tmp/vulcan_e2/g0/index_swapped), extended to every stage of the T1 v2_sep path (bench/t1/ours_chimeric.py):
  1. MS1  `python -m vulcan.ms1` (the committed default = E7b's exact fast kernels) run SEPARATELY for each side into
          its own fresh directory; the feature and hypothesis tables must be bit-identical (file sha256 and every
          column) -- they are label-free inputs, so this checks they cannot depend on anything but the store.
  2. search  vulcan.chimeric.search.search_chimeric with each side's OWN precursor groups (from its own MS1 table),
          E2 production SearchParams, report_psms 5: chimeric PSM tables bit-identical (every column, NaN-aware).
  3. features  FeatureBuilderV2.run_v2 (the T1 path's ScorerV2 builder), each side with a FRESH prediction cache
          (AlphaPeptDeep, pinned weights): every v2 feature column bit-identical, same NCE.
  4. scores  frozen v1 (static_gbm_pred) and v2 (static_v2) applied: bit-identical.
Refuses to run unless the working tree is clean (A3.1.3); records HEAD, the frozen tree hashes and the model sha256.

    python scripts/a3/g0_freeze.py --out results/2026-10-28_A3_freeze/g0.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

STATIC, VAR = {"C": 57.021464}, {"M": 15.9949}
V1 = os.path.join(ROOT, "results", "2026-10-27_M3_scoring", "models", "static_gbm_pred.pkl")
V2 = os.path.join(ROOT, "results", "2026-10-28_E7_chimeric_v2", "models", "static_v2.pkl")
STORE = f"/tmp/vulcan_t1_{os.getuid()}/D1/stage/Velos005137.vspec.arrow"
SIDES = (("base", "/tmp/vulcan_e2/index/pfu_contam_E2001"), ("swapped", "/tmp/vulcan_e2/g0/index_swapped"))
FROZEN = ["vulcan", "bench", "scripts/a3", "scripts/e7"]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def git(*a):
    return subprocess.run(["git", "-C", ROOT, *a], capture_output=True, text=True, check=True).stdout.strip()


def arrow_cols_equal(pa_a, pa_b):
    if pa_a.column_names != pa_b.column_names or pa_a.num_rows != pa_b.num_rows:
        return {"_schema_or_rows": False}
    out = {}
    for c in pa_a.column_names:
        x, y = pa_a.column(c).to_numpy(zero_copy_only=False), pa_b.column(c).to_numpy(zero_copy_only=False)
        if x.dtype.kind == "f":
            out[c] = x.shape == y.shape and bool(np.array_equal(x.view(np.uint8), y.view(np.uint8)))
        else:
            out[c] = bool(np.array_equal(x, y))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", default=f"/tmp/vulcan_a3freeze_{os.getuid()}/g0")
    a = ap.parse_args()
    clean = git("status", "--porcelain", "--untracked-files=no") == "" and \
        git("status", "--porcelain", "--untracked-files=all", "--", "vulcan", "bench", "scripts", "config", "tests") == ""
    if not clean:
        sys.exit("g0_freeze: working tree not clean (A3.1.3)")
    head = git("rev-parse", "HEAD")
    import pyarrow as pa
    import torch

    from vulcan.chimeric import search as CS
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    from vulcan.predict.apd import APD, model_sha256
    from vulcan.predict.cache import PredCache
    from vulcan.rescore.features_v2 import FeatureBuilderV2, model_psm, read_ms1_features
    from vulcan.rescore.static import StaticModel

    shutil.rmtree(a.work, ignore_errors=True)
    os.makedirs(a.work)
    rec = {"git_sha": head, "clean_tree": clean, "frozen_tree_hashes": {p: git("rev-parse", f"HEAD:{p}") for p in FROZEN},
           "store": STORE, "store_sha256": sha(STORE), "indexes": dict(SIDES),
           "models": {"v1": sha(V1), "v2": sha(V2)}, "alphapeptdeep_weights_sha256": model_sha256(),
           "hostname": os.uname().nodename, "start": time.strftime("%Y-%m-%dT%H:%M:%S")}
    assert rec["models"]["v1"].startswith("31bfe685") and rec["models"]["v2"].startswith("3f2c69cc")
    ms = {"v1": StaticModel.load(V1), "v2": StaticModel.load(V2)}
    p = E.SearchParams(report_psms=5)
    sp = E.load_store(STORE)
    side = {}
    for name, idx in SIDES:
        t0 = time.perf_counter()
        md = os.path.join(a.work, name, "ms1")
        os.makedirs(md)
        subprocess.run([sys.executable, "-m", "vulcan.ms1", STORE, md, "--threads", "8"], check=True, cwd=ROOT,
                       env={**os.environ, "PYTHONPATH": ROOT}, stdout=subprocess.DEVNULL)
        run = os.path.basename(STORE)[:-len(".vspec.arrow")]
        hyp, feat = os.path.join(md, run + ".hypotheses.arrow"), os.path.join(md, run + ".features.arrow")
        groups = CS.precursor_groups(CS.read_hypotheses(hyp), sp)
        mf = read_ms1_features(feat)
        l1 = L1Index(idx)
        gi = E.GpuIndex(l1, p, STATIC, VAR)
        tbl, st = CS.search_chimeric(sp, gi, l1, p, groups)
        del gi
        torch.cuda.empty_cache()
        psm = {c: tbl.column(c).to_numpy(zero_copy_only=False) for c in tbl.column_names}
        pred = APD()
        caches = {}
        croot = os.path.join(a.work, name, "predcache")       # fresh prediction cache per side

        def fac(nce, caches=caches, croot=croot, pred=pred):
            if nce not in caches:
                caches[nce] = PredCache(croot, pred, nce, "QE", STATIC, VAR)
            return caches[nce]
        f, info = FeatureBuilderV2(fac, STATIC, VAR).run_v2(l1, psm, sp, mf)
        mp = model_psm(psm)
        side[name] = {"hyp": hyp, "feat": feat, "tbl": tbl, "f": f, "s": {k: m.apply(mp, f) for k, m in ms.items()},
                      "nce": info["nce"], "st": st, "wall_s": time.perf_counter() - t0}
        del pred, caches
        torch.cuda.empty_cache()
        print(name, tbl.num_rows, info["nce"], round(side[name]["wall_s"], 1), flush=True)
    b, s = side["base"], side["swapped"]

    def read_arrow(path):
        with pa.OSFile(path, "rb") as src:
            return pa.ipc.open_file(src).read_all()
    ms1 = {}
    for k in ("hyp", "feat"):
        cols = arrow_cols_equal(read_arrow(b[k]), read_arrow(s[k]))
        ms1[k] = {"file_sha256": [sha(b[k]), sha(s[k])], "file_identical": sha(b[k]) == sha(s[k]),
                  "columns_identical": all(cols.values()), "columns_not_identical": [c for c, v in cols.items() if not v],
                  "n_rows": read_arrow(b[k]).num_rows}
    pcols = arrow_cols_equal(b["tbl"], s["tbl"])
    res = {
        "ms1_hypotheses_bit_identical": ms1["hyp"]["file_identical"] and ms1["hyp"]["columns_identical"],
        "ms1_features_bit_identical": ms1["feat"]["file_identical"] and ms1["feat"]["columns_identical"],
        "ms1": ms1,
        "chimeric_psm_tables_bit_identical": bool(pcols) and all(pcols.values()),
        "psm_columns_not_identical": [c for c, v in pcols.items() if not v],
        "arrow_equals_raw": bool(b["tbl"].equals(s["tbl"])),
        "n_psms": b["tbl"].num_rows, "n_secondary_psms": int((b["tbl"].column("group").to_numpy() > 0).sum()),
        "n_secondary_groups": b["st"].get("n_secondary_groups"),
        "nce_same": b["nce"] == s["nce"], "nce": b["nce"],
        "features_bit_identical": {k: bool(np.array_equal(np.asarray(b["f"][k]).view(np.uint8),
                                                          np.asarray(s["f"][k]).view(np.uint8))) for k in b["f"]},
        "scores_bit_identical": {k: bool(np.array_equal(np.asarray(b["s"][k]).view(np.uint8),
                                                        np.asarray(s["s"][k]).view(np.uint8))) for k in b["s"]},
        "n_feature_columns": len(b["f"]), "wall_s": {k: v["wall_s"] for k, v in side.items()},
    }
    res["pass"] = bool(res["ms1_hypotheses_bit_identical"] and res["ms1_features_bit_identical"] and
                       res["chimeric_psm_tables_bit_identical"] and res["nce_same"] and
                       all(res["features_bit_identical"].values()) and all(res["scores_bit_identical"].values()))
    rec["G0"] = res
    rec["end"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump(rec, open(a.out, "w"), indent=1, default=float)
    print("G0", res["pass"], {k: res[k] for k in ("ms1_hypotheses_bit_identical", "ms1_features_bit_identical",
                                                  "chimeric_psm_tables_bit_identical", "nce_same")}, flush=True)
    sys.exit(0 if res["pass"] else 1)


if __name__ == "__main__":
    main()
