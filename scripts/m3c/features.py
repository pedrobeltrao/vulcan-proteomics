#!/usr/bin/env python3
"""IT-CID L3 features (label-blind; vulcan.rescore.itcid + Prosit CID/TMT-CID via Koina) for PSM directories.

    python scripts/m3c/features.py --psms DIR --family human_tmt [--index-kind prod]

For every <run>.psms.parquet in DIR writes <run>.feat.parquet (same row order) and feat_info.json.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

import pyarrow as pa  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

STORE_DIRS = [os.path.join(C.WORK, "store_train"), "/tmp/vulcan_m6/d4/store", "/tmp/vulcan_m5lc/store"]


def find_store(run):
    for d in STORE_DIRS:
        p = os.path.join(d, run + ".vspec.arrow")
        if os.path.exists(p):
            return p
    raise FileNotFoundError(run)


_caches = {}


def builder(family):
    from vulcan.predict.cache import PredCache
    from vulcan.predict.koina import KoinaProsit
    from vulcan.rescore.itcid import ItcidFeatureBuilder
    st, var = C.mods(family)
    kind = "tmt" if C.is_tmt(family) else "lf"
    pred = KoinaProsit(kind, st, var)

    def factory(nce):
        k = (kind, nce)
        if k not in _caches:
            _caches[k] = PredCache(os.path.join(C.WORK, "predcache"), pred, nce, "IT", st, var)
        return _caches[k]
    return ItcidFeatureBuilder(factory, st, var), pred


def run_dir(d, family, kind="prod", force=False):
    from vulcan.gpu.engine import load_store
    fb, pred = builder(family)
    l1 = C.load_l1(family, kind)
    infos = {"predictor": pred.version, "runs": {}}
    ip = os.path.join(d, "feat_info.json")
    if os.path.exists(ip):
        infos = json.load(open(ip))
    for f in sorted(glob.glob(os.path.join(d, "*.psms.parquet"))):
        run = os.path.basename(f)[:-len(".psms.parquet")]
        outp = os.path.join(d, run + ".feat.parquet")
        if os.path.exists(outp) and not force:
            continue
        t = pq.read_table(f)
        psm = {k: t.column(k).to_numpy() for k in ("pf_row", "charge", "spec_idx", "rank", "hyperscore", "rt")}
        t1 = time.perf_counter()
        sp = load_store(find_store(run))
        t_load = time.perf_counter() - t1
        feats, info = fb.run(l1, psm, sp)
        info["timing_s"]["load_store"] = t_load
        info["n_psms"] = int(len(psm["pf_row"]))
        pq.write_table(pa.table(feats), outp + ".partial")
        os.replace(outp + ".partial", outp)
        infos["runs"][run] = info
        json.dump(infos, open(ip, "w"), indent=1, default=float)
        print(run, json.dumps({k: info[k] for k in ("rt_mad_min", "frac_pred_missing", "cache_stats")}, default=float),
              flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--psms", required=True)
    ap.add_argument("--family", required=True)
    ap.add_argument("--index-kind", default="prod")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    run_dir(a.psms, a.family, a.index_kind, a.force)


if __name__ == "__main__":
    main()
