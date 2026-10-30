#!/usr/bin/env python3
"""E7b MS1 equality record (dev only): the fast exact MS1 path vs (a) the stored E7 outputs and (b) the reference
kernels in the same process.

For every run of D0, D1, D8h (<work>/stage/<run>.vspec.arrow):
  (a) vulcan.ms1 features + hypotheses with the current default code vs the files E7 wrote with the old code
      (<work>/ms1/<run>.features.arrow / .hypotheses.arrow): every column bitwise (values + validity), schema
      metadata equal;
  (b) Hills arrays (detect_hills reference=True vs fused) and group_isotopes (exact_reference=True vs fast) bitwise.
Writes --out JSON.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402


def read_arrow(path):
    import pyarrow as pa
    import pyarrow.ipc as ipc
    with pa.OSFile(path, "rb") as f:
        return ipc.open_file(f).read_all()


def col_equal(a, b) -> bool:
    import pyarrow.compute as pc
    if a.type != b.type or len(a) != len(b):
        return False
    if len(a) == 0:
        return True
    if not pc.all(pc.equal(pc.is_null(a), pc.is_null(b))).as_py():
        return False
    va = a.combine_chunks() if hasattr(a, "combine_chunks") else a
    vb = b.combine_chunks() if hasattr(b, "combine_chunks") else b
    try:
        na = va.to_numpy(zero_copy_only=False)
        nb = vb.to_numpy(zero_copy_only=False)
    except Exception:  # noqa: BLE001
        return va.equals(vb)
    if na.dtype.kind in "fiub":
        return na.tobytes() == nb.tobytes()
    return bool(np.array_equal(na, nb))


def table_equal(t, r) -> dict:
    cols = {c: (c in r.column_names) and col_equal(t.column(c), r.column(c)) for c in t.column_names}
    return {"columns_equal": all(cols.values()) and t.column_names == r.column_names,
            "metadata_equal": t.schema.metadata == r.schema.metadata, "n_rows": [t.num_rows, r.num_rows],
            "differing_columns": [c for c, v in cols.items() if not v]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--datasets", nargs="+", default=["D0", "D1", "D8h"])
    ap.add_argument("--work-root", default=f"/tmp/vulcan_t1_{os.getuid()}")
    a = ap.parse_args()
    from vulcan import hw
    hw.pin_thread_pools(8)
    from vulcan.ms1 import hills as HL
    from vulcan.ms1 import isotopes as IS
    from vulcan.ms1.features import FeatureParams, find_features, load_ms1
    from vulcan.ms1.hypotheses import build_hypotheses
    out = {"hw": hw.fingerprint(), "loadavg_start": os.getloadavg(), "runs": {}}
    for ds in a.datasets:
        W = os.path.join(a.work_root, ds)
        for f in sorted(os.listdir(os.path.join(W, "stage"))):
            if not f.endswith(".vspec.arrow"):
                continue
            run = f[:-len(".vspec.arrow")]
            store = os.path.join(W, "stage", f)
            t0 = time.time()
            fr = find_features(store, FeatureParams())
            hyp, _ = build_hypotheses(store, fr)
            rec = {"dataset": ds,
                   "features_vs_E7": table_equal(fr.features, read_arrow(os.path.join(W, "ms1", run + ".features.arrow"))),
                   "hypotheses_vs_E7": table_equal(hyp, read_arrow(os.path.join(W, "ms1", run + ".hypotheses.arrow")))}
            ms1 = load_ms1(store)
            p = FeatureParams()
            lens = np.diff(ms1.off)
            scan = np.repeat(np.arange(len(lens), dtype=np.int32), lens)
            kw = dict(htol_ppm=p.htol_ppm, minlh=p.minlh, hvf=p.hvf, gap=p.gap, min_intensity=p.min_intensity)
            hf = HL.detect_hills(scan, ms1.mz, ms1.intensity, ms1.ms1_rt, **kw)
            hr = HL.detect_hills(scan, ms1.mz, ms1.intensity, ms1.ms1_rt, reference=True, **kw)
            rec["hills_fast_vs_reference"] = all(
                getattr(hf, k).dtype == getattr(hr, k).dtype and getattr(hf, k).tobytes() == getattr(hr, k).tobytes()
                for k in ("off", "p_scan", "p_mz", "p_int", "p_src", "h_mz", "h_start", "h_end", "h_apex_scan",
                          "h_apex_int", "h_area", "h_len"))
            gf = IS.group_isotopes(hf)
            gr = IS.group_isotopes(hr, exact_reference=True)
            rec["isotopes_fast_vs_reference"] = bool(gf["fast_path"]) and all(
                gf[k].dtype == gr[k].dtype and gf[k].tobytes() == gr[k].tobytes()
                for k in gf if isinstance(gf[k], np.ndarray)) and gf["n_candidates"] == gr["n_candidates"]
            rec["n_hills"] = int(hf.n)
            rec["n_features"] = fr.features.num_rows
            rec["wall_s"] = time.time() - t0
            rec["ok"] = (rec["features_vs_E7"]["columns_equal"] and rec["features_vs_E7"]["metadata_equal"]
                         and rec["hypotheses_vs_E7"]["columns_equal"] and rec["hypotheses_vs_E7"]["metadata_equal"]
                         and rec["hills_fast_vs_reference"] and rec["isotopes_fast_vs_reference"])
            out["runs"][run] = rec
            print(ds, run[:30], rec["ok"], rec["features_vs_E7"]["differing_columns"],
                  rec["hypotheses_vs_E7"]["differing_columns"], rec["hills_fast_vs_reference"],
                  rec["isotopes_fast_vs_reference"], flush=True)
    out["all_ok"] = all(r["ok"] for r in out["runs"].values())
    out["n_runs"] = len(out["runs"])
    out["loadavg_end"] = os.getloadavg()
    json.dump(out, open(a.out, "w"), indent=1, default=float)
    print("ALL_OK", out["all_ok"], out["n_runs"])


if __name__ == "__main__":
    main()
