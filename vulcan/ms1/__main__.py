"""CLI: MS1 features and the precursor-hypothesis table for one store.

    python -m vulcan.ms1 STORE OUTDIR [--calib CALIB.json] [--device cpu|cuda] [--threads N] [--no-hypotheses]

Writes OUTDIR/<run>.features.arrow (FEATURE_SCHEMA), OUTDIR/<run>.hypotheses.arrow (HYPOTHESIS_SCHEMA) and
OUTDIR/<run>.ms1.json (timings incl. process CPU seconds, counts, hypothesis stats). Schemas: vulcan/ms1/SCHEMA.md.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m vulcan.ms1")
    ap.add_argument("store")
    ap.add_argument("outdir")
    ap.add_argument("--calib", default=None, help="vulcan.calib JSON (precursor offsets + isotope priors)")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--no-hypotheses", action="store_true")
    a = ap.parse_args(argv)
    from .. import hw
    n = hw.pin_thread_pools(a.threads)
    t_cpu0, t0 = time.process_time(), time.perf_counter()
    from ..io import store
    from .features import FeatureParams, find_features, write_features
    from .hypotheses import build_hypotheses, hypothesis_stats, write_table
    os.makedirs(a.outdir, exist_ok=True)
    run = store.run_meta(a.store)["run_name"]
    t_imp = time.perf_counter() - t0
    fr = find_features(a.store, FeatureParams(), device=a.device)
    write_features(fr.features, os.path.join(a.outdir, f"{run}.features.arrow"))
    t_feat = time.perf_counter() - t0
    c_feat = time.process_time() - t_cpu0
    out = {"run": run, "threads": n, "device": a.device, "n_features": fr.features.num_rows,
           "n_hills": fr.params["n_hills"], "timings_s": fr.timings, "import_s": t_imp,
           "features_wall_s": t_feat, "features_process_cpu_s": c_feat}
    if not a.no_hypotheses:
        calib = priors = None
        if a.calib:
            from ..calib.calibration import Calibration
            calib = Calibration.from_json(a.calib)
            priors = calib.extra.get("isotope_priors")
        t1 = time.perf_counter()
        tbl, st = build_hypotheses(a.store, fr, priors=priors, calib=calib)
        write_table(tbl, os.path.join(a.outdir, f"{run}.hypotheses.arrow"))
        out["hypotheses_wall_s"] = time.perf_counter() - t1
        out["hypothesis_stats"] = hypothesis_stats(tbl, st)
    out["total_wall_s"] = time.perf_counter() - t0
    out["total_process_cpu_s"] = time.process_time() - t_cpu0
    with open(os.path.join(a.outdir, f"{run}.ms1.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items() if k != "hypothesis_stats"}, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
