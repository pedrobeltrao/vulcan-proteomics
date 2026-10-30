#!/usr/bin/env python3
"""M2b ingest profile of the direct Thermo reader (vulcan/io/thermo_direct.py + VulcanRaw .NET helper).

One fresh 1-thread process (like a cohort decode worker) converts each RAW `--reps` times, from a local NVMe copy, with
.NET section profiling on. Rep 0 of the first file carries the per-process fixed costs (runtime load, JIT, lazy Python
imports); later reps are the steady state of a persistent worker. Reports core-s per 1e5 MS2 by section.

    python scripts/m2b/ingest_profile.py --out results/.../ingest_profile.json RAW [RAW ...]
"""
from __future__ import annotations

import argparse
import json
import os
import resource
import shutil
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)


def _cs():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--stage", default=f"/tmp/vulcan_m2b_{os.getuid()}/ingest_profile")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("raws", nargs="+")
    a = ap.parse_args()
    os.environ["VULCAN_THREADS"] = "1"
    from vulcan.hw import pin_thread_pools
    pin_thread_pools(1)
    from vulcan.io.store import configure_threads
    configure_threads(1)
    os.makedirs(a.stage, exist_ok=True)
    c0, t0 = _cs(), time.perf_counter()
    from vulcan.io import thermo_direct as T
    T.load_runtime()
    res = {"runtime_load": {"core_s": _cs() - c0, "wall_s": time.perf_counter() - t0}, "runs": []}
    for raw in a.raws:
        name = os.path.basename(raw)[:-4]
        local = os.path.join(a.stage, name + ".raw")
        shutil.copyfile(raw, local)
        st = os.path.join(a.stage, name + ".vspec.arrow")
        for rep in range(a.reps):
            c, t = _cs(), time.perf_counter()
            m = T.convert(local, st, do_prefetch=False, profile=True, source_path=os.path.abspath(raw))
            c, t = _cs() - c, time.perf_counter() - t
            n2 = max(1, m["n_ms2"])
            prof = m.get("dotnet_profile_s") or {}
            r = {"run": name, "rep": rep, "instrument": m.get("instrument_model"), "n_spectra": m["n_spectra"],
                 "n_ms2": m["n_ms2"], "n_peaks": m["n_peaks"], "raw_bytes": m["source_size"], "core_s": c,
                 "wall_s": t, "core_s_per_1e5_ms2": c * 1e5 / n2, "cpu_s": m["cpu_s"],
                 "dotnet_profile_s": prof,
                 "python_side_decode_s": m["cpu_s"]["decode"] - sum(prof.values()),
                 "per_1e5_ms2": {**{k: v * 1e5 / n2 for k, v in prof.items()},
                                 "python_batch_copy": (m["cpu_s"]["decode"] - sum(prof.values())) * 1e5 / n2,
                                 "arrow_write": m["cpu_s"]["arrow_write"] * 1e5 / n2}}
            res["runs"].append(r)
            print(json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()
                              if k not in ("cpu_s", "dotnet_profile_s")}), flush=True)
        for f in (local, st, st + ".json"):
            if os.path.exists(f):
                os.remove(f)
    steady = [r for r in res["runs"] if r["rep"] > 0]
    tot = sum(r["n_ms2"] for r in steady)
    res["steady_core_s_per_1e5_ms2"] = sum(r["core_s"] for r in steady) * 1e5 / tot
    keys = steady[0]["per_1e5_ms2"].keys()
    res["steady_sections_per_1e5_ms2"] = {k: sum(r["per_1e5_ms2"][k] * r["n_ms2"] for r in steady) / tot for k in keys}
    first = res["runs"][0]
    res["first_conversion_core_s_per_1e5_ms2"] = first["core_s_per_1e5_ms2"]
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)
    print("steady", round(res["steady_core_s_per_1e5_ms2"], 2), {k: round(v, 2) for k, v in
                                                                  res["steady_sections_per_1e5_ms2"].items()})


if __name__ == "__main__":
    main()
