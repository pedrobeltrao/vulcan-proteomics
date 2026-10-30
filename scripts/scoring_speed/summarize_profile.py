#!/usr/bin/env python3
"""Print median per-sub-stage wall/core-s over measured runs (excludes warm-up) of a profile JSON."""
import json
import sys

import numpy as np

d = json.load(open(sys.argv[1]))
runs = [r for r in d["runs"].values() if r["tag"] != "warmup"]
keys = sorted({k for r in runs for k in r["sub_stages"]})
print(f"{'stage':42s} {'med wall':>9s} {'med core':>9s} {'sum wall':>9s}")
for k in keys:
    w = np.array([r["sub_stages"].get(k, {"wall_s": 0})["wall_s"] for r in runs])
    c = np.array([r["sub_stages"].get(k, {}).get("cpu_core_s", np.nan) for r in runs])
    print(f"{k:42s} {np.median(w):9.3f} {np.median(c):9.3f} {w.sum():9.2f}")
for k in ("scoring_total", "end_to_end_total"):
    w = np.array([r[k]["wall_s"] for r in runs]); c = np.array([r[k]["cpu_core_s"] for r in runs])
    print(f"{k:42s} {np.median(w):9.3f} {np.median(c):9.3f} {w.sum():9.2f}  max {w.max():.2f}")
print({k: (round(v["wall_s"], 2) if isinstance(v, dict) else v) for k, v in d["one_time"].items()})
