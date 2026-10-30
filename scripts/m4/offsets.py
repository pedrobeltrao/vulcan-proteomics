#!/usr/bin/env python3
"""Mass-offset list for D0 (DESIGN L2 mass-offset mode): Unimod deltas (vulcan.calib.pilot.KNOWN_OFFSETS, minus
isotope/label entries that the closed hypotheses or a label-free sample make irrelevant) + label-blind pilot
histogram peaks (results/2026-10-27_M2_ms1_pilot, Sage-as-searcher pilot on the 10% subsample). |delta| < 0.5 Da is
covered by closed hypotheses (isotope errors 0..2) and is excluded. Deduplicated within 0.01 Da."""
import json
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
from vulcan.calib.pilot import KNOWN_OFFSETS, NEUTRON  # noqa: E402

skip = ("isotope", "SILAC", "TMT", "iTRAQ")
uni = {k: v for k, v in KNOWN_OFFSETS.items() if not k.startswith(skip) and abs(v) >= 0.5}
m = json.load(open(os.path.join(ROOT, "results/2026-10-27_M2_ms1_pilot/metrics.json")))
pil = []
for run, d in m["pilot"].items():
    for pk in d.get("offset_histogram_top", []):
        x = pk["delta_da"]
        if abs(x) < 0.5 or (abs(x - round(x / NEUTRON) * NEUTRON) <= 0.02 and abs(round(x / NEUTRON)) <= 3):
            continue
        pil.append(x)
allv = sorted(list(uni.values()) + pil)
out = []
for v in allv:
    if not out or abs(v - out[-1]) > 0.01:
        out.append(round(v, 5))
json.dump(out, open(sys.argv[1], "w"))
print(len(uni), "unimod;", len(pil), "pilot peaks;", len(out), "offsets:", out)
