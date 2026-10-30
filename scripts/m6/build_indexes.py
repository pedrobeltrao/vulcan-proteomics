#!/usr/bin/env python3
"""Build the M6 L1 indexes (prod = targets + paired decoys, no entrapment):
  D9: ProteoBench HYE FASTA (incl. its Cont_ contaminants), CAM static, Met-ox variable (default params).
  D4: human + contaminants (same FASTA as M1e/E2), CAM + TMT-K static.  The TMT peptide-N-term static mod is not
      expressible in the index (vulcan_core: residue statics only); scripts/m6/tmt.py adds it at load time
      (mass32 + TMT for every row, mask bit 63 + nterm delta in the kernels)."""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402
from vulcan.index.build import build  # noqa: E402

which = sys.argv[1:] or ["D9", "D4"]
for w in which:
    out, fa, par = ((C.D9_INDEX, C.D9_FASTA, {}) if w == "D9" else
                    (C.D4_INDEX, C.D4_FASTA, {"static_mods": {"C": 57.021464, "K": C.TMT}}))
    if os.path.exists(os.path.join(out, "meta.json")):
        continue
    t0 = time.perf_counter()
    info = build([fa], out, par)
    print(w, round(time.perf_counter() - t0, 1), "s", json.dumps(info, default=str)[:400], flush=True)
