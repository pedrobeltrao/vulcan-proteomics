#!/usr/bin/env python3
"""Sage 0.14.7 open search on D0 (same peptide-level DB, same MGFs as E2): precursor_tol da [-500, 150]
(Sage: candidate mass in [obs + lo, obs + hi], i.e. delta = obs - calc in [-150, +500] Da, as ours),
isotope_errors [0, 0], report_psms 5, timed with scripts/baselines/timed.py --lock (32 threads).
Sage's `wide_window` flag is the DIA/chimeric isolation-window mode (+-2.4 Da default without isolation info), not an
open search, so the open comparison uses the wide Da tolerance, as Sage's docs recommend for open search."""
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

from vulcan.gpu import sage_bridge as SB  # noqa: E402
from vulcan.gpu.host import SearchParams  # noqa: E402

mode = sys.argv[1] if len(sys.argv) > 1 else "open"
threads = int(sys.argv[2]) if len(sys.argv) > 2 else 32
rep = int(sys.argv[3]) if len(sys.argv) > 3 else 0
p = SearchParams(report_psms=5)
out = os.path.join(C.WORK, "sage", f"D0_{mode}_t{threads}_r{rep}")
os.makedirs(out, exist_ok=True)
cfg = SB.sage_config(os.path.join(C.E2WORK, "sage", "human_contam_prod.peptides.fasta"), C.D0_MGF, out, p,
                     C.STATIC, C.VAR_CLOSED, 2)
if mode == "open":
    cfg["precursor_tol"] = {"da": [-500.0, 150.0]}
    cfg["isotope_errors"] = [0, 0]
SB.write_config(cfg, os.path.join(out, "config.json"))
subprocess.run([sys.executable, os.path.join(C.ROOT, "scripts/baselines/timed.py"), "--label", f"M4_sage_{mode}_D0_t{threads}",
                "--log", os.path.join(C.RES, "timings.jsonl"), "--lock", "--threads", str(threads),
                "--stdout", os.path.join(out, "sage.log"),
                "--meta", json.dumps({"dataset": "D0", "mode": mode, "tool": "sage-0.14.7", "threads": threads,
                                      "rep": rep, "outdir": out}), "--",
                os.path.join(C.ROOT, "tools/bin/sage"), os.path.join(out, "config.json"),
                "--disable-telemetry-i-dont-want-to-improve-sage"], check=True, cwd=C.ROOT)
