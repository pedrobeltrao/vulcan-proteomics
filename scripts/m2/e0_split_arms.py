#!/usr/bin/env python3
"""E0 follow-up: split the tolerance effect into precursor-only and fragment-only tightening (calibrated spectra).

    python scripts/m2/e0_split_arms.py --out results/2026-10-27_M2_ms1_pilot

Arms (calibrated spectra; everything else as in pilot_e0.py PROD):
  F  pilot precursor tolerance, prod fragment tolerance (+-20 ppm)
  G  prod precursor tolerance (+-10 ppm), pilot fragment tolerance
  H  pilot precursor tolerance, fragment at the k=6 core rule (mu +- 6 sd)
Results are merged into <out>/pilot_e0.json under e0.
"""
import argparse
import json
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vulcan import hw  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--work", default="/tmp/vulcan_ms1")
    a = ap.parse_args()
    hw.pin_thread_pools()
    from pilot_e0 import FASTA_DIR, PROD, RUNS, calibrated, counts
    from vulcan.calib.calibration import Calibration
    from vulcan.calib.searcher import SageSearcher, SearchConfig, load_batch
    mpath = os.path.join(a.out, "pilot_e0.json")
    allres = json.load(open(mpath))
    for key, (run, ds, fasta) in RUNS.items():
        if run not in allres:
            continue
        cal = Calibration.from_json(os.path.join(a.out, f"{run}.calib.json"))
        if cal.fragment_unit != "ppm":
            continue
        tol = cal.tolerances
        mu_sd = allres[run]["pilot"]["tolerance_info"]["fragment_after"]["core"]
        f6 = (mu_sd["mu"] - 6 * mu_sd["sd"], mu_sd["mu"] + 6 * mu_sd["sd"])
        full = calibrated(load_batch(os.path.join(a.work, "store", f"{run}.vspec.arrow"), 2), cal)
        arms = {"F_cal_pilotprec_prodfrag": dict(PROD, precursor_ppm=tuple(tol["precursor_ppm_after"])),
                "G_cal_prodprec_pilotfrag": dict(PROD, fragment_ppm=tuple(tol["fragment_ppm_after"])),
                "H_cal_pilotprec_frag6sd": dict(PROD, precursor_ppm=tuple(tol["precursor_ppm_after"]),
                                                fragment_ppm=f6)}
        e0 = allres[run]["e0"]
        base = e0["A_raw_prodtol"]
        for name, kw in arms.items():
            cfg = SearchConfig(fasta=os.path.join(FASTA_DIR, fasta), **kw)
            psms, acct = SageSearcher().search(full, cfg, os.path.join(a.work, "pilot", run, "e0", name))
            v = {"tolerances": {"precursor_ppm": list(cfg.precursor_ppm), "fragment_ppm": list(cfg.fragment_ppm)},
                 **counts(psms), "sage_wall_s": acct["wall_s"], "sage_core_s": acct["core_s"]}
            v["delta_peptides_vs_A_pct"] = 100.0 * (v["peptides_q01"] / base["peptides_q01"] - 1)
            v["delta_psms_vs_A_pct"] = 100.0 * (v["psms_q01"] / base["psms_q01"] - 1)
            e0[name] = v
            print(run, name, v, flush=True)
        with open(mpath, "w") as fh:
            json.dump(allres, fh, indent=1, default=float)


if __name__ == "__main__":
    main()
