#!/usr/bin/env python3
"""M2: MS1 features -> hypothesis table -> L0.5 pilot (Sage interim searcher) -> calibrated hypotheses -> E0 first read.

    python scripts/m2/pilot_e0.py --out results/2026-10-27_M2_ms1_pilot [--runs 01A,01B,Velos] [--work /tmp/vulcan_ms1]

E0 arms (Sage 0.14.7 on MGFs written from the same store; full MS2 of each run; prod search settings otherwise):
  A  raw spectra,        prod tolerances (precursor +-10 ppm, fragment +-20 ppm)      = baseline
  B  calibrated spectra, pilot tolerances after recalibration                         = pilot pass on
  C  calibrated spectra, prod tolerances                                              (calibration only)
  D  raw spectra,        pilot tolerances before recalibration                        (tolerance only)
  E  calibrated spectra, background-subtracted 99%-coverage tolerances (alternative tolerance rule)
Pilot tolerances = mu +- 4 sd of the sigma-clipped anchor-error core (vulcan.calib.calibration.core_tolerance).
Counted at Sage peptide q <= 1% and PSM (spectrum) q <= 1%, targets only. Labels are used here for evaluation (L4)
only; nothing upstream (features, hypotheses, calibration curves) reads them, except the pilot's anchor threshold.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
from vulcan import hw  # noqa: E402

RUNS = {
    "01A": ("b1906_293T_proteinID_01A_QE3_122212", "D0", "human_contam.fasta"),
    "01B": ("b1937_293T_proteinID_01B_QE3_122212", "D0", "human_contam.fasta"),
    "Velos": ("Velos005137", "D1", "pfu_contam.fasta"),
}
FASTA_DIR = os.path.join(ROOT, "data", "cache", "M1e", "entrapment", "inputs")
PROD = dict(precursor_ppm=(-10.0, 10.0), fragment_ppm=(-20.0, 20.0), isotope_errors=(0, 2), precursor_charge=(2, 4),
            variable_mods={"M": [15.9949]}, max_variable_mods=2)


def counts(psms):
    import pyarrow.compute as pc
    t = psms.filter(pc.invert(psms.column("is_decoy")))
    q = t.column("q_value").to_numpy()
    pq = t.column("peptide_q").to_numpy()
    pep = np.asarray(t.column("peptide").to_pylist(), dtype=object)
    return {"psms_q01": int((q <= 0.01).sum()), "peptides_q01": int(len(set(pep[pq <= 0.01].tolist()))),
            "n_rank1": int(psms.num_rows)}


def calibrated(batch, cal):
    from dataclasses import replace
    lens = np.diff(batch.offsets)
    prt = np.repeat(batch.rt, lens)
    return replace(batch, precursor_mz=cal.correct_precursor(batch.precursor_mz, batch.rt),
                   mz=cal.correct_fragment(np.asarray(batch.mz, np.float64), prt))


def priors_from_pilot(by_source: dict, isotope_errors=(-1, 0, 1, 2, 3), pseudo=0.5) -> dict:
    out = {}
    for src in ("ms1_feature", "instrument", "fallback"):
        h = by_source.get(src, {}) or {}
        c = np.array([h.get(str(e), 0) for e in isotope_errors], float) + pseudo
        out[src] = {e: float(v) for e, v in zip(isotope_errors, c / c.sum())}
    return out


def hyp_coverage(hyp, psms, tol_ppm=10.0) -> dict:
    """Confident target PSMs of a full-run search whose (charge, peptide mass) is a row of the hypothesis table."""
    import pyarrow.compute as pc
    from vulcan.calib.pilot import hypothesis_priors
    t = psms.filter(pc.and_(pc.invert(psms.column("is_decoy")), pc.less_equal(psms.column("q_value"), 0.01)))
    r = hypothesis_priors(hyp, t, tol_ppm=tol_ppm)
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--runs", default="01A,01B,Velos")
    ap.add_argument("--work", default="/tmp/vulcan_ms1")
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--skip-arms", default="")
    a = ap.parse_args()
    nthr = hw.pin_thread_pools(a.threads)
    import pyarrow as pa  # noqa: F401
    from vulcan.calib.pilot import PilotParams, run_pilot
    from vulcan.calib.searcher import SageSearcher, SearchConfig, load_batch
    from vulcan.ms1.features import find_features
    from vulcan.ms1.hypotheses import build_hypotheses, hypothesis_stats, write_table

    os.makedirs(a.out, exist_ok=True)
    mpath = os.path.join(a.out, "pilot_e0.json")
    allres = json.load(open(mpath)) if os.path.exists(mpath) else {}
    for key in a.runs.split(","):
        run, ds, fasta = RUNS[key]
        fasta = os.path.join(FASTA_DIR, fasta)
        st = os.path.join(a.work, "store", f"{run}.vspec.arrow")
        wd = os.path.join(a.work, "pilot", run)
        os.makedirs(wd, exist_ok=True)
        res = {"run": run, "dataset": ds, "fasta": os.path.basename(fasta), "threads": nthr}
        t0 = time.perf_counter()
        c0 = time.process_time()
        fr = find_features(st)
        res["features"] = {"n": fr.features.num_rows, "wall_s": time.perf_counter() - t0,
                           "process_cpu_s": time.process_time() - c0}
        hyp0, st0 = build_hypotheses(st, fr)
        # pilot pass (seeded subsample, wide closed search, own decoys only for the anchor threshold)
        pr = run_pilot(st, fasta, os.path.join(wd), SageSearcher(annotate=True), PilotParams(), hypotheses=hyp0)
        cal = pr.calibration
        pri = priors_from_pilot(pr.report["priors"].get("by_source", {}))
        cal.extra["isotope_priors"] = {k: {str(e): v for e, v in d.items()} for k, d in pri.items()}
        cal.to_json(os.path.join(a.out, f"{run}.calib.json"))
        res["pilot"] = pr.report
        res["pilot_timings"] = pr.timings
        # calibrated hypothesis table with pilot isotope priors
        t1 = time.perf_counter()
        hyp, sth = build_hypotheses(st, fr, priors=pri, calib=cal)
        res["hypotheses_wall_s"] = time.perf_counter() - t1
        write_table(hyp, os.path.join(wd, f"{run}.hypotheses.arrow"))
        res["hypothesis_stats"] = hypothesis_stats(hyp, sth)
        # E0 arms
        full = load_batch(st, 2)
        fullc = calibrated(full, cal)
        tol = cal.tolerances
        fkey = "fragment_ppm_after" if cal.fragment_unit == "ppm" else "fragment_Da_after"
        fkey_b = fkey.replace("after", "before")
        arms = {
            "A_raw_prodtol": (full, dict(PROD)),
            "B_cal_pilottol": (fullc, dict(PROD, precursor_ppm=tuple(tol["precursor_ppm_after"]),
                                           fragment_ppm=tuple(tol[fkey]))),
            "C_cal_prodtol": (fullc, dict(PROD)),
            "D_raw_pilottol": (full, dict(PROD, precursor_ppm=tuple(tol["precursor_ppm_before"]),
                                          fragment_ppm=tuple(tol[fkey_b]))),
            "E_cal_cov99tol": (fullc, dict(PROD, precursor_ppm=tuple(tol["precursor_ppm_after_cov99"]),
                                           fragment_ppm=tuple(tol[fkey + "_cov99"]))),
        }
        if cal.fragment_unit != "ppm":   # low-res MS2: keep the prod fragment setting (ppm calibration not applied)
            for k in ("B_cal_pilottol", "D_raw_pilottol", "E_cal_cov99tol"):
                arms[k][1]["fragment_ppm"] = PROD["fragment_ppm"]
        res["e0"] = {}
        for name, (b, kw) in arms.items():
            if name in a.skip_arms.split(","):
                continue
            cfg = SearchConfig(fasta=fasta, **kw)
            psms, acct = SageSearcher().search(b, cfg, os.path.join(wd, "e0", name))
            res["e0"][name] = {"tolerances": {"precursor_ppm": list(cfg.precursor_ppm),
                                              "fragment_ppm": list(cfg.fragment_ppm)}, **counts(psms),
                               "sage_wall_s": acct["wall_s"], "sage_core_s": acct["core_s"]}
            if name == "A_raw_prodtol":
                res["hypothesis_coverage_of_baseline_ids"] = hyp_coverage(hyp, psms)
            print(run, name, res["e0"][name], flush=True)
        base = res["e0"].get("A_raw_prodtol")
        if base:
            for k, v in res["e0"].items():
                v["delta_peptides_vs_A_pct"] = 100.0 * (v["peptides_q01"] / base["peptides_q01"] - 1)
                v["delta_psms_vs_A_pct"] = 100.0 * (v["psms_q01"] / base["psms_q01"] - 1)
        res["total_wall_s"] = time.perf_counter() - t0
        allres[run] = res
        with open(mpath, "w") as fh:
            json.dump(allres, fh, indent=1, default=float)
        print(json.dumps({k: res[k] for k in ("run", "features", "e0")}, indent=1, default=float), flush=True)


if __name__ == "__main__":
    main()
