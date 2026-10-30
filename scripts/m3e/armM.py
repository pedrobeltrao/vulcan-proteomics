#!/usr/bin/env python3
"""Arm M on dev (DESIGN §6.3 (M), A1.5): D7 yeast trypsin-only (TrP) runs LM2486/2487/2490 (PXD031616; yeast is
absent from scorer and AlphaPeptDeep training, D4 leakage check) searched against yeast + contaminants + P37840 with
human entrapment subsampled to r = 1 (scripts/m3e/armM_prepare.py; seeds 2001-2003 = 3 replicates).
The 3 runs form one experiment (as D0's 2 runs).  Frozen static_gbm_pred (A1.1) and hyperscore (reference);
production estimator (paired double competition, peptidoform units); combined estimator FDP_C = N_E(1+1/r)/(N_T+N_E)
with r on the effective search space (peptide level, = 1.0); verdict per §6.3 table under both CI conventions.
    --search   (GPU)  vulcan.gpu.run-equivalent search, report_psms = 5
    --features (pred env) scripts/m3/features.py
    --eval
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m3"))
import m3_common as C  # noqa: E402

import numpy as np  # noqa: E402

W = os.environ.get("VULCAN_M3E_WORK", "/tmp/vulcan_m3e")
RUNS = ["LM2486_DDA", "LM2487_DDA", "LM2490_DDA"]
SEEDS = (2001, 2002, 2003)
ALPHAS = (0.001, 0.01, 0.05)
RES = os.path.join(C.ROOT, "results", "2026-10-27_M3e_hardening")


def idx(s):
    return os.path.join(W, "index", f"yeast_M{s}")


def pdir(s):
    return os.path.join(W, "psms", f"D7_M{s}")


def search():
    import pyarrow.parquet as pq
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    p = E.SearchParams(report_psms=C.REPORT)
    stores = [E.load_store(os.path.join(C.WORK, "store_null", r + ".vspec.arrow")) for r in RUNS]
    for s in SEEDS:
        os.makedirs(pdir(s), exist_ok=True)
        l1 = L1Index(idx(s))
        gi = E.GpuIndex(l1, p, C.STATIC, C.VAR)
        for sp in stores:
            out = os.path.join(pdir(s), sp.run + ".psms.parquet")
            if not os.path.exists(out):
                tbl, st, _ = E.search(sp, gi, l1, p)
                pq.write_table(tbl, out)
                print(s, sp.run, tbl.num_rows, flush=True)
        del gi
        import torch
        torch.cuda.empty_cache()


def features():
    py = os.path.join(C.ROOT, "conda_envs", "pred", "bin", "python")
    subprocess.run([py, os.path.join(C.ROOT, "scripts", "m3", "features.py"), "--psms", *[pdir(s) for s in SEEDS],
                    "--index", *[idx(s) for s in SEEDS], "--dataset", "D7"], check=True)


def evaluate(ci_conventions=("a", "b")):
    import m3_data as MD
    from vulcan.rescore.static import StaticModel
    from vulcan.stats import entrapment as ENT
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import double_competition
    model = StaticModel.load(os.path.join(C.RES, "models", "static_gbm_pred.pkl"))
    dbi = json.load(open(os.path.join(W, "armM", "armM_db.json")))
    r_eff = 1.0
    out = {"runs": RUNS, "r_peptide_eff": {s: dbi["seeds"][str(s)]["r_peptide_eff"] for s in SEEDS},
           "r_per_mass_decile": {s: dbi["seeds"][str(s)]["r_per_mass_decile"] for s in SEEDS},
           "r_peptidoform": {}, "scorers": {}}
    per = {}
    for s in SEEDS:
        runs = MD.load_runs(pdir(s))
        psm = MD.concat([p for _, p, _ in runs])
        lab = MD.labels(idx(s))
        nb = json.load(open(os.path.join(idx(s), "meta.json")))["n_peptidoforms_by_family"]
        out["r_peptidoform"][s] = nb["E"] / nb["T"]
        rows = psm["pf_row"]
        spec = np.array([f"{a}:{b}" for a, b in zip(psm["run"], psm["spec_idx"])])
        unit = psm["pf_id"]
        tb = hash_tiebreak(unit, 11)
        sc = {"hyperscore": psm["hyperscore"].astype(float),
              "static_gbm_pred": np.concatenate([model.apply(p, f) for _, p, f in runs])}
        for name, score in sc.items():
            r = double_competition(spec, unit, lab["tdc_pair"][rows], score, lab["is_decoy"][rows], seed=11, tiebreak=tb)
            per.setdefault(name, {})[s] = (r, lab["is_entrapment"][rows][r.row], unit[r.row])
    for name, d in per.items():
        rec = {"ids_target_side": {}, "entrapment": {}}
        for a in ALPHAS:
            reps = []
            for s in SEEDS:
                r, ent, u = d[s]
                ts = ~r.is_decoy
                acc = (r.q <= a) & ts
                reps.append(ENT.Replicate(ENT.combined_counts(u[ts], ent[ts], acc[ts]), s, 1001))
                rec["ids_target_side"].setdefault(str(a), {})[s] = {"N_T": int((acc & ~ent).sum()),
                                                                    "N_E": int((acc & ent).sum())}
            rec["entrapment"][str(a)] = {cv: ENT.evaluate_arm(reps, a, "combined", r=r_eff, ci_convention=cv)
                                         for cv in ci_conventions}
            for cv in ci_conventions:
                e = rec["entrapment"][str(a)][cv]
                print(name, a, cv, round(e["mean_fdp"], 5), [round(x, 5) for x in e["ci95"]], e["verdict"],
                      "LB", round(e["mean_lb"], 5), flush=True)
        out["scorers"][name] = rec
    os.makedirs(RES, exist_ok=True)
    json.dump(out, open(os.path.join(RES, "armM.json"), "w"), indent=1, default=float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--search", action="store_true")
    ap.add_argument("--features", action="store_true")
    ap.add_argument("--eval", action="store_true")
    a = ap.parse_args()
    if a.search:
        search()
    if a.features:
        features()
    if a.eval:
        evaluate()


if __name__ == "__main__":
    main()
