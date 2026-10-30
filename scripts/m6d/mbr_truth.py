#!/usr/bin/env python3
"""M6d part 3: fresh MBR false-transfer truth check with the M6 MBR settings UNCHANGED (MBRParams() defaults: donors at
ID q <= 0.1%, claimed-feature exclusion, per-acceptor-run transfer q; fixed a-priori transfer score).  The experiment is
every run of the set (donors + acceptors, IDs pooled over all runs by scripts/m6/ids.py); acceptors are runs known to
lack a proteome, so any transfer of a precursor unique to that proteome into them is false.  Same truth definition and
estimators as M6 (lfq_d9.py, `mbr_truth_human_only`).
    scripts/m2b/vpy.sh scripts/m6d/mbr_truth.py DATASET OUT.json
"""
from __future__ import annotations

import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import datasets as DS  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from vulcan.quant import lfq as Q  # noqa: E402


def main():
    name, out = sys.argv[1], sys.argv[2]
    cfg = DS.get(name)
    X = cfg.extra
    runs = X["donor_runs"] + X["acceptor_runs"]
    acc_runs = X["acceptor_runs"]
    absent = set(X["absent_species"])           # species absent from the acceptor runs
    present = set(X["present_species"])         # species present in the acceptor runs
    ids_dir = X["ids_all"]
    ids = pd.read_parquet(os.path.join(ids_dir, "ids.parquet"))
    pep = pd.read_parquet(os.path.join(ids_dir, "peptides.parquet")).set_index("pf_id")
    acc = ids[ids["pass_0.01"]]
    t0 = time.perf_counter()
    ms1, tabs = {}, {}
    for run in runs:
        r = Q.load_run_ms1(cfg.ms1, cfg.store(run), run)
        ms1[run] = r
        tabs[run] = Q.precursor_table(Q.link_psms(acc[acc["run"] == run], r), r)
    dec = pep.index[pep["is_decoy"]].to_numpy()
    al = Q.Alignment({r: t[~t["pf_id"].isin(dec)] for r, t in tabs.items()})
    dd = ids[ids["pass_0.001"]]
    dk = set(zip(dd["pf_id"], dd["charge"].astype(np.int8)))
    P = Q.MBRParams()
    tr, info = Q.match_between_runs(tabs, ms1, al, dec, donor_keys=dk, params=P)
    wall = time.perf_counter() - t0
    tgt = ~tr["is_decoy_donor"].to_numpy()
    sp = pep["species"].reindex(tr["pf_id"].to_numpy()).to_numpy()
    cont = pep["contaminant"].reindex(tr["pf_id"].to_numpy()).to_numpy().astype(bool)
    inA = tr["run"].isin(acc_runs).to_numpy() & tgt & ~cont
    ab = np.isin(sp, list(absent))
    pr = np.isin(sp, list(present))
    res = {"dataset": name, "donor_runs": X["donor_runs"], "acceptor_runs": acc_runs, "absent_species": sorted(absent),
           "present_species": sorted(present), "mbr_params": P.__dict__, "donor_alpha": 0.001, "id_alpha": 0.01,
           "alignment": al.summary(), "mbr_wall_s_incl_link": wall, "candidates_in_acceptors": int(inA.sum()),
           "absent_candidates_in_acceptors": int((inA & ab).sum()), "present_candidates_in_acceptors": int((inA & pr).sum())}
    St, Sr, Sm = tr["St"].to_numpy(), tr["Sd_rt"].to_numpy(), tr["Sd_mz"].to_numpy()
    for tq in (0.01, 0.05):
        q = (tr["transfer_q"] <= tq).to_numpy()
        n_acc = int((inA & q & (ab | pr)).sum())
        n_ab = int((inA & q & ab).sum())
        n_ab_c, n_pr_c = int((inA & ab).sum()), int((inA & pr).sum())
        thr = np.nanmin(St[q]) if q.any() else np.inf
        rec = {"accepted_in_acceptors": n_acc, "accepted_absent_species (false)": n_ab,
               "false_fraction_lower_bound": n_ab / max(n_acc, 1),
               "combined_est_FDP": n_ab * (1 + n_pr_c / max(n_ab_c, 1)) / max(n_acc, 1),
               "absent_null_acceptance_rate": n_ab / max(n_ab_c, 1),
               "present_acceptance_rate": int((inA & q & pr).sum()) / max(n_pr_c, 1),
               "decoy_wins_in_absent_null_subset_rt_mz": [
                   int((inA & ab & np.isfinite(Sr) & (Sr >= St) & (Sr >= thr)).sum()),
                   int((inA & ab & np.isfinite(Sm) & (Sm >= St) & (Sm >= thr)).sum())],
               "absent_target_wins_above_thr": int((inA & ab & np.isfinite(St) & (St > Sr) & (St > Sm) & (St >= thr)).sum()),
               "per_acceptor_run": {}}
        for r in acc_runs:
            m = tr["run"].to_numpy() == r
            na = int((m & inA & q & (ab | pr)).sum())
            nf = int((m & inA & q & ab).sum())
            rec["per_acceptor_run"][r] = {"accepted": na, "false": nf, "lb": nf / max(na, 1)}
        # transfers into the donor runs (where both proteomes are present) for scale
        rec["accepted_in_donor_runs"] = int((~tr["run"].isin(acc_runs).to_numpy() & tgt & q).sum())
        res[f"q{tq}"] = rec
    # MS2-level false IDs in the acceptor runs (identification, not MBR) for context
    hi = acc[acc["run"].isin(acc_runs) & ~acc["is_decoy"]]
    hsp = pep["species"].reindex(hi["pf_id"].to_numpy()).to_numpy()
    hct = pep["contaminant"].reindex(hi["pf_id"].to_numpy()).to_numpy().astype(bool)
    res["ms2_psm_winners_in_acceptors"] = {"n": int(len(hi)), "absent_species": int((np.isin(hsp, list(absent)) & ~hct).sum())}
    json.dump(res, open(out, "w"), indent=1, default=float)
    print(json.dumps({k: res[k] for k in ("candidates_in_acceptors", "absent_candidates_in_acceptors")}),
          json.dumps({k: {kk: vv for kk, vv in res[k].items() if kk != "per_acceptor_run"} for k in ("q0.01", "q0.05")}))


if __name__ == "__main__":
    main()
