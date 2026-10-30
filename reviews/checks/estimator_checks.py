#!/usr/bin/env python3
"""Independent review (M3, item 3: estimators).  Dev data only (D0/D1 under /tmp/vulcan_m3); no sealed access.

Checks
  A. Re-ranking confound: static_gbm_pred / static_linear_pred restricted to the hyperscore rank-1 candidate of
     each spectrum (what Sage-LDA scores) vs the top-5 re-ranked list, same estimator (paired double competition).
  B. Arm C units: r at peptidoform level (the unit counted) vs r used (peptide level, M1e), and a per-mass-decile
     combined estimator sum_d N_E,d (1 + 1/r_d) / (N_T + N_E).
  C. Entrapment CI: bootstrap-only CI (no added between-seed variance) and the resulting cwc call.

    conda_envs/pred/bin/python reviews/checks/estimator_checks.py > reviews/checks/estimator_checks.out
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))
sys.path.insert(0, ROOT)
import evaluate as EV  # noqa: E402
import m3_common as C  # noqa: E402

import numpy as np  # noqa: E402

from vulcan.rescore.static import StaticModel  # noqa: E402
from vulcan.stats import entrapment as ENT  # noqa: E402

out = {}
models = {s: StaticModel.load(os.path.join(C.RES, "models", s + ".pkl"))
          for s in ("static_gbm_pred", "static_linear_pred")}

for ds in ("D0", "D1"):
    res = out.setdefault(ds, {})
    prod = EV.Search(ds, "prod")
    rank1 = prod.psm["rank"] == 1 if prod.psm["rank"].min() == 1 else prod.psm["rank"] == 0
    for sc, m in models.items():
        s = prod.static_score(m)
        full = prod.compete(s)
        s1 = np.where(rank1, s, -np.inf)
        r1 = prod.compete(s1)
        hs = prod.psm["hyperscore"].astype(float)
        fr = full.row[(full.q <= 0.01) & ~full.is_decoy]
        res[sc] = {"ids_1pct_top5_rerank": int(((full.q <= 0.01) & ~full.is_decoy).sum()),
                   "ids_1pct_rank1_only": int(((r1.q <= 0.01) & ~r1.is_decoy & np.isfinite(r1.score)).sum()),
                   "ids_0.1pct_top5_rerank": int(((full.q <= 0.001) & ~full.is_decoy).sum()),
                   "ids_0.1pct_rank1_only": int(((r1.q <= 0.001) & ~r1.is_decoy & np.isfinite(r1.score)).sum()),
                   "frac_accepted_winner_not_hyperscore_rank1": float((~rank1[fr]).mean())}
    del prod
    # B. arm C units / per-decile r
    cs = EV.Search(ds, "C")
    lab_meta = json.load(open(os.path.join(cs.index, "meta.json")))["n_peptidoforms_by_family"]
    cj = json.load(open(os.path.join(C.M1E_DB, f"{EV.FAM[ds]}.C.json")))
    r_pep = cj["r"]
    r_pf = lab_meta["E"] / lab_meta["T"]
    res["armC_r_used_peptide"] = r_pep
    res["armC_r_peptidoform"] = r_pf
    res["armC_r_per_mass_decile"] = cj["r_per_mass_decile"]
    edges = np.asarray(cj["decile_mass_edges_da"])
    from vulcan.index.build import L1Index
    l1 = L1Index(cs.index)
    for sc in ("hyperscore", "static_gbm_pred"):
        w = EV.scored(cs, sc, models if sc != "hyperscore" else {})
        t = ~w["is_decoy"]
        acc = t & (w["q"] <= 0.01)
        rows = w["row"][acc]
        ent = cs.is_ent[rows]
        mass = np.asarray(l1.mass32)[cs.psm["pf_row"][rows]] if hasattr(l1, "mass32") else None
        n_t, n_e = int((~ent).sum()), int(ent.sum())
        d = {"N_T": n_t, "N_E": n_e, "FDP_C_global_r_pep": n_e * (1 + 1 / r_pep) / (n_t + n_e),
             "FDP_C_global_r_pf": n_e * (1 + 1 / r_pf) / (n_t + n_e)}
        if mass is not None:
            dec = np.searchsorted(edges, mass)
            rd = np.asarray(cj["r_per_mass_decile"])
            ne_d = np.bincount(dec[ent], minlength=10)
            d["N_E_per_decile"] = ne_d.tolist()
            d["FDP_C_per_decile_r_pep"] = float((ne_d * (1 + 1 / rd)).sum() / (n_t + n_e))
        res[f"armC_{sc}_1pct"] = d
    del cs

# C. bootstrap-only CI from eval.json
ev = json.load(open(os.path.join(C.RES, "eval.json")))
cc = {}
for ds, d in ev.items():
    for sc, r in d.items():
        for a in ("0.01", "0.05"):
            P = r["entrapment"][a]["P"]
            f = np.asarray(P["fdp_per_replicate"])
            bw = f.var(ddof=1) / f.size
            sb = float(np.sqrt(max(P["se"] ** 2 - bw, 0)))
            lo = P["mean_fdp"] - 1.959964 * sb
            cc[f"{ds}/{sc}/{a}"] = {"mean": P["mean_fdp"], "lo_reported": P["ci95"][0], "lo_bootstrap_only": lo,
                                    "cwc_reported": P["consistent_with_control"], "cwc_bootstrap_only": bool(lo <= float(a))}
out["armP_ci_decomposition"] = cc
print(json.dumps(out, indent=1, default=float))
