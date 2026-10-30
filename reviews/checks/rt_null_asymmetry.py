#!/usr/bin/env python3
"""Review check (M3 independent review): is the RT feature a label-correlated artifact on null spectra, and does the
G2 mixture exercise it the way production does?  Dev data only (/tmp/vulcan_m3 work dir).

(1) G2 mixture: per null provenance, RT range of null spectra vs the host D0 run, and target-vs-decoy candidate
    rt_abs_resid / pred_sa (all top-5 candidates).
(2) In-situ null proxy on D0_prod/D0_E2001: spectra whose best hyperscore is in the bottom 40% of the run (almost all
    unidentifiable), candidate feature means by family (T, E, dec(T), dec(E)), and the PSM-level target win rate
    under hyperscore vs each frozen scorer.
"""
import os, sys, json
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))
import m3_common as C, m3_data as MD
from vulcan.rescore.static import StaticModel
from vulcan.stats._util import hash_tiebreak
from vulcan.stats.tdc import psm_competition

STATIC = ["static_linear_base", "static_gbm_base", "static_linear_pred", "static_gbm_pred"]
ms = {s: StaticModel.load(os.path.join(C.RES, "models", s + ".pkl")) for s in STATIC}
out = {}

def auc(pos, neg):
    from scipy.stats import mannwhitneyu
    if len(pos) == 0 or len(neg) == 0: return None
    return float(mannwhitneyu(pos, neg).statistic / (len(pos) * len(neg)))

# (1) G2 mixture
lab = MD.labels(C.index_dir("human_contam", "prod"))
g = {}
for run, psm, feat in MD.load_runs(C.psm_dir("G2_mix")):
    prov = np.load(os.path.join(C.WORK, "store_null", run + ".provenance.npy"))
    pv = prov[psm["spec_idx"]]
    dec = lab["is_decoy"][psm["pf_row"]]
    r = {}
    for k, name in ((0, "real"), (1, "offset"), (2, "pfu"), (3, "yeast")):
        m = pv == k
        if not m.any(): continue
        rr = {"rt_min_q01_q50_q99": np.quantile(psm["rt"][m], [0.01, 0.5, 0.99]).round(2).tolist(),
              "n_cand": int(m.sum())}
        for f in ("rt_abs_resid", "pred_sa"):
            rr[f + "_mean_T"] = float(feat[f][m & ~dec].mean()); rr[f + "_mean_D"] = float(feat[f][m & dec].mean())
            rr[f + "_AUC_T_vs_D"] = auc(feat[f][m & ~dec], feat[f][m & dec])
        r[name] = rr
    g[run] = r
out["G2_mixture"] = g

# (2) in-situ null proxy
for tag, idx in (("D0_prod", "prod"), ("D0_E2001", "E2001"), ("D1_prod", "prod")):
    fam = "human_contam" if tag.startswith("D0") else "pfu_contam"
    lab = MD.labels(C.index_dir(fam, idx))
    res = {}
    for run, psm, feat in MD.load_runs(C.psm_dir(tag)):
        spec = psm["spec_idx"]
        r1 = psm["rank"] == 1
        best = np.full(spec.max() + 1, -np.inf); np.maximum.at(best, spec, psm["hyperscore"])
        thr = np.quantile(psm["hyperscore"][r1], 0.40)
        low = best[spec] <= thr
        famv = lab["family"][psm["pf_row"]]
        dec = lab["is_decoy"][psm["pf_row"]]
        rr = {"n_low_spectra": int(np.unique(spec[low]).size)}
        for f in ("rt_abs_resid", "rt_resid_z", "pred_sa", "pred_entropy"):
            rr[f] = {int(c): float(feat[f][low & (famv == c)].mean()) for c in np.unique(famv[low])}
            rr[f + "_AUC_nondecoy_vs_decoy"] = auc(feat[f][low & ~dec], feat[f][low & dec])
        rr["rt_pred_AUC_T_vs_decT_all_candidates"] = auc(feat["rt_pred"][(famv == 0)], feat["rt_pred"][(famv == 2)])
        tb = hash_tiebreak(psm["pf_id"], 11)
        sub = np.flatnonzero(low)
        def twin(score):
            w = psm_competition(spec[sub], score[sub], tiebreak=tb[sub])
            w = sub[w]
            return {"target_win": float((~dec[w]).mean()),
                    "entrap_win_among_nondecoy": float(lab["is_entrapment"][psm["pf_row"][w]][~dec[w]].mean()) if idx != "prod" else None,
                    "n": int(w.size)}
        rr["win_fixed"] = twin(psm["hyperscore"].astype(float))
        for s, m in ms.items():
            rr["win_" + s] = twin(m.apply(psm, feat))
        res[run] = rr
    out[tag] = res
print(json.dumps(out, indent=1, default=float))
json.dump(out, open(os.path.join(os.path.dirname(__file__), "rt_null_asymmetry.json"), "w"), indent=1, default=float)
