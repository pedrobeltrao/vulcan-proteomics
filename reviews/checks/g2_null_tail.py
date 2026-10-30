#!/usr/bin/env python3
"""Review check: what G2 (PSM-level win rate over ALL null spectra) does not look at -- the score TAIL.
On the G2 mixture runs (production search, top-5), run the production estimator (paired double competition) per
scorer and count, above the q<=1% / 5% threshold, null-provenance (offset + homolog-filtered wrong-species) winners
by label.  A label-neutral score on nulls gives N_null_target ~= N_null_decoy above threshold.  Also: how often each
static scorer changes the PSM winner (G2's only lever).  Dev only.
"""
import os, sys, json
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))
import m3_common as C, m3_data as MD
import gates as GT
from vulcan.rescore.static import StaticModel
from vulcan.stats._util import hash_tiebreak
from vulcan.stats.tdc import double_competition
from scipy.stats import binomtest

STATIC = ["static_linear_base", "static_gbm_base", "static_linear_pred", "static_gbm_pred"]
ms = {s: StaticModel.load(os.path.join(C.RES, "models", s + ".pkl")) for s in STATIC}
lab = MD.labels(C.index_dir("human_contam", "prod"))
acc = {s: {a: {"T_null": 0, "D_null": 0, "T_all": 0} for a in (0.01, 0.05)} for s in ["hyperscore"] + STATIC}
for run, psm, feat in MD.load_runs(C.psm_dir("G2_mix")):
    prov = np.load(os.path.join(C.WORK, "store_null", run + ".provenance.npy"))
    rows, spec = psm["pf_row"], psm["spec_idx"]
    dec = lab["is_decoy"][rows]; tseq = lab["target_seq"][rows]; pv = prov[spec]
    bad = np.zeros(prov.size, bool)
    for src in (2, 3):
        m = pv == src
        b = GT.homolog_excluded_seqs(set(tseq[m].tolist()), src)
        bad[spec[m & np.isin(tseq, list(b))]] = True
    tb = hash_tiebreak(psm["pf_id"], 11)
    key = np.array([f"{run}:{x}" for x in spec])
    pair = lab["tdc_pair"][rows]
    scores = {"hyperscore": psm["hyperscore"].astype(float)}
    scores.update({s: m.apply(psm, feat) for s, m in ms.items()})
    for s, sc in scores.items():
        r = double_competition(key, psm["pf_id"], pair, sc, dec, seed=11, tiebreak=tb)
        sp_w = spec[r.row]
        nul = (pv[r.row] > 0) & ~bad[sp_w]
        for a in (0.01, 0.05):
            thr_ok = r.q <= a
            # decoys "above threshold": decoy winners with score >= min accepted target score
            tmin = r.score[thr_ok & ~r.is_decoy].min()
            above = r.score >= tmin
            acc[s][a]["T_null"] += int((above & nul & ~r.is_decoy).sum())
            acc[s][a]["D_null"] += int((above & nul & r.is_decoy).sum())
            acc[s][a]["T_all"] += int((thr_ok & ~r.is_decoy).sum())
for s in acc:
    for a in acc[s]:
        d = acc[s][a]; n = d["T_null"] + d["D_null"]
        d["binom_p_two_sided"] = float(binomtest(d["T_null"], n, 0.5).pvalue) if n else None
print(json.dumps(acc, indent=1))
json.dump(acc, open(os.path.join(os.path.dirname(__file__), "g2_null_tail.json"), "w"), indent=1)
