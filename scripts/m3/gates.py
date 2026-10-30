#!/usr/bin/env python3
"""M3 gates for learned / predicted components (DESIGN §6.4): G0 (label swap), G2 (differential null test), G5 (leak
AUC on G2 nulls), plus G1's diagnostic sequence-only AUC.  Dev only.

G0  The E2 swap index (/tmp/vulcan_e2/g0/index_swapped: FASTA order shuffled, pair list permuted, target/decoy roles
    swapped; D1 + E2001) vs the original index: the D1 PSM tables (top-5) are searched against both, the L3
    predicted features and every frozen static score are computed on both, and must be bit-identical.
G2  On the G2 mixture (scripts/m3/g2_nulls.py), per null spectrum: target win = the PSM-level winner (all hypotheses
    and candidates of the spectrum, hash tie-break) is a target.  fixed = hyperscore; learned = each static model
    (re-ranks the top-5) and RESET (PSM competition is never redone with the learned score, so its PSM-level
    win rate is identical to fixed by construction; reported).  vulcan.stats.g2.g2_differential, class OT-HCD.
G5  leak_auc(learned score of the null winners, target win) <= 0.55 for each learned scorer (RESET: scores of the
    null competitions from a RESET fit on the whole mixture).

    conda_envs/pred/bin/python scripts/m3/gates.py [--g0] [--g2]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402
import m3_data as MD  # noqa: E402

import numpy as np  # noqa: E402

from vulcan.rescore.static import BASE_FEATURES, PRED_FEATURES, StaticModel, feature_matrix  # noqa: E402
from vulcan.stats import g2 as G2  # noqa: E402
from vulcan.stats._util import hash_tiebreak  # noqa: E402
from vulcan.stats.reset import reset  # noqa: E402
from vulcan.stats.tdc import double_competition, psm_competition  # noqa: E402

STATIC = ["static_linear_base", "static_gbm_base", "static_linear_pred", "static_gbm_pred"]
SEED = 11


def models():
    return {s: StaticModel.load(os.path.join(C.RES, "models", s + ".pkl")) for s in STATIC}


# ------------------------------------------------------------------------------------------------ G0
def g0():
    import pyarrow.parquet as pq
    import torch
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    from vulcan.predict.apd import APD
    from vulcan.predict.cache import PredCache
    from vulcan.rescore.features import FeatureBuilder
    base = C.index_dir("pfu_contam", "E2001")
    alt = os.path.join(C.E2.WORK, "g0", "index_swapped")
    sp = E.load_store(C.E2.store_path(C.DEV["D1"]["raw"][0]))
    p = E.SearchParams(report_psms=C.REPORT)
    pred = APD()
    caches = {}
    ms = models()
    outs = []
    for idx in (base, alt):
        l1 = L1Index(idx)
        gi = E.GpuIndex(l1, p, C.STATIC, C.VAR)
        tbl, _, _ = E.search(sp, gi, l1, p)
        del gi
        torch.cuda.empty_cache()
        psm = {c: tbl.column(c).to_numpy(zero_copy_only=False) for c in MD.PSM_COLS}
        # fresh cache dir per side so that no prediction is shared between the two
        root = os.path.join(C.WORK, "g0_cache", os.path.basename(idx))

        def fac(nce, root=root):
            if (root, nce) not in caches:
                caches[(root, nce)] = PredCache(root, pred, nce, "QE", C.STATIC, C.VAR)
            return caches[(root, nce)]
        feat, info = FeatureBuilder(fac, C.STATIC, C.VAR).run(l1, psm, sp)
        scores = {s: m.apply(psm, feat) for s, m in ms.items()}
        outs.append((tbl, feat, scores, info["nce"]))
    (t0, f0, s0, n0), (t1, f1, s1, n1) = outs
    res = {"psm_tables_bit_identical": bool(t0.equals(t1)), "nce_same": n0 == n1,
           "features_bit_identical": {k: bool(np.array_equal(f0[k], f1[k])) for k in f0},
           "static_scores_bit_identical": {k: bool(np.array_equal(s0[k], s1[k])) for k in s0},
           "n_psms": t0.num_rows}
    res["pass"] = res["psm_tables_bit_identical"] and res["nce_same"] and all(res["features_bit_identical"].values()) \
        and all(res["static_scores_bit_identical"].values())
    return res


# ------------------------------------------------------------------------------------------------ G2 / G5
_FOREIGN = {}


def foreign_sets(src):
    """Foreign (I/L-collapsed) tryptic peptides + their isobaric and Hamming-1 keys: src 2 = P. furiosus + contaminants
    (from the pfu_contam L1 index), src 3 = yeast (data/fasta/yeast_UP000002311.fasta) + contaminants + the spiked
    human alpha-synuclein (P37840) of D7."""
    if src in _FOREIGN:
        return _FOREIGN[src]
    sys.path.insert(0, os.path.join(C.ROOT, "scripts", "entrapment"))
    from peptides import isobaric_key
    if src == 2:
        lab = MD.labels(C.index_dir("pfu_contam", "prod"))
        F = set(lab["target_seq"][lab["family"] == 0].tolist())
    else:
        from vulcan.index._core import core
        text = "".join(open(os.path.join(C.ROOT, "data", "fasta", f)).read()
                       for f in ("yeast_UP000002311.fasta", "contaminants_frankenfield2022.fasta"))
        hum = open(os.path.join(C.ROOT, "data", "fasta", "human_sp_canonical.fasta")).read().split("\n>")
        text += "".join(">" + e.lstrip(">") + "\n" for e in hum if "|P37840|" in e.split("\n")[0])
        F = set(core.digest_fasta(text, 2, 7, 50, "KR", "P", True)["keys"])
    fiso = {isobaric_key(f) for f in F}
    fkeys = set()
    for f in F:
        for i in range(len(f)):
            fkeys.add(hash((f[:i], f[i + 1:], len(f))))
    _FOREIGN[src] = (F, fiso, fkeys, isobaric_key)
    return _FOREIGN[src]


def homolog_excluded_seqs(cand_seqs, src):
    """Human target peptides that are I/L-identical, isobaric (N<->GG, Q<->AG) or single-substitution neighbours of
    a foreign (or contaminant) tryptic peptide."""
    F, fiso, fkeys, isobaric_key = foreign_sets(src)
    return {s for s in cand_seqs if s in F or isobaric_key(s) in fiso
            or any(hash((s[:i], s[i + 1:], len(s))) in fkeys for i in range(len(s)))}


def g2():
    ms = models()
    lab = MD.labels(C.index_dir("human_contam", "prod"))
    d = C.psm_dir("G2_mix")
    runs = MD.load_runs(d)
    rec = {"win_fixed": [], "score_fixed": [], "prov": [], "win": {s: [] for s in STATIC}, "score": {s: [] for s in STATIC},
           "reset_win": [], "reset_score": []}
    excl_counts = {}
    for run, psm, feat in runs:
        prov = np.load(os.path.join(C.WORK, "store_null", run + ".provenance.npy"))
        rows = psm["pf_row"]
        dec = lab["is_decoy"][rows]
        tseq = lab["target_seq"][rows]
        spec = psm["spec_idx"]
        tb = hash_tiebreak(psm["pf_id"], SEED)
        pv = prov[spec]
        # homolog filter for species nulls (label-symmetric: pair target sequence of ANY candidate)
        bad_spec = np.zeros(prov.size, bool)
        excl_counts[run] = {}
        for src in (2, 3):
            sp_null = pv == src
            bad = homolog_excluded_seqs(set(tseq[sp_null].tolist()), src)
            bad_spec[spec[sp_null & np.isin(tseq, list(bad))]] = True
            excl_counts[run][f"source{src}"] = {"spectra_excluded": int(np.unique(spec[sp_null & np.isin(tseq, list(bad))]).size),
                                                "homolog_peptides": len(bad)}
        null_row = (pv > 0) & ~bad_spec[spec]
        hs = psm["hyperscore"].astype(float)
        w0 = psm_competition(spec, hs, tiebreak=tb)
        w0 = w0[null_row[w0]]
        order_spec = spec[w0]
        rec["win_fixed"].append(~dec[w0])
        rec["score_fixed"].append(hs[w0])
        rec["prov"].append(prov[order_spec])
        for s, m in ms.items():
            sc = m.apply(psm, feat)
            w = psm_competition(spec, sc, tiebreak=tb)
            w = w[null_row[w]]
            # align by spectrum with the fixed winners
            o = np.argsort(spec[w])
            pos = np.searchsorted(spec[w][o], order_spec)
            ww = w[o][pos]
            assert np.array_equal(spec[ww], order_spec)
            rec["win"][s].append(~dec[ww])
            rec["score"][s].append(sc[ww])
        # RESET on the whole mixture (production: hyperscore double competition winners, LDA over base+pred)
        pair = lab["tdc_pair"][rows]
        key = np.array([f"{run}:{x}" for x in spec])
        r = double_competition(key, psm["pf_id"], pair, hs, dec, seed=SEED, tiebreak=tb)
        X = feature_matrix(psm, feat, BASE_FEATURES + PRED_FEATURES)
        rs = reset(X[r.row], r.is_decoy, init_score=r.score, seed=SEED)
        # learned score for every PSM row via the final linear map is not exposed; score the null PSM winners that
        # survived to the pair stage (subset) for G5
        nr = null_row[r.row]
        rec["reset_win"].append(~r.is_decoy[nr])
        rec["reset_score"].append(rs.score[nr])
    cat = lambda x: np.concatenate(x)  # noqa: E731
    wf, sf, prov = cat(rec["win_fixed"]), cat(rec["score_fixed"]), cat(rec["prov"])
    out = {"n_null_competitions": int(wf.size), "n_offset": int((prov == 1).sum()), "n_species_pfu": int((prov == 2).sum()),
           "n_species_yeast": int((prov == 3).sum()),
           "homolog_filter": excl_counts, "p_target_win_fixed": float(wf.mean()),
           "budget_g2_sample_size": G2.g2_sample_size(), "scorers": {}}
    for s in STATIC:
        wl, sl = cat(rec["win"][s]), cat(rec["score"][s])
        res = G2.g2_differential(wf, wl, sf)
        by_src = {src: G2.g2_differential(wf[prov == k], wl[prov == k], sf[prov == k])["pass"] for src, k in
                  (("offset", 1), ("pfu", 2), ("yeast", 3)) if (prov == k).sum() > 100}
        out["scorers"][s] = {"g2": res, "g2_pass_by_source": by_src, "g5": G2.leak_auc(sl, wl),
                             "g5_fixed_reference": G2.leak_auc(sf, wf), "g2_absolute_diag": G2.g2_absolute(wl, sl),
                             "frac_winner_changed": float(np.mean(wl != wf))}
        print(s, "G2", res["pass"], round(res["strata"][0]["overall"]["diff"], 4), round(res["strata"][0]["overall"]["ucb"], 4),
              "top", round(res["strata"][0]["top"]["ucb"], 4), "G5", out["scorers"][s]["g5"], flush=True)
    rw, rsc = cat(rec["reset_win"]), cat(rec["reset_score"])
    out["scorers"]["reset_lda_pred"] = {
        "g2": "identical to fixed by construction (RESET never redoes PSM/peptide competition, §6.4); diff = 0",
        "g5": G2.leak_auc(rsc, rw), "n_null_pair_winners": int(rw.size)}
    print("reset G5", out["scorers"]["reset_lda_pred"]["g5"])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--g0", action="store_true")
    ap.add_argument("--g2", action="store_true")
    ap.add_argument("--out", default=os.path.join(C.RES, "gates.json"))
    a = ap.parse_args()
    out = json.load(open(a.out)) if os.path.exists(a.out) else {}
    if a.g0:
        out["G0"] = g0()
        print("G0", out["G0"]["pass"])
    if a.g2:
        out["G2_G5"] = g2()
    json.dump(out, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
