#!/usr/bin/env python3
"""Amendment A1 gates on the frozen default scorer (static_gbm_pred, A1.1), with and without the RT features.  Dev only.

Pipelines compared on the same G2 null mixtures (scripts/m3e/g2mix.py; production search, prod index, top-5 PSMs,
L3 features from scripts/m3/features.py):
  fixed       hyperscore (the fixed-score pipeline)
  gbm         static_gbm_pred, unchanged (production; RT features from the RT-mapped nulls)
  gbm_noRT    the same frozen model with the two RT features (rt_abs_resid, rt_resid_z) replaced, for every PSM of a
              run, by that run's median over all PSMs (label-blind; the RT channel carries no information).  This is
              an ablation of the frozen model, NOT a new scorer (A1.3: "reported with and without the RT features").
Per mixture run and pipeline: PSM-level winner of every spectrum (all hypotheses x top-5 candidates, hash tie-break
seed 11) and the production estimator (paired double competition over the run, peptidoform units) -> the
pipeline's own threshold tau_alpha = lowest accepted target score at q <= alpha.
Null spectra: offsets (prov 1, no filter, as specified) + species nulls (prov >= 2) that pass the A1.3 homolog filter:
  a species-null spectrum is dropped if, for ANY of the three pipelines, the matched peptide or its pair partner
  (target or decoy sequence; label-symmetric) shares a 6-mer with / is within Hamming <= 2 of its source proteome or
  the contaminants (I/L collapsed; vulcan.gates.homolog).  One common clean set for all pipelines (paired tests).
Gates (vulcan.gates.stats, vulcan.stats.g2):
  G2 (original, unchanged): g2_differential over clean nulls with a candidate (margins 0.01 / 0.02 top decile)
  tail-G2 (A1.2): at the 1% threshold (gated); 5% reported
  G5: leak AUC of the learned score over clean null winners (<= 0.55), CI and top decile
  G1: sequence-only classifier (scripts/m3e/seqclf.py, trained on training-run candidates) logit of the winning
      candidate vs learned score, top score decile of clean null winners: |rho| 99% CI inside (-0.05, 0.05)
    conda_envs/pred/bin/python scripts/m3e/gates_a1.py --out results/.../gates_a1.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m3"))
sys.path.insert(0, HERE)
import m3_common as C  # noqa: E402
import m3_data as MD  # noqa: E402

import numpy as np  # noqa: E402

from vulcan.gates import stats as GS  # noqa: E402
from vulcan.gates.homolog import HomologFilter, read_fasta_seqs  # noqa: E402
from vulcan.rescore.static import StaticModel  # noqa: E402
from vulcan.stats import g2 as G2  # noqa: E402
from vulcan.stats._util import hash_tiebreak  # noqa: E402
from vulcan.stats.tdc import double_competition, psm_competition  # noqa: E402

W = os.environ.get("VULCAN_M3E_WORK", "/tmp/vulcan_m3e")
MIX = os.path.join(W, "null_mix")
PSM = os.path.join(W, "psms", "G2mix")
F = os.path.join(C.ROOT, "data", "fasta")
CONTAM = os.path.join(F, "contaminants_frankenfield2022.fasta")
SRC_FASTA = {2: ["pfu_UP000001013.fasta"], 3: ["yeast_UP000002311.fasta"], 4: ["dn_UP000008524.fasta"],
             5: ["dn_UP000001939.fasta", "dn_UP000008816.fasta"], 6: ["dn_UP000001688.fasta"], 7: ["dn_UP000001417.fasta"]}
PROV_NAME = {1: "offset", 2: "pfu", 3: "yeast", 4: "tbrucei", 5: "saureus", 6: "synechococcus", 7: "bcereus"}
RT_FEATS = ("rt_abs_resid", "rt_resid_z")
PIPES = ("fixed", "gbm", "gbm_noRT")
SEED = 11


def neutral_rt(feat):
    f = dict(feat)
    for k in RT_FEATS:
        f[k] = np.full_like(feat[k], np.median(feat[k]))
    return f


def homolog_filters():
    contam = read_fasta_seqs([CONTAM])
    snca = [s for h, s in zip(*_fasta_with_headers(os.path.join(F, "human_sp_canonical.fasta"))) if "|P37840|" in h]
    out = {}
    for pv, fs in SRC_FASTA.items():
        extra = contam + (snca if pv == 3 else [])
        out[pv] = HomologFilter(read_fasta_seqs([os.path.join(F, f) for f in fs]) + extra)
    return out


def _fasta_with_headers(path):
    hs, ss, cur = [], [], []
    for line in open(path):
        if line.startswith(">"):
            if cur:
                ss.append("".join(cur))
            hs.append(line)
            cur = []
        else:
            cur.append(line.strip())
    ss.append("".join(cur))
    return hs, ss


def _load(d, names):
    import pyarrow.parquet as pq
    out = []
    for run in names:
        t = pq.read_table(os.path.join(d, run + ".psms.parquet"))
        psm = {c: t.column(c).to_numpy(zero_copy_only=False) for c in MD.PSM_COLS}
        ft = pq.read_table(os.path.join(d, run + ".feat.parquet"))
        out.append((run, psm, {c: ft.column(c).to_numpy() for c in ft.column_names}))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-runs", type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()
    from vulcan.index.build import L1Index
    model = StaticModel.load(os.path.join(C.RES, "models", "static_gbm_pred.pkl"))
    idx = C.index_dir("human_contam", "prod")
    lab = MD.labels(idx)
    l1 = L1Index(idx)
    import seqclf
    clf = seqclf.load()
    hfs = homolog_filters()
    print("filters ready", round(time.time() - t0), flush=True)
    seq_cache = {}

    def seqs(rows):
        out = []
        for r in rows:
            r = int(r)
            if r not in seq_cache:
                seq_cache[r] = l1.sequence(r)
            out.append(seq_cache[r])
        return out

    # pass 1: winners, thresholds, and the candidate rows to test per null source
    rec = {k: [] for k in ("prov", "srcrun", "run", "pairrows")}
    per = {p: {k: [] for k in ("win_t", "score", "tail1", "tail5", "logit")} for p in PIPES}
    thr = {}
    import glob as _g
    have = sorted(os.path.basename(f)[:-len(".feat.parquet")] for f in _g.glob(os.path.join(PSM, "*.feat.parquet")))
    allr = sorted(os.path.basename(f)[:-len(".psms.parquet")] for f in _g.glob(os.path.join(PSM, "*.psms.parquet")))
    if not a.max_runs and have != allr:
        raise SystemExit(f"features missing for {len(allr) - len(have)} mixtures")
    runs = [x for x in MD.load_runs_subset(PSM, have[:a.max_runs] if a.max_runs else have)] if hasattr(MD, "load_runs_subset") \
        else _load(PSM, have[:a.max_runs] if a.max_runs else have)
    if a.max_runs:
        runs = runs[:a.max_runs]
    need = {src: set() for src in hfs}
    for run, psm, feat in runs:
        prov = np.load(os.path.join(MIX, run + ".prov.npy"))
        srun = np.load(os.path.join(MIX, run + ".srcrun.npy"))
        rows, spec = psm["pf_row"], psm["spec_idx"]
        dec = lab["is_decoy"][rows]
        tb = hash_tiebreak(psm["pf_id"], SEED)
        key = spec.astype(np.int64)
        scores = {"fixed": psm["hyperscore"].astype(float), "gbm": model.apply(psm, feat),
                  "gbm_noRT": model.apply(psm, neutral_rt(feat))}
        wins, thr[run] = {}, {}
        for p, sc in scores.items():
            w = psm_competition(spec, sc, tiebreak=tb)            # one winner row per spectrum with candidates
            wins[p] = w[np.argsort(spec[w], kind="stable")]
            r = double_competition(key, psm["pf_id"], lab["tdc_pair"][rows], sc, dec, seed=SEED, tiebreak=tb)
            thr[run][p] = {"1pct": GS.tail_threshold(r.score, r.is_decoy, r.q, 0.01),
                           "5pct": GS.tail_threshold(r.score, r.is_decoy, r.q, 0.05),
                           "n_accepted_1pct": int(((r.q <= 0.01) & ~r.is_decoy).sum())}
        sp_w = spec[wins["fixed"]]
        for p in PIPES:
            assert np.array_equal(spec[wins[p]], sp_w)
        pv = prov[sp_w]
        null = pv >= 1
        # label-symmetric: the winning row and its pair partner, for every pipeline
        pr = np.stack([rows[wins[p]] for p in PIPES] + [l1.partner[rows[wins[p]]] for p in PIPES], 1)[null]
        for src in hfs:
            m = pv[null] == src
            need[src].update(pr[m].ravel().tolist())
        rec["prov"].append(pv[null])
        rec["srcrun"].append(srun[sp_w][null])
        rec["run"].append(np.full(null.sum(), run))
        rec["pairrows"].append(pr)
        for p in PIPES:
            w = wins[p]
            s = scores[p][w]
            per[p]["win_t"].append(~dec[w][null])
            per[p]["score"].append(s[null])
            per[p]["tail1"].append((s >= thr[run][p]["1pct"])[null])
            per[p]["tail5"].append((s >= thr[run][p]["5pct"])[null])
            per[p]["logit"].append(clf.logit_rows(l1, rows[w][null], feat["rt_pred"][w][null]))
        print(run, "nulls", int(null.sum()), {p: thr[run][p]["n_accepted_1pct"] for p in PIPES}, round(time.time() - t0),
              flush=True)
    # pass 2: A1.3 homolog filter, one call per null source over every candidate pair seen
    near_rows = {}
    for src, hf in hfs.items():
        rr = np.array(sorted(x for x in need[src] if x >= 0), np.int64)
        near_rows[src] = rr[hf.near(seqs(rr))] if rr.size else rr
        print("homolog", PROV_NAME[src], rr.size, "near", near_rows[src].size, round(time.time() - t0), flush=True)
    pv_all = np.concatenate(rec["prov"])
    pr_all = np.concatenate(rec["pairrows"])
    near = np.zeros(pv_all.size, bool)
    for src, nr in near_rows.items():
        m = pv_all == src
        near[m] = np.isin(pr_all[m], nr).any(axis=1)
    rec["near"] = [near]
    rec["clean"] = [(pv_all >= 1) & ~near]
    del rec["pairrows"]
    cat = lambda x: np.concatenate(x)  # noqa: E731
    R = {k: cat(v) for k, v in rec.items()}
    P = {p: {k: cat(v) for k, v in d.items()} for p, d in per.items()}
    cl = R["clean"]
    out = {"n_mixture_runs": len(runs), "thresholds": thr,
           "null_counts": {"n_null_with_candidate": int(cl.size), "n_clean": int(cl.sum()),
                           "by_source": {PROV_NAME[int(u)]: {"n": int((R["prov"] == u).sum()),
                                                             "clean": int((cl & (R["prov"] == u)).sum())}
                                         for u in np.unique(R["prov"])},
                           "n_null_projects_species": int(len([u for u in np.unique(R["prov"]) if u >= 2]))},
           "learned": {}}
    F_ = P["fixed"]
    for p in ("gbm", "gbm_noRT"):
        L = P[p]
        res = {}
        res["G2"] = G2.g2_differential(F_["win_t"][cl], L["win_t"][cl], F_["score"][cl])
        res["G2_by_source"] = {PROV_NAME[int(u)]: G2.g2_differential(F_["win_t"][cl & (R["prov"] == u)],
                                                                    L["win_t"][cl & (R["prov"] == u)],
                                                                    F_["score"][cl & (R["prov"] == u)])["pass"]
                               for u in np.unique(R["prov"]) if (cl & (R["prov"] == u)).sum() > 100}
        res["frac_winner_changed"] = float(np.mean(L["win_t"][cl] != F_["win_t"][cl]))
        res["tail_G2_1pct"] = GS.tail_g2(L["tail1"][cl], L["win_t"][cl], F_["tail1"][cl], F_["win_t"][cl])
        res["tail_G2_5pct_diag"] = GS.tail_g2(L["tail5"][cl], L["win_t"][cl], F_["tail5"][cl], F_["win_t"][cl])
        res["tail_by_source_1pct"] = {PROV_NAME[int(u)]: {"n_tail": int((L["tail1"] & cl & (R["prov"] == u)).sum()),
                                                          "n_tail_target": int((L["tail1"] & L["win_t"] & cl & (R["prov"] == u)).sum()),
                                                          "n_tail_fixed": int((F_["tail1"] & cl & (R["prov"] == u)).sum()),
                                                          "n_tail_fixed_target": int((F_["tail1"] & F_["win_t"] & cl & (R["prov"] == u)).sum())}
                                      for u in np.unique(R["prov"])}
        res["tail_unfiltered_1pct_diag"] = GS.tail_g2(L["tail1"], L["win_t"], F_["tail1"], F_["win_t"])
        res["G5"] = GS.g5_gate(L["score"][cl], L["win_t"][cl])
        res["G5_fixed_reference"] = GS.g5_gate(F_["score"][cl], F_["win_t"][cl])
        res["G1"] = GS.g1_gate(L["score"][cl], L["logit"][cl])
        res["G1_fixed_reference"] = GS.g1_gate(F_["score"][cl], F_["logit"][cl])
        # POST HOC diagnostic (added after seeing the tail result; never changes a verdict): the same statistics
        # without the offset nulls and the S. aureus source, whose tail target winners are human skin/sweat proteins
        # (diag_source_tail.json), i.e. real human material rather than nulls.
        ex = cl & (R["prov"] != 1) & (R["prov"] != 5)
        res["POST_HOC_diag_excluding_offset_and_saureus"] = {
            "tail_G2_1pct": GS.tail_g2(L["tail1"][ex], L["win_t"][ex], F_["tail1"][ex], F_["win_t"][ex]),
            "G1": GS.g1_gate(L["score"][ex], L["logit"][ex]), "G1_fixed_reference": GS.g1_gate(F_["score"][ex], F_["logit"][ex]),
            "G5": GS.g5_gate(L["score"][ex], L["win_t"][ex])}
        res["G1_by_source_top_decile_rho"] = {PROV_NAME[int(u)]: GS.g1_gate(L["score"][cl & (R["prov"] == u)],
                                                                           L["logit"][cl & (R["prov"] == u)])["rho"]
                                              for u in np.unique(R["prov"]) if (cl & (R["prov"] == u)).sum() > 500}
        out["learned"][p] = res
        t = res["tail_G2_1pct"]
        print(p, "G2", res["G2"]["pass"], "tail", t["verdict"], t["n_tail"], round(t["frac_target"], 3),
              round(t["frac_target_ucb"], 3), round(t["diff"], 3), round(t["diff_ucb"], 3),
              "G5", round(res["G5"]["all"]["auc"], 4), "G1", round(res["G1"]["rho"], 4), res["G1"]["ci"], flush=True)
    out["seq_classifier"] = clf.info
    out["wall_s"] = time.time() - t0
    json.dump(out, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
