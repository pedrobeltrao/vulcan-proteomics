#!/usr/bin/env python3
"""Gates for the frozen cross-encoder ce_s (DESIGN §6.4; dev only).

  --g0   label swap + DB permutation: D1 searched against pfu_contam_E2001 and the E2 swap index
         (/tmp/vulcan_e2/g0/index_swapped); CE scores must be bit-identical (plus PSM tables, as in M3).
  --g2   on the M3e G2 null mixtures (scripts/m3e/g2mix.py; production search, top-5) with the M3e A1.3 clean-null
         construction (homolog filter over the winners/pair partners of every compared pipeline):
           G2-original  vulcan.stats.g2.g2_differential(fixed=hyperscore, learned=CE)  [GATE]
           G5           leak AUC of the CE score over clean null winners <= 0.55          [GATE]
           tail-G2, G1  computed on these nulls as DIAGNOSTICS ONLY: these are the D12 nulls on which A1 gates are
                        recorded as failed; verdicts need A2's fresh nulls (DN2), so they are 'pending'.
         static_gbm_pred is recomputed on the same clean set as the reference learned pipeline.
    python scripts/m7/gates_ce.py --g0 --g2 --out results/2026-10-28_M7_encoders/gates_ce.json
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3e"))
import m3_common as C  # noqa: E402
import m3_data as MD  # noqa: E402

import numpy as np  # noqa: E402

import e6_eval as E6  # noqa: E402

RES = E6.RES
PIPES = ("fixed", "gbm", "ce")


def g0(model_path):
    import torch
    from vulcan.encoders.score import load_model, score_table
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    model = load_model(model_path)
    base = C.index_dir("pfu_contam", "E2001")
    alt = os.path.join(C.E2.WORK, "g0", "index_swapped")
    sp = E.load_store(C.E2.store_path(C.DEV["D1"]["raw"][0]))
    p = E.SearchParams(report_psms=C.REPORT)
    outs = []
    for idx in (base, alt):
        l1 = L1Index(idx)
        gi = E.GpuIndex(l1, p, C.STATIC, C.VAR)
        tbl, _, _ = E.search(sp, gi, l1, p)
        del gi
        torch.cuda.empty_cache()
        psm = {c: tbl.column(c).to_numpy(zero_copy_only=False) for c in ("pf_row", "charge", "spec_idx", "pf_id")}
        s, _ = score_table(model, l1, psm, sp)
        s2, _ = score_table(model, l1, psm, sp)          # repeat: run-to-run determinism on a shared GPU
        outs.append((tbl, s, s2))
    (t0, a0, a1), (t1, b0, b1) = outs
    res = {"psm_tables_bit_identical": bool(t0.equals(t1)), "ce_scores_bit_identical": bool(np.array_equal(a0, b0)),
           "ce_scores_repeat_identical": bool(np.array_equal(a0, a1) and np.array_equal(b0, b1)),
           "n_psms": int(t0.num_rows), "model_sha256": model.sha256}
    res["pass"] = res["psm_tables_bit_identical"] and res["ce_scores_bit_identical"]
    return res


def g2():
    import gates_a1 as GA
    from vulcan.gates import stats as GS
    from vulcan.index.build import L1Index
    from vulcan.stats import g2 as G2
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import double_competition, psm_competition
    import seqclf
    t0 = time.time()
    model = E6.models()["static_gbm_pred"]
    idx = C.index_dir("human_contam", "prod")
    lab = MD.labels(idx)
    l1 = L1Index(idx)
    clf = seqclf.load()
    hfs = GA.homolog_filters()
    seq_cache = {}

    def seqs(rows):
        out = []
        for r in rows:
            r = int(r)
            if r not in seq_cache:
                seq_cache[r] = l1.sequence(r)
            out.append(seq_cache[r])
        return out
    have = sorted(os.path.basename(f)[:-len(".feat.parquet")] for f in glob.glob(os.path.join(GA.PSM, "*.feat.parquet")))
    runs = GA._load(GA.PSM, have)
    rec = {k: [] for k in ("prov", "run", "pairrows")}
    per = {p: {k: [] for k in ("win_t", "score", "tail1", "tail5", "logit")} for p in PIPES}
    thr = {}
    need = {src: set() for src in hfs}
    for run, psm, feat in runs:
        prov = np.load(os.path.join(GA.MIX, run + ".prov.npy"))
        rows, spec = psm["pf_row"], psm["spec_idx"]
        dec = lab["is_decoy"][rows]
        tb = hash_tiebreak(psm["pf_id"], GA.SEED)
        key = spec.astype(np.int64)
        scores = {"fixed": psm["hyperscore"].astype(float), "gbm": model.apply(psm, feat),
                  "ce": np.load(E6.ce_path(GA.PSM, run)).astype(np.float64)}
        wins, thr[run] = {}, {}
        for p, sc in scores.items():
            w = psm_competition(spec, sc, tiebreak=tb)
            wins[p] = w[np.argsort(spec[w], kind="stable")]
            r = double_competition(key, psm["pf_id"], lab["tdc_pair"][rows], sc, dec, seed=GA.SEED, tiebreak=tb)
            thr[run][p] = {"1pct": GS.tail_threshold(r.score, r.is_decoy, r.q, 0.01),
                           "5pct": GS.tail_threshold(r.score, r.is_decoy, r.q, 0.05),
                           "n_accepted_1pct": int(((r.q <= 0.01) & ~r.is_decoy).sum())}
        sp_w = spec[wins["fixed"]]
        pv = prov[sp_w]
        null = pv >= 1
        pr = np.stack([rows[wins[p]] for p in PIPES] + [l1.partner[rows[wins[p]]] for p in PIPES], 1)[null]
        for src in hfs:
            need[src].update(pr[pv[null] == src].ravel().tolist())
        rec["prov"].append(pv[null]); rec["run"].append(np.full(null.sum(), run)); rec["pairrows"].append(pr)
        for p in PIPES:
            w = wins[p]
            s = scores[p][w]
            per[p]["win_t"].append(~dec[w][null]); per[p]["score"].append(s[null])
            per[p]["tail1"].append((s >= thr[run][p]["1pct"])[null]); per[p]["tail5"].append((s >= thr[run][p]["5pct"])[null])
            per[p]["logit"].append(clf.logit_rows(l1, rows[w][null], feat["rt_pred"][w][null]))
        print(run, int(null.sum()), {p: thr[run][p]["n_accepted_1pct"] for p in PIPES}, round(time.time() - t0), flush=True)
    near_rows = {}
    for src, hf in hfs.items():
        rr = np.array(sorted(x for x in need[src] if x >= 0), np.int64)
        near_rows[src] = rr[hf.near(seqs(rr))] if rr.size else rr
    pv_all = np.concatenate(rec["prov"])
    pr_all = np.concatenate(rec["pairrows"])
    near = np.zeros(pv_all.size, bool)
    for src, nr in near_rows.items():
        m = pv_all == src
        near[m] = np.isin(pr_all[m], nr).any(axis=1)
    cl = (pv_all >= 1) & ~near
    P = {p: {k: np.concatenate(v) for k, v in d.items()} for p, d in per.items()}
    F_ = P["fixed"]
    out = {"n_mixture_runs": len(runs), "n_null_with_candidate": int(cl.size), "n_clean": int(cl.sum()), "learned": {}}
    for p in ("gbm", "ce"):
        L = P[p]
        res = {"G2": G2.g2_differential(F_["win_t"][cl], L["win_t"][cl], F_["score"][cl]),
               "frac_winner_changed_vs_hyperscore": float(np.mean(L["win_t"][cl] != F_["win_t"][cl])),
               "G5": GS.g5_gate(L["score"][cl], L["win_t"][cl]),
               "DIAG_tail_G2_1pct_D12_nulls": GS.tail_g2(L["tail1"][cl], L["win_t"][cl], F_["tail1"][cl], F_["win_t"][cl]),
               "DIAG_G1_D12_nulls": GS.g1_gate(L["score"][cl], L["logit"][cl]),
               "null_target_win_rate": float(L["win_t"][cl].mean())}
        ex = cl & (pv_all != 1) & (pv_all != 5)
        res["DIAG_excl_offset_saureus"] = {"tail_G2_1pct": GS.tail_g2(L["tail1"][ex], L["win_t"][ex], F_["tail1"][ex], F_["win_t"][ex]),
                                           "G1": GS.g1_gate(L["score"][ex], L["logit"][ex])}
        out["learned"][p] = res
        t = res["DIAG_tail_G2_1pct_D12_nulls"]
        print(p, "G2", res["G2"]["pass"], "G5", round(res["G5"]["all"]["auc"], 4), "| diag tail", t["verdict"], t["n_tail"],
              round(t["frac_target"], 3), "G1", round(res["DIAG_G1_D12_nulls"]["rho"], 4), flush=True)
    out["fixed_null_target_win_rate"] = float(F_["win_t"][cl].mean())
    out["wall_s"] = time.time() - t0
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--g0", action="store_true")
    ap.add_argument("--g2", action="store_true")
    ap.add_argument("--model", default=os.path.join(RES, "models", "ce_s.pt"))
    ap.add_argument("--out", default=os.path.join(RES, "gates_ce.json"))
    a = ap.parse_args()
    out = json.load(open(a.out)) if os.path.exists(a.out) else {}
    if a.g0:
        out["G0"] = g0(a.model)
        print("G0", out["G0"], flush=True)
        json.dump(out, open(a.out, "w"), indent=1, default=float)
    if a.g2:
        out["G2_G5"] = g2()
        json.dump(out, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
