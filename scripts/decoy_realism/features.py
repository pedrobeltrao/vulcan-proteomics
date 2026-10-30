#!/usr/bin/env python3
"""3g step 2a-2c (PLAN.md): per-feature E-vs-dec(E) asymmetry on E9 null competitions.

2a  tail: features of the best-evidence PSM of every E / dec(E) global pair winner at/above the 5% (and 1%) cut
    (tables from e9_tail.py), AUC(E winners vs dec(E) winners) per feature, bootstrap CI; arms C and P; frozen score and
    hyperscore tails.
2b  bulk: spectrum-level winners (by the frozen score) that are E or dec(E), first --bulk-runs runs, AUC per feature
    overall and by decile of the winner score.
2c  within-spectrum paired: spectra whose top-5 holds both an E row and the dec(E) row of the same pair; sign rate
    of feature(E) - feature(dec E) per feature (ties excluded) and of the frozen score.
    conda_envs/pe/bin/python scripts/decoy_realism/features.py
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402
import common as C  # noqa: E402
from vulcan.rescore.static import BASE_FEATURES, PRED_FEATURES, feature_matrix  # noqa: E402
from vulcan.stats import global_fdr as G  # noqa: E402
from vulcan.stats._util import hash_tiebreak  # noqa: E402

PSM_COLS = ["run", "spec_idx", "rank", "pf_row", "pf_id", "charge", "isotope_error", "hyperscore", "delta_next",
            "delta_best", "matched_b", "matched_y", "matched_peaks", "matched_intensity_pct", "average_ppm",
            "delta_mass_ppm", "longest_b", "longest_y", "longest_y_pct", "peptide_len", "missed_cleavages", "n_var",
            "scored_candidates", "poisson", "rt", "summed_intensity"]
NAMES = BASE_FEATURES + PRED_FEATURES
EXTRA = ["rt_pred", "summed_intensity"]
ALLF = NAMES + EXTRA
WORK = "/tmp/vulcan_e9"
LAB = {}


def run_table(arm, run):
    od = os.path.join(WORK, "out", arm)
    t = pq.read_table(os.path.join(od, run + ".psms.parquet"), columns=PSM_COLS)
    psm = {c: t.column(c).to_numpy(zero_copy_only=False) for c in PSM_COLS}
    ft = pq.read_table(os.path.join(od, run + ".feat.parquet"))
    feat = {c: ft.column(c).to_numpy() for c in ft.column_names}
    X = feature_matrix(psm, feat, NAMES)
    X = np.column_stack([X, feat["rt_pred"].astype(float), psm["summed_intensity"].astype(float)])
    sc = pq.read_table(os.path.join(od, run + ".scored.parquet"), columns=["score"]).column("score").to_numpy()
    return psm, X, sc.astype(float)


def key(spec, row):
    return (np.asarray(spec, np.int64) << 26) | np.asarray(row, np.int64)


def job_tail(args):
    arm, run, specs, rows = args
    psm, X, _ = run_table(arm, run)
    k = key(psm["spec_idx"], psm["pf_row"])
    o = np.argsort(k)
    q = key(specs, rows)
    pos = o[np.searchsorted(k[o], q)]
    assert np.all(k[pos] == q)
    return X[pos]


def job_bulk(args):
    arm, run = args
    lab = LAB[arm]
    if not os.path.exists(os.path.join(WORK, "out", arm, run + ".psms.parquet")):
        return None
    psm, X, sc = run_table(arm, run)
    fam = lab["family"][psm["pf_row"]]
    w = G.psm_winners(psm["spec_idx"], sc, hash_tiebreak(psm["pf_id"], C.SEED))
    m = np.isin(fam[w], (1, 3))
    bulk = (X[w][m], sc[w][m], fam[w][m] == 1)
    # within-spectrum pairs (all candidates of the spectrum, top-5)
    pair = lab["tdc_pair"][psm["pf_row"]]
    ent = np.isin(fam, (1, 3))
    idx = np.flatnonzero(ent)
    kk = (psm["spec_idx"][idx].astype(np.int64) << 0)
    df = pd.DataFrame({"i": idx, "spec": kk, "pair": pair[idx], "isE": fam[idx] == 1})
    gpb = df.groupby(["spec", "pair"])["isE"].agg(["sum", "count"])
    both = gpb[(gpb["sum"] == 1) & (gpb["count"] == 2)].index
    if len(both):
        d2 = df.set_index(["spec", "pair"]).loc[both].reset_index()
        e = d2[d2.isE].sort_values(["spec", "pair"]).i.values
        dd = d2[~d2.isE].sort_values(["spec", "pair"]).i.values
        diff = X[e] - X[dd]
        sdiff = sc[e] - sc[dd]
        # is either member the spectrum winner?
        win = np.zeros(len(psm["spec_idx"]), bool)
        win[w] = True
        pw = win[e] | win[dd]
    else:
        diff, sdiff, pw = np.zeros((0, X.shape[1])), np.zeros(0), np.zeros(0, bool)
    return bulk, (diff, sdiff, pw)


def auc_ci(pos, neg, n_boot=300, seed=C.SEED):
    a = C.auc(pos, neg)
    rng = np.random.default_rng(seed)
    rs = [C.auc(pos[rng.integers(0, len(pos), len(pos))], neg[rng.integers(0, len(neg), len(neg))])
          for _ in range(n_boot)]
    return [a, float(np.quantile(rs, .025)), float(np.quantile(rs, .975))]


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--bulk-runs", type=int, default=200)
    ap.add_argument("--procs", type=int, default=16)
    a = ap.parse_args()
    man = C.manifest()
    runs = [r["run"] for r in man]
    res = {"features": ALLF}
    for arm in ("C", "P"):
        LAB[arm] = C.labels(C.IDX[arm])
    for arm in ("C", "P"):
        for col in ("score", "hyperscore"):
            tab = pd.read_parquet(os.path.join(C.SCR, f"e9tail_{arm}_{col}.parquet"))
            jobs = [(arm, runs[r], g.spec_idx.values, g.pf_row.values) for r, g in tab.groupby("run")]
            order = np.concatenate([g.index.values for _, g in tab.groupby("run")])
            with ProcessPoolExecutor(a.procs) as ex:
                Xs = list(ex.map(job_tail, jobs, chunksize=2))
            X = np.zeros((len(tab), len(ALLF)))
            X[order] = np.concatenate(Xs)
            np.save(os.path.join(C.SCR, f"e9tail_{arm}_{col}.X.npy"), X)
            out = {}
            for lvl, msk in (("0.05", np.ones(len(tab), bool)), ("0.01", tab.above1.values)):
                e = tab.E_won.values & msk
                d = ~tab.E_won.values & msk
                out[lvl] = {f: auc_ci(X[e, j], X[d, j]) for j, f in enumerate(ALLF)}
                out[lvl]["_score"] = auc_ci(tab.score.values[e], tab.score.values[d])
                out[lvl]["_n"] = [int(e.sum()), int(d.sum())]
            res[f"tail_{arm}_{col}"] = out
            print(arm, col, "tail done", flush=True)
        # bulk + within-spectrum pairs
        with ProcessPoolExecutor(a.procs) as ex:
            outs = list(ex.map(job_bulk, [(arm, r) for r in runs[:a.bulk_runs]], chunksize=2))
        outs = [o for o in outs if o is not None]
        Xb = np.concatenate([o[0][0] for o in outs])
        sb = np.concatenate([o[0][1] for o in outs])
        eb = np.concatenate([o[0][2] for o in outs])
        q = np.quantile(sb, np.linspace(0, 1, 11))
        dec = np.clip(np.searchsorted(q, sb, side="right") - 1, 0, 9)
        bulk = {"n": int(len(sb)), "coin_all": C.coin(eb), "by_decile": {}}
        bulk["auc_all"] = {f: C.auc(Xb[eb, j], Xb[~eb, j]) for j, f in enumerate(ALLF)}
        for k in range(10):
            m = dec == k
            bulk["by_decile"][k] = {"coin": C.coin(eb[m]), "score_lo": float(q[k]),
                                    "auc": {f: C.auc(Xb[m & eb, j], Xb[m & ~eb, j]) for j, f in enumerate(ALLF)}}
        res[f"bulk_{arm}"] = bulk
        D = np.concatenate([o[1][0] for o in outs])
        S = np.concatenate([o[1][1] for o in outs])
        PW = np.concatenate([o[1][2] for o in outs])
        pr = {"n_pairs": int(len(S)), "n_pairs_with_winner": int(PW.sum())}
        for nm, msk in (("all", np.ones(len(S), bool)), ("winner", PW)):
            dd = {}
            for j, f in enumerate(ALLF + ["_score"]):
                v = S[msk] if f == "_score" else D[msk, j]
                nz = v[v != 0]
                dd[f] = {"frac_E_greater": float((nz > 0).mean()) if len(nz) else None, "n_nontie": int(len(nz)),
                         "ci95": C.wilson(int((nz > 0).sum()), len(nz)), "median_diff": float(np.median(v)) if len(v) else None}
            pr[nm] = dd
        res[f"paired_{arm}"] = pr
        print(arm, "bulk/paired done", bulk["coin_all"], pr["n_pairs"], flush=True)
        json.dump(res, open(os.path.join(C.SCR, "features.json"), "w"), indent=1, default=float)
    json.dump(res, open(os.path.join(C.SCR, "features.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
