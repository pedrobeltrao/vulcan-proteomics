#!/usr/bin/env python3
"""3g step 2d (PLAN.md): causal ablation through the FROZEN static_gbm_pred (not retrained).  For feature group g,
every candidate's g-columns are set to that run's median over all its PSMs (M3e "without RT" convention, label-blind),
the frozen model is re-applied, then spectrum winners -> best evidence -> global paired competition -> the arm's own
1%/5% cuts -> E-vs-dec(E) coin.  Arm C (Arabidopsis) and arm P (shuffled control), E9 N = 1000 prefix.
    OMP_NUM_THREADS=1 conda_envs/pe/bin/python scripts/decoy_realism/ablation.py
"""
import json, os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
import numpy as np  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402
import common as C  # noqa: E402
from vulcan.rescore.static import StaticModel  # noqa: E402
from vulcan.stats import global_fdr as G  # noqa: E402
from vulcan.stats._util import hash_tiebreak  # noqa: E402

MODEL = os.path.join(C.ROOT, "results", "2026-10-27_M3_scoring", "models", "static_gbm_pred.pkl")
PSM_COLS = ["run", "spec_idx", "rank", "pf_row", "pf_id", "charge", "isotope_error", "hyperscore", "delta_next",
            "delta_best", "matched_b", "matched_y", "matched_peaks", "matched_intensity_pct", "average_ppm",
            "delta_mass_ppm", "longest_b", "longest_y", "longest_y_pct", "peptide_len", "missed_cleavages", "n_var",
            "scored_candidates", "poisson", "rt"]
PRED_MS2 = ["pred_sa", "pred_sa_sqrt", "pred_pearson", "pred_entropy", "pred_frac_obs", "pred_top6", "pred_obs_tic"]
RT = ["rt_abs_resid", "rt_resid_z"]
GROUPS = {
    "none": [],
    "pred_ms2": PRED_MS2, "rt": RT, "pred_all": PRED_MS2 + RT,
    **{f"only_{f}": [f] for f in PRED_MS2},
    "hs_family": ["hyperscore", "poisson"],
    "ladder": ["longest_b", "longest_y", "longest_y_pct", "matched_b", "matched_y", "matched_peaks"],
    "delta": ["delta_next", "delta_best"],
    "masscharge": ["average_ppm", "delta_mass_ppm", "isotope_error", "charge"],
    "other": ["matched_intensity_pct", "scored_candidates", "peptide_len", "missed_cleavages", "n_var"],
    "all_but_pred": ["hyperscore", "poisson", "longest_b", "longest_y", "longest_y_pct", "matched_b", "matched_y",
                     "matched_peaks", "delta_next", "delta_best", "matched_intensity_pct"],
}
GROUPS_P = ["none", "pred_ms2", "rt", "pred_all", "hs_family", "ladder", "all_but_pred"]
_M = None


def job(args):
    global _M
    if _M is None:
        _M = StaticModel.load(MODEL)
    arm, run, i, groups = args
    od = os.path.join("/tmp/vulcan_e9", "out", arm)
    f = os.path.join(od, run + ".psms.parquet")
    if not os.path.exists(f):
        return None
    t = pq.read_table(f, columns=PSM_COLS)
    psm = {c: t.column(c).to_numpy(zero_copy_only=False) for c in PSM_COLS}
    ft = pq.read_table(os.path.join(od, run + ".feat.parquet"))
    feat = {c: ft.column(c).to_numpy() for c in ft.column_names}
    stored = pq.read_table(os.path.join(od, run + ".scored.parquet"), columns=["score"]).column("score").to_numpy()
    tb = hash_tiebreak(psm["pf_id"], C.SEED)
    out = {}
    for g in groups:
        p2, f2 = dict(psm), dict(feat)
        for c in GROUPS[g]:
            src = p2 if c in p2 else f2
            v = src[c]
            src[c] = np.full_like(v, np.median(v))
        sc = _M.apply(p2, f2).astype(np.float32).astype(np.float64)
        if g == "none":
            out["_maxdiff_vs_stored"] = float(np.max(np.abs(sc - stored.astype(np.float64)))) if len(sc) else 0.0
        w = G.psm_winners(psm["spec_idx"], sc, tb)
        out[g] = (psm["pf_row"][w].astype(np.int64), psm["pf_id"][w], sc[w])
    out["_run"] = i
    return out


def main():
    man = C.manifest()
    runs = [r["run"] for r in man]
    res = {}
    for arm, groups in (("C", list(GROUPS)), ("P", GROUPS_P)):
        t0 = time.time()
        lab = C.labels(C.IDX[arm])
        fam = lab["family"]
        with ProcessPoolExecutor(int(os.environ.get("ABL_PROCS", "12"))) as ex:
            outs = [o for o in ex.map(job, [(arm, r, i, groups) for i, r in enumerate(runs)], chunksize=2)
                    if o is not None]
        res[f"{arm}_maxdiff_vs_stored"] = max(o["_maxdiff_vs_stored"] for o in outs)
        print(arm, "scored", time.time() - t0, res[f"{arm}_maxdiff_vs_stored"], flush=True)
        for g in groups:
            d = {"run": np.concatenate([np.full(len(o[g][0]), o["_run"], np.int32) for o in outs]),
                 "pf_row": np.concatenate([o[g][0] for o in outs]), "pf_id": np.concatenate([o[g][1] for o in outs]),
                 "score": np.concatenate([o[g][2] for o in outs])}
            gl, cut = C.compete(d, lab)
            f = fam[d["pf_row"][gl.row]]
            r = {"n_targets_1pct": int(np.sum((gl.score >= cut[0.01]) & (f == 0)))}
            for al, c in cut.items():
                m = (gl.score >= c) & np.isin(f, (1, 3))
                r[str(al)] = C.coin(f[m] == 1)
                r[str(al)]["win_E"] = None
                np.save(os.path.join(C.SCR, f"abl_{arm}_{g}_{al}.npy"), (f[m] == 1))
            res[f"{arm}_{g}"] = r
            print(arm, g, r["0.05"]["coin"], r["0.01"]["coin"], r["n_targets_1pct"], flush=True)
            json.dump(res, open(os.path.join(C.SCR, "ablation.json"), "w"), indent=1, default=float)
        del outs


if __name__ == "__main__":
    main()
