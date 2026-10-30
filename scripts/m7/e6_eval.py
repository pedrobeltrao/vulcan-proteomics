#!/usr/bin/env python3
"""E6: frozen cross-encoder (ce_s) vs frozen static_gbm_pred, RESET-LDA-pred and hyperscore on dev D0, D1, the D8
holdout (first 10 runs) and D7 arm M.  Production estimator, unchanged (scripts/m3/evaluate.py code path: top-5 PSM
table, paired double competition, peptidoform units, picked protein; RESET pooled peptide-level only).

  --score   (GPU) CE logits for every PSM table of the existing M3/M3e searches -> $W/ce/<dir>/<run>.ce.npy
  --eval    IDs at own q, §6.3 verdicts (P, C, M; CI convention a decides, b reported), matched-FDP IDs and ratios

The CE score enters through the feature dict (feat["ce_s"]) so every downstream function is the M3/M3e one.
    python scripts/m7/e6_eval.py --score --eval --model results/2026-10-28_M7_encoders/models/ce_s.pt
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
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3e"))
import m3_common as C  # noqa: E402
import m3_data as MD  # noqa: E402

import numpy as np  # noqa: E402

W = os.environ.get("VULCAN_M7_WORK", "/tmp/vulcan_m7")
M3E = "/tmp/vulcan_m3e"
RES = os.path.join(ROOT, "results", "2026-10-28_M7_encoders")
KINDS = ("prod", "E2001", "E2002", "E2003", "C")
ALPHAS = (0.001, 0.01, 0.05)
GRID = np.concatenate([np.linspace(0, 0.002, 41)[1:], np.linspace(0.002, 0.08, 161)[1:]])
SCORERS = ("hyperscore", "static_gbm_pred", "ce_s", "reset_lda_pred")
N_BOOT = 2000


def psm_dirs():
    d = {}
    for ds in ("D0", "D1"):
        for k in KINDS:
            d[f"{ds}_{k}"] = (os.path.join(C.WORK, "psms", f"{ds}_{k}"), C.index_dir("human_contam" if ds == "D0" else "pfu_contam", k))
    for k in KINDS:
        d[f"D8_{k}"] = (os.path.join(M3E, "psms", f"D8_{k}"), C.index_dir("human_contam", k))
    for s in (2001, 2002, 2003):
        d[f"D7_M{s}"] = (os.path.join(M3E, "psms", f"D7_M{s}"), os.path.join(M3E, "index", f"yeast_M{s}"))
    d["G2mix"] = (os.path.join(M3E, "psms", "G2mix"), C.index_dir("human_contam", "prod"))
    return d


def find_spectra(run):
    import pickle
    from vulcan.gpu.engine import load_store
    for p in (f"/tmp/vulcan_e2/store/{run}.vspec.arrow", f"/tmp/vulcan_m3/store_null/{run}.vspec.arrow",
              f"{M3E}/store_d8/{run}.vspec.arrow"):
        if os.path.exists(p):
            return load_store(p)
    p = f"{M3E}/null_mix/{run}.spectra.pkl"
    if os.path.exists(p):
        return pickle.load(open(p, "rb"))
    raise FileNotFoundError(run)


def ce_path(d, run):
    return os.path.join(W, "ce", os.path.basename(d.rstrip("/")), run + ".ce.npy")


def score(model_path, tags=None):
    import pyarrow.parquet as pq
    from vulcan.encoders.score import load_model, score_table
    from vulcan.index.build import L1Index
    model = load_model(model_path)
    timing = {}
    for tag, (d, idx) in psm_dirs().items():
        if tags and tag not in tags:
            continue
        files = sorted(glob.glob(os.path.join(d, "*.psms.parquet")))
        todo = [f for f in files if not os.path.exists(ce_path(d, os.path.basename(f)[:-13]))]
        if not todo:
            continue
        l1 = L1Index(idx)
        cache = {}
        for f in todo:
            run = os.path.basename(f)[:-len(".psms.parquet")]
            t = pq.read_table(f, columns=["pf_row", "charge", "spec_idx"])
            psm = {k: t.column(k).to_numpy() for k in ("pf_row", "charge", "spec_idx")}
            if run not in cache:
                cache = {run: find_spectra(run)}
            s, info = score_table(model, l1, psm, cache[run])
            os.makedirs(os.path.dirname(ce_path(d, run)), exist_ok=True)
            np.save(ce_path(d, run), s)
            timing[f"{tag}/{run}"] = info
            print(tag, run, {k: round(v, 3) if isinstance(v, float) else v for k, v in info.items()}, flush=True)
    p = os.path.join(W, "ce", "timing.json")
    old = json.load(open(p)) if os.path.exists(p) else {}
    old.update(timing)
    old["model_sha256"] = model.sha256
    json.dump(old, open(p, "w"), indent=1)


# ------------------------------------------------------------------------------------------------ evaluation glue
class CEModel:
    names = ["ce_s"]

    def apply(self, psm, feat):
        return feat["ce_s"].astype(np.float64)


def patch_loader():
    orig = MD.load_runs

    def load_runs(d, with_feat=True):
        out = orig(d, with_feat)
        for run, p, f in out:
            if f is not None:
                f["ce_s"] = np.load(ce_path(d, run))
                assert len(f["ce_s"]) == len(p["pf_row"])
        return out
    MD.load_runs = load_runs


def models():
    from vulcan.rescore.static import StaticModel
    return {"static_gbm_pred": StaticModel.load(os.path.join(C.RES, "models", "static_gbm_pred.pkl")), "ce_s": CEModel()}


def grid_counts(searches, ws, est, r=1.0):
    """Per grid q: prod IDs and per-replicate entrapment count tuples. Returns ids (G,), counts (J, G, 4)."""
    import evaluate as EV
    from vulcan.stats import entrapment as ENT
    w = ws["prod"]
    tq = w["q"][~w["is_decoy"]]
    ids = np.array([(tq <= g).sum() for g in GRID], float)
    keys = [k for k in ws if k.startswith("E")] if est == "paired" else ["C"]
    cnt = np.zeros((len(keys), len(GRID), 4))
    for j, k in enumerate(keys):
        for gi, g in enumerate(GRID):
            c = EV.peptide_entrapment(searches[k], ws[k], g, est)
            cnt[j, gi] = c.totals()
    return ids, cnt


def fdp_curve(cnt, est, r=1.0):
    from vulcan.stats import entrapment as ENT
    if est == "paired":
        return ENT.paired(cnt[..., 0], cnt[..., 1], cnt[..., 2], cnt[..., 3]).mean(0)
    return ENT.combined(cnt[..., 0], cnt[..., 1], r).mean(0)


def matched(fdp, ids, f):
    over = np.flatnonzero(fdp > f)
    if len(over) == 0:
        return float(ids[-1]), float(GRID[-1])
    k = int(over[0])
    if k == 0:
        return 0.0, 0.0
    w = (f - fdp[k - 1]) / (fdp[k] - fdp[k - 1])
    return float(ids[k - 1] + w * (ids[k] - ids[k - 1])), float(GRID[k - 1] + w * (GRID[k] - GRID[k - 1]))


def mclass(sc):
    return {"hyperscore": "fixed", "ce_s": "sequence_capacity"}.get(sc, "sequence_capacity")


def eval_experiment(searches, sc, mods, r_c, with_protein=True):
    """One experiment (D0: 2 runs pooled; D1; one D8 run). Returns record + grid data."""
    import evaluate as EV
    from vulcan.stats import entrapment as ENT
    ws = {k: EV.scored(s, sc, mods) for k, s in searches.items()}
    w = ws["prod"]
    rec = {"ids": {}, "ids_unique_sequences": {}, "proteins": {}, "entrapment": {}}
    tseq = searches["prod"].lab["target_seq"][searches["prod"].psm["pf_row"][w["row"]]]
    for a in ALPHAS:
        acc = (w["q"] <= a) & ~w["is_decoy"]
        rec["ids"][str(a)] = int(acc.sum())
        rec["ids_unique_sequences"][str(a)] = int(np.unique(tseq[acc]).size)
    if with_protein and sc != "reset_lda_pred":
        pl = EV.protein_level(searches["prod"], w, ALPHAS)
        rec["proteins"] = {str(a): v["n_groups"] for a, v in pl.items()}
    else:
        rec["proteins"] = "n/a (RESET: pooled peptide level only, §6.4)" if sc == "reset_lda_pred" else "skipped"
    reps = {}
    for a in ALPHAS:
        reps[a] = {"P": [ENT.Replicate(EV.peptide_entrapment(searches[k], ws[k], a, "paired"), int(k[1:]), 1001)
                         for k in ("E2001", "E2002", "E2003")],
                   "C": [ENT.Replicate(EV.peptide_entrapment(searches["C"], ws["C"], a, "combined"), -1, 1001)]}
    gP = grid_counts(searches, ws, "paired")
    gC = grid_counts(searches, ws, "combined")
    if sc == "reset_lda_pred":
        rec["reset"] = {k: v["extra"] for k, v in ws.items()}
    return rec, reps, gP, gC


def verdicts(reps, r_c):
    from vulcan.stats import entrapment as ENT
    out = {}
    for a, rr in reps.items():
        out[str(a)] = {cv: {"P": ENT.evaluate_arm(rr["P"], a, "paired", n_boot=N_BOOT, ci_convention=cv),
                            "C": ENT.evaluate_arm(rr["C"], a, "combined", r=r_c, n_boot=N_BOOT, ci_convention=cv)}
                       for cv in ("a", "b")}
    return out


def slim(v):
    keep = ("mean_fdp", "ci95", "mean_lb", "lb_ci95", "mean_n_t", "verdict", "fdp_per_replicate")
    return {k: v[k] for k in keep}


def eval_d0d1(mods, out):
    import evaluate as EV
    for ds in ("D0", "D1"):
        fam = EV.FAM[ds]
        t0 = time.time()
        searches = {k: EV.Search(ds, k) for k in KINDS}
        r_c = C.c_arm_r(fam)
        res = out.setdefault(ds, {"r_C": r_c, "scorers": {}})
        curves = {}
        for sc in SCORERS:
            rec, reps, gP, gC = eval_experiment(searches, sc, mods, r_c)
            v = verdicts(reps, r_c)
            rec["entrapment"] = {a: {cv: {arm: slim(x) for arm, x in d.items()} for cv, d in vv.items()} for a, vv in v.items()}
            curves[sc] = (gP, gC)
            fP = fdp_curve(gP[1], "paired")
            fC = fdp_curve(gC[1], "combined", r_c)
            rec["matched"] = {f"P_{f}": matched(fP, gP[0], f) for f in (0.01, 0.05)}
            rec["matched"].update({f"C_{f}": matched(fC, gC[0], f) for f in (0.01, 0.05)})
            res["scorers"][sc] = rec
            e = rec["entrapment"]["0.01"]["a"]
            print(ds, sc, rec["ids"], "P", round(e["P"]["mean_fdp"], 4), e["P"]["verdict"], "C", round(e["C"]["mean_fdp"], 4),
                  e["C"]["verdict"], "matchedP1%", round(rec["matched"]["P_0.01"][0]), round(time.time() - t0), flush=True)
        for comp in ("static_gbm_pred", "reset_lda_pred", "hyperscore"):
            res.setdefault("ratio_ce_vs", {})[comp] = {
                k: (res["scorers"]["ce_s"]["matched"][k][0] / res["scorers"][comp]["matched"][k][0]
                    if res["scorers"][comp]["matched"][k][0] else None) for k in ("P_0.01", "P_0.05", "C_0.01", "C_0.05")}
            res["ratio_ce_vs"][comp]["own_q_0.01"] = res["scorers"]["ce_s"]["ids"]["0.01"] / res["scorers"][comp]["ids"]["0.01"]
        print(ds, "ratios", json.dumps(res["ratio_ce_vs"]), flush=True)


def eval_d8(mods, out):
    import d8 as D8
    from vulcan.rescore.static import feature_matrix
    from vulcan.stats import entrapment as ENT
    if not hasattr(D8.RunSearch, "matrix"):     # RESET needs the feature matrix (as evaluate.Search.matrix)
        D8.RunSearch.matrix = lambda self, names: np.vstack([feature_matrix(p, f, names) for _, p, f in self.runs])
    labs = {k: MD.labels(C.index_dir("human_contam", k)) for k in KINDS}
    r_c = C.c_arm_r("human_contam")
    runs = D8.runs()
    res = out.setdefault("D8", {"runs": runs, "r_C": r_c, "per_run": {}, "pooled": {}})
    pooled_reps = {sc: {a: {"P": {s: [] for s in (2001, 2002, 2003)}, "C": []} for a in ALPHAS} for sc in SCORERS}
    G = {sc: {"ids": [], "P": [], "C": []} for sc in SCORERS}
    t0 = time.time()
    for run in runs:
        searches = {k: D8.RunSearch(run, k, labs[k]) for k in KINDS}
        for sc in SCORERS:
            rec, reps, gP, gC = eval_experiment(searches, sc, mods, r_c)
            res["per_run"].setdefault(run, {})[sc] = {"ids": rec["ids"], "proteins": rec["proteins"]}
            G[sc]["ids"].append(gP[0]); G[sc]["P"].append(gP[1]); G[sc]["C"].append(gC[1])
            for a in ALPHAS:
                for rp in reps[a]["P"]:
                    c = rp.counts
                    pooled_reps[sc][a]["P"][rp.entrap_seed].append(ENT.EntrapmentCounts(
                        key=np.array([f"{run}|{x}" for x in c.key]), t=c.t, e=c.e, e_only=c.e_only, e_gt_t=c.e_gt_t))
                c = reps[a]["C"][0].counts
                pooled_reps[sc][a]["C"].append(ENT.EntrapmentCounts(key=np.array([f"{run}|{x}" for x in c.key]), t=c.t,
                                                                    e=c.e, e_only=c.e_only, e_gt_t=c.e_gt_t))
        print("D8", run, {sc: res["per_run"][run][sc]["ids"]["0.01"] for sc in SCORERS}, round(time.time() - t0), flush=True)

    def cat(lst):
        return ENT.EntrapmentCounts(key=np.concatenate([c.key for c in lst]), t=np.concatenate([c.t for c in lst]),
                                    e=np.concatenate([c.e for c in lst]), e_only=np.concatenate([c.e_only for c in lst]),
                                    e_gt_t=np.concatenate([c.e_gt_t for c in lst]))
    for sc in SCORERS:
        rec = {"ids": {str(a): int(sum(res["per_run"][r][sc]["ids"][str(a)] for r in runs)) for a in ALPHAS}, "entrapment": {}}
        for a in ALPHAS:
            rp = [ENT.Replicate(cat(pooled_reps[sc][a]["P"][s]), s, 1001) for s in (2001, 2002, 2003)]
            rc = [ENT.Replicate(cat(pooled_reps[sc][a]["C"]), -1, 1001)]
            rec["entrapment"][str(a)] = {cv: {"P": slim(ENT.evaluate_arm(rp, a, "paired", n_boot=N_BOOT, ci_convention=cv)),
                                              "C": slim(ENT.evaluate_arm(rc, a, "combined", r=r_c, n_boot=N_BOOT, ci_convention=cv))}
                                         for cv in ("a", "b")}
        ids = np.stack(G[sc]["ids"])            # (R, G)
        P = np.stack(G[sc]["P"])                # (R, J, G, 4)
        Cc = np.stack(G[sc]["C"])
        G[sc] = (ids, P, Cc)
        rec["matched"] = {f"P_{f}": matched(fdp_curve(P.sum(0), "paired"), ids.sum(0), f) for f in (0.01, 0.05)}
        rec["matched"].update({f"C_{f}": matched(fdp_curve(Cc.sum(0), "combined", r_c), ids.sum(0), f) for f in (0.01, 0.05)})
        res["pooled"][sc] = rec
        e = rec["entrapment"]["0.01"]["a"]
        print("D8 pooled", sc, rec["ids"], "P", round(e["P"]["mean_fdp"], 4), e["P"]["verdict"], "C", round(e["C"]["mean_fdp"], 4),
              e["C"]["verdict"], "matchedP1%", round(rec["matched"]["P_0.01"][0]), flush=True)
    # per-run bootstrap of matched-FDP ratios (runs resampled; pooled FDP curve recomputed in every draw)
    rng = np.random.default_rng(11)
    R = len(runs)
    draws = rng.multinomial(R, np.full(R, 1 / R), size=N_BOOT).astype(float)
    res["ratio_ce_vs"] = {}
    for comp in ("static_gbm_pred", "reset_lda_pred", "hyperscore"):
        rr = {}
        for arm, f in (("P", 0.01), ("P", 0.05), ("C", 0.01), ("C", 0.05)):
            def m_ids(sc, wv):
                ids, P, Cc = G[sc]
                cnt = np.tensordot(wv, P if arm == "P" else Cc, 1)
                return matched(fdp_curve(cnt, "paired" if arm == "P" else "combined", r_c), wv @ ids, f)[0]
            point = m_ids("ce_s", np.ones(R)) / m_ids(comp, np.ones(R))
            bs = np.array([m_ids("ce_s", wv) / max(m_ids(comp, wv), 1e-9) for wv in draws])
            rr[f"{arm}_{f}"] = {"ratio": point, "ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]}
        rr["own_q_0.01"] = res["pooled"]["ce_s"]["ids"]["0.01"] / res["pooled"][comp]["ids"]["0.01"]
        res["ratio_ce_vs"][comp] = rr
    print("D8 ratios", json.dumps(res["ratio_ce_vs"]), flush=True)


def eval_armM(mods, out):
    import armM as AM
    from vulcan.stats import entrapment as ENT
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import double_competition
    res = out.setdefault("D7_armM", {"runs": AM.RUNS, "r": 1.0, "scorers": {}})
    per = {}
    psm_seq = {}
    for s in AM.SEEDS:
        runs = MD.load_runs(AM.pdir(s))
        psm = MD.concat([p for _, p, _ in runs])
        lab = MD.labels(AM.idx(s))
        rows = psm["pf_row"]
        spec = np.array([f"{a}:{b}" for a, b in zip(psm["run"], psm["spec_idx"])])
        unit = psm["pf_id"]
        tb = hash_tiebreak(unit, 11)
        sc = {"hyperscore": psm["hyperscore"].astype(float),
              "static_gbm_pred": np.concatenate([mods["static_gbm_pred"].apply(p, f) for _, p, f in runs]),
              "ce_s": np.concatenate([mods["ce_s"].apply(p, f) for _, p, f in runs])}
        for name, score in sc.items():
            r = double_competition(spec, unit, lab["tdc_pair"][rows], score, lab["is_decoy"][rows], seed=11, tiebreak=tb)
            per.setdefault(name, {})[s] = (r, lab["is_entrapment"][rows][r.row], unit[r.row], lab["target_seq"][rows][r.row])
    train_pos = training_positive_seqs()
    for name, d in per.items():
        rec = {"ids_target_side": {}, "entrapment": {}, "memorization": {}}
        for a in ALPHAS:
            reps = []
            for s in AM.SEEDS:
                r, ent, u, tseq = d[s]
                ts = ~r.is_decoy
                acc = (r.q <= a) & ts
                reps.append(ENT.Replicate(ENT.combined_counts(u[ts], ent[ts], acc[ts]), s, 1001))
                rec["ids_target_side"].setdefault(str(a), {})[s] = {"N_T": int((acc & ~ent).sum()), "N_E": int((acc & ent).sum())}
                if a == 0.01:
                    es = np.array([x.replace("I", "L") for x in tseq[acc & ent]])
                    allent = np.array([x.replace("I", "L") for x in tseq[ts & ent]])
                    rec["memorization"][s] = {"n_accepted_entrapment": int(len(es)),
                                              "frac_accepted_entrapment_in_training_positives": float(np.mean([x in train_pos for x in es])) if len(es) else None,
                                              "frac_competing_entrapment_in_training_positives": float(np.mean([x in train_pos for x in allent])) if len(allent) else None}
            rec["entrapment"][str(a)] = {cv: slim(ENT.evaluate_arm(reps, a, "combined", r=1.0, n_boot=N_BOOT, ci_convention=cv))
                                         for cv in ("a", "b")}
        res["scorers"][name] = rec
        e = rec["entrapment"]["0.01"]["a"]
        print("armM", name, round(e["mean_fdp"], 5), [round(x, 5) for x in e["ci95"]], e["verdict"], "LB", round(e["mean_lb"], 5),
              rec["memorization"], flush=True)


_TP = None


def training_positive_seqs():
    """I/L-collapsed stripped sequences of all training positives (memorization diagnostic only)."""
    global _TP
    if _TP is None:
        z = np.load(os.path.join(W, "dataset", "train.npz"))
        tok = z["cand_tok"][z["cand_label"] == 1]
        from vulcan.encoders.model import OX_TOKEN
        tok = np.where(tok == OX_TOKEN, ord("M") - 64, tok)
        _TP = {bytes((t[t > 0] + 64).astype(np.uint8)).decode().replace("I", "L") for t in tok}
    return _TP


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--score", action="store_true")
    ap.add_argument("--eval", action="store_true")
    ap.add_argument("--model", default=os.path.join(RES, "models", "ce_s.pt"))
    ap.add_argument("--tags", nargs="*")
    ap.add_argument("--sets", nargs="+", default=["D0D1", "D8", "armM"])
    ap.add_argument("--out", default=os.path.join(RES, "e6.json"))
    a = ap.parse_args()
    if a.score:
        score(a.model, a.tags)
    if a.eval:
        import evaluate as EV
        patch_loader()
        EV.STATIC = list(EV.STATIC) + ["ce_s"]
        mods = models()
        out = json.load(open(a.out)) if os.path.exists(a.out) else {}
        out["model"] = a.model
        out["grid"] = GRID.tolist()
        if "D0D1" in a.sets:
            eval_d0d1(mods, out)
            json.dump(out, open(a.out, "w"), indent=1, default=float)
        if "armM" in a.sets:
            eval_armM(mods, out)
            json.dump(out, open(a.out, "w"), indent=1, default=float)
        if "D8" in a.sets:
            eval_d8(mods, out)
            json.dump(out, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
