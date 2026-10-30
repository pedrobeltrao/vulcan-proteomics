#!/usr/bin/env python3
"""Amendment A2 step 5: gates on DN2 (DESIGN §6.0b; statistics in vulcan.gates.a2_stats).  Runs ONCE.

Phases (each refuses to run unless the previous phase's record is committed and unmodified):
  --finalize  (after step 4, before any statistic; label-free)  For every mixture and every compared pipeline
              (hyperscore, static_gbm_pred, static_gbm_pred without RT) the PSM-level winner of every null spectrum;
              final A1.3 (c): a null spectrum is dropped if, for ANY pipeline, the winner or its pair partner
              (l1.partner) is A1.3-near its project's reference (own proteome ∪ cRAP ∪ Frankenfield ∪ Panel H ∪
              Panel B).  Applied to A, B candidates and FS (stricter for FS); A1 offsets unfiltered (as A1).
              Writes <work>/final_membership.npz + final_membership.json (sha256, label-free counts only).
  --gates     (step 5)  Verifies config/a2.yaml, the frozen scorer and seq-classifier hashes, and that
              results/<date>_A2_admission/final_membership.json is committed with the same sha256; then computes:
                tail-G2 (A2.5) on A; reported: per project, leave-one-project-out, without RT, B, FS, A1 offsets
                G1 (A2.6) on A with the hyperscore control; FS if evaluable; reported: within-label rho, B, offsets
                G2-original and G5 on A ∪ B; G5 also on A and FS (FS: evaluable iff >= 200 target and >= 200 decoy winners; DESIGN §6.0b point 1)
  --g3        (step 5)  exhaustive vs production top-K over ALL DN2 mixtures with the A2 hypotheses (A2.8); pools A and
              A ∪ B.
The FDR-claim verdict = tail-G2(A) pass AND G1(A) pass AND G3(A) pass AND G3(A∪B) pass AND G2-original pass AND every
G5 pass AND no evaluable FS gate fails.  Not evaluable blocks (P3).
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import pickle
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3e"))

import numpy as np  # noqa: E402

from vulcan.gates import a2_config  # noqa: E402
from vulcan.gates import a2_stats as A2S  # noqa: E402
from vulcan.gates import stats as GS  # noqa: E402

W = os.environ.get("VULCAN_A2_WORK", "/tmp/vulcan_a2")
MIX = os.path.join(W, "mix")
PSM = os.path.join(W, "psms", "A2mix")
PIPES = ("fixed", "gbm", "gbm_noRT")
POOL_CODE = {20: "A", 21: "B", 22: "FS", 1: "OFF"}


# ------------------------------------------------------------------------------------------------ pure evaluation
def evaluate(R: dict, P: dict, cfg: dict) -> dict:
    """R: per null spectrum 'pool' (A/B/FS/OFF), 'project'; P[pipe]: 'win_t', 'score', 'tail1', 'logit' aligned with R.
    Returns every gate result and the overall verdict."""
    from vulcan.stats import g2 as G2
    g = cfg["gates"]
    pool, proj = np.asarray(R["pool"]), np.asarray(R["project"])
    F = P["fixed"]
    out = {"n_by_pool": {k: int((pool == k).sum()) for k in ("A", "B", "FS", "OFF")}, "learned": {}}

    def tail(m, L):
        return A2S.tail_g2_a2(proj[m], L["tail1"][m], L["win_t"][m], F["tail1"][m], F["win_t"][m],
                              B=g["tail_g2"]["bootstrap_B"], seed=g["tail_g2"]["bootstrap_seed"])

    def g1(m, L):
        return A2S.g1_gate_a2(L["score"][m], L["logit"][m], F["score"][m], F["logit"][m],
                              n_min=g["g1"]["n_min_top"])

    for p in ("gbm", "gbm_noRT"):
        L = P[p]
        mA, mB, mFS, mOFF = (pool == "A"), (pool == "B"), (pool == "FS"), (pool == "OFF")
        mAB = mA | mB
        res = {"tail_G2_A": tail(mA, L), "G1_A": g1(mA, L)}
        res["G2_original_AB"] = G2.g2_differential(F["win_t"][mAB], L["win_t"][mAB], F["score"][mAB])
        res["G5_AB"] = GS.g5_gate(L["score"][mAB], L["win_t"][mAB], alarm=g["g2_original_g5"]["g5_alarm"])
        res["G5_A"] = GS.g5_gate(L["score"][mA], L["win_t"][mA], alarm=g["g2_original_g5"]["g5_alarm"])
        fs_t, fs_g1 = tail(mFS, L), g1(mFS, L)
        nT, nD = int(L["win_t"][mFS].sum()), int((~L["win_t"][mFS]).sum())
        fs_g5 = GS.g5_gate(L["score"][mFS], L["win_t"][mFS], alarm=g["g2_original_g5"]["g5_alarm"]) \
            if (nT >= g["g2_original_g5"]["fs_g5_min_winners"] and nD >= g["g2_original_g5"]["fs_g5_min_winners"]) else {"pass": None, "note": "not evaluable", "nT": nT, "nD": nD}
        res["FS"] = {"tail_G2": fs_t, "G1": fs_g1, "G5": fs_g5}
        # reported only
        res["reported"] = {
            "tail_G2_B": tail(mB, L), "tail_G2_OFF": tail(mOFF, L), "G1_B": g1(mB, L), "G1_OFF": g1(mOFF, L),
            "G1_A_within_label": A2S.fisher_pooled_rho(
                [(L["score"][mA & L["win_t"]], L["logit"][mA & L["win_t"]]),
                 (L["score"][mA & ~L["win_t"]], L["logit"][mA & ~L["win_t"]])]),
            "tail_G2_A_per_project": {str(u): tail(mA & (proj == u), L)["a1_statistics"] for u in np.unique(proj[mA])},
            "tail_G2_A_leave_one_project_out": {str(u): tail(mA & (proj != u), L)["verdict"]
                                                for u in np.unique(proj[mA])}}
        fs_fail = any(x.get("verdict") == "fail" for x in (fs_t, fs_g1)) or fs_g5.get("pass") is False
        res["blocks"] = {"tail_G2_A": res["tail_G2_A"]["verdict"] != "pass", "G1_A": res["G1_A"]["verdict"] != "pass",
                         "G2_original_AB": not res["G2_original_AB"]["pass"], "G5_AB": not res["G5_AB"]["pass"],
                         "G5_A": not res["G5_A"]["pass"], "FS_evaluable_fail": bool(fs_fail)}
        out["learned"][p] = res
    return out


def overall(gates: dict, g3: dict | None) -> dict:
    """FDR-claim verdict for the frozen scorer (gbm; the noRT ablation is reported)."""
    b = dict(gates["learned"]["gbm"]["blocks"])
    if g3 is None:
        b["G3"] = True
    else:
        b["G3_A"] = g3["A"]["verdict"] != "pass"
        b["G3_AB"] = g3["AB"]["verdict"] != "pass"
    return {"blocked_by": sorted(k for k, v in b.items() if v), "fdr_claims_allowed": not any(b.values())}


# ------------------------------------------------------------------------------------------------ guards
def committed_unmodified(path: str) -> bool:
    rel = os.path.relpath(path, ROOT)
    ls = subprocess.run(["git", "-C", ROOT, "ls-files", "--error-unmatch", rel], capture_output=True)
    st = subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "--", rel], capture_output=True, text=True)
    return ls.returncode == 0 and not st.stdout.strip()


def require_final_membership(record: str) -> dict:
    if not committed_unmodified(record):
        raise SystemExit(f"gates_a2: {record} is not committed/unmodified: commit the final membership first (P2)")
    rec = json.load(open(record))
    got = a2_config.sha256_file(os.path.join(W, "final_membership.npz"))
    if got != rec["sha256"]:
        raise SystemExit("gates_a2: final membership on disk differs from the committed record")
    return rec


# ------------------------------------------------------------------------------------------------ data phases
def _pipelines(record_dir):  # pragma: no cover (GPU/data)
    """Per mixture: winners of every pipeline (label-free part)."""
    import m3_common as C
    import m3_data as MD
    from gates_a1 import _load, neutral_rt
    from vulcan.index.build import L1Index
    from vulcan.rescore.static import StaticModel
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import psm_competition
    model = StaticModel.load(os.path.join(C.RES, "models", "static_gbm_pred.pkl"))
    l1 = L1Index(C.index_dir("human_contam", "prod"))
    names = sorted(os.path.basename(f)[:-len(".feat.parquet")] for f in glob.glob(os.path.join(PSM, "*.feat.parquet")))
    for run, psm, feat in _load(PSM, names):
        prov = np.load(os.path.join(MIX, run + ".prov.npy"))
        srun = np.load(os.path.join(MIX, run + ".srcrun.npy"))
        spec = psm["spec_idx"]
        tb = hash_tiebreak(psm["pf_id"], 11)
        scores = {"fixed": psm["hyperscore"].astype(float), "gbm": model.apply(psm, feat),
                  "gbm_noRT": model.apply(psm, neutral_rt(feat))}
        wins = {}
        for p, sc in scores.items():
            w = psm_competition(spec, sc, tiebreak=tb)
            wins[p] = w[np.argsort(spec[w], kind="stable")]
        yield run, psm, feat, prov, srun, scores, wins, l1, MD


def finalize(cfg):  # pragma: no cover (GPU/data)
    import own_search as OS
    from vulcan.gates.homolog import HomologFilter
    info = json.load(open(os.path.join(MIX, "mix_info.json")))
    src_proj = {s["id"]: s["accession"] for s in info["sources"]}
    refs = {}
    for acc in sorted(set(src_proj.values())):
        upid = json.load(open(os.path.join(W, "admission", acc, "counts.json")))["upid"]
        _, src = OS.build_own_index(upid, cfg)
        refs[acc] = HomologFilter(src["own"] + src["contam"])
    keep_all, keys = [], []
    for run, psm, feat, prov, srun, scores, wins, l1, _ in _pipelines(None):
        rows = psm["pf_row"]
        sp_w = psm["spec_idx"][wins["fixed"]]
        pv = prov[sp_w]
        null = np.isin(pv, [20, 21, 22])
        pr = np.stack([rows[wins[p]] for p in PIPES] + [l1.partner[rows[wins[p]]] for p in PIPES], 1)
        near = np.zeros(sp_w.size, bool)
        for acc, hf in refs.items():
            m = null & np.isin(srun[sp_w], [k for k, v in src_proj.items() if v == acc])
            if m.any():
                rr = np.unique(pr[m])
                bad = rr[hf.near([l1.sequence(int(r)) for r in rr])]
                near[m] = np.isin(pr[m], bad).any(axis=1)
        keep = np.ones(prov.size, bool)
        keep[sp_w[null & near]] = False
        keep_all.append(keep)
        keys.append(run)
    np.savez(os.path.join(W, "final_membership.npz"), **{k: v for k, v in zip(keys, keep_all)})
    rec = {"sha256": a2_config.sha256_file(os.path.join(W, "final_membership.npz")), "mixtures": keys,
           "n_dropped_final_c": int(sum((~k).sum() for k in keep_all))}
    print(json.dumps(rec))
    return rec


def collect(cfg):  # pragma: no cover (GPU/data)
    import seqclf
    clf = seqclf.load()
    fm = np.load(os.path.join(W, "final_membership.npz"))
    info = json.load(open(os.path.join(MIX, "mix_info.json")))
    src_proj = {s["id"]: s["accession"] for s in info["sources"]}
    R = {"pool": [], "project": []}
    P = {p: {k: [] for k in ("win_t", "score", "tail1", "logit")} for p in PIPES}
    from vulcan.stats.tdc import double_competition
    from vulcan.stats._util import hash_tiebreak
    lab = None
    for run, psm, feat, prov, srun, scores, wins, l1, MD in _pipelines(None):
        lab = lab or MD.labels(l1.path)
        rows = psm["pf_row"]
        dec = lab["is_decoy"][rows]
        tb = hash_tiebreak(psm["pf_id"], 11)
        sp_w = psm["spec_idx"][wins["fixed"]]
        keep = fm[run][sp_w]
        pv = prov[sp_w]
        m = keep & np.isin(pv, [1, 20, 21, 22])
        R["pool"].append(np.array([POOL_CODE[int(x)] for x in pv[m]]))
        R["project"].append(np.array([src_proj.get(int(s), "host-offset") for s in srun[sp_w][m]]))
        for p in PIPES:
            w = wins[p]
            r = double_competition(psm["spec_idx"].astype(np.int64), psm["pf_id"], lab["tdc_pair"][rows], scores[p],
                                   dec, seed=11, tiebreak=tb)
            thr = GS.tail_threshold(r.score, r.is_decoy, r.q, 0.01)
            s = scores[p][w]
            P[p]["win_t"].append(~dec[w][m])
            P[p]["score"].append(s[m])
            P[p]["tail1"].append((s >= thr)[m])
            P[p]["logit"].append(clf.logit_rows(l1, rows[w][m], feat["rt_pred"][w][m]))
    return ({k: np.concatenate(v) for k, v in R.items()},
            {p: {k: np.concatenate(v) for k, v in d.items()} for p, d in P.items()})


def search_exhaustive_hyp(sp, gi, l1, p, hyp, spec_subset, batch: int = 20_000_000):  # pragma: no cover (GPU)
    """vulcan.gpu.engine.search_exhaustive with given (A2) hypotheses; returns the K2 debug dict like engine.search
    (consumed by scripts/m2/e2_correctness._rescored, as scripts/m3e/g3.py)."""
    import torch
    from vulcan.gpu import engine as E
    st = E.Stats()
    peaks = E.k0_preprocess(sp, p, gi.dev, st)
    o_h, o_r, o_c, scored, matched, d_spec, d_maxfz = E.run_k1(gi, hyp, peaks, p, st)
    h, r = E.window_candidates(hyp, gi.dev, spec_subset)
    spec = d_spec[h].long()
    parts, i, n = [], 0, len(h)
    while i < n:
        j = min(n, i + batch)
        if j < n:
            j = int(torch.searchsorted(spec, spec[j - 1], side="right"))
        parts.append(E.run_k2(gi, h[i:j], r[i:j], d_spec, d_maxfz, peaks, p))
        i = j
    k2 = {k: torch.cat([q[k] for q in parts]) for k in parts[0]}
    return None, st, {"peaks": peaks, "hyp": hyp, "h": h, "r": r, "k2": k2}


def g3(cfg):  # pragma: no cover (GPU/data)
    import m3_common as C
    import m3_data as MD
    import torch
    sys.path.insert(0, os.path.join(ROOT, "scripts", "m2"))
    from e2_correctness import _best_per_label, _rescored
    from vulcan.gates.a2_fs import a2_hypotheses
    from vulcan.gpu import engine as E
    from vulcan.gpu.host import pre_hypotheses
    from vulcan.index.build import L1Index
    idx = C.index_dir("human_contam", "prod")
    l1 = L1Index(idx)
    isdec = MD.labels(idx)["is_decoy"]
    p0 = E.SearchParams()
    gi = E.GpuIndex(l1, p0, C.STATIC, C.VAR)
    fm = np.load(os.path.join(W, "final_membership.npz"))
    acc = {k: [] for k in ("bt_x", "bd_x", "bt_k", "bd_k", "pool")}
    for f in sorted(glob.glob(os.path.join(MIX, "*.spectra.pkl"))):
        name = os.path.basename(f)[:-len(".spectra.pkl")]
        sp = pickle.load(open(f, "rb"))
        prov = np.load(os.path.join(MIX, name + ".prov.npy"))
        single = np.load(os.path.join(MIX, name + ".single.npy"))
        valid, _ = pre_hypotheses(sp, gi, p0)
        hyp = a2_hypotheses(sp, valid, gi.mass32_np, p0, single)
        n = len(sp)
        sel = np.flatnonzero(np.isin(prov, [20, 21]) & fm[name])
        _, _, dx = search_exhaustive_hyp(sp, gi, l1, p0, hyp, sel)
        sx, hx, rx = _rescored(dx, l1, p0)
        btx, bdx = _best_per_label(sx, hx, isdec[rx], n)
        _, _, dk = E.search(sp, gi, l1, p0, pre_hyp=(valid, hyp))
        sk, hk, rk = _rescored(dk, l1, p0)
        btk, bdk = _best_per_label(sk, hk, isdec[rk], n)
        for k, v in (("bt_x", btx), ("bd_x", bdx), ("bt_k", btk), ("bd_k", bdk)):
            acc[k].append(v[sel])
        acc["pool"].append(np.where(prov[sel] == 20, "A", "B"))
        torch.cuda.empty_cache()
    A = {k: np.concatenate(v) for k, v in acc.items()}
    mA = A["pool"] == "A"
    gg = cfg["gates"]["g3"]
    return {"A": A2S.g3_gate_a2(*(np.where(mA, A[k], -np.inf) for k in ("bt_x", "bd_x", "bt_k", "bd_k")),
                                n_min=gg["n_min_per_label_per_bin"], tol=gg["tol"]),
            "AB": A2S.g3_gate_a2(A["bt_x"], A["bd_x"], A["bt_k"], A["bd_k"], n_min=gg["n_min_per_label_per_bin"],
                                 tol=gg["tol"])}


def main():  # pragma: no cover
    ap = argparse.ArgumentParser()
    ap.add_argument("--finalize", action="store_true")
    ap.add_argument("--gates", action="store_true")
    ap.add_argument("--g3", action="store_true")
    ap.add_argument("--record", required=True, help="results/<date>_A2_admission/final_membership.json")
    ap.add_argument("--out")
    a = ap.parse_args()
    cfg = a2_config.load()
    a2_config.check_panel_files(cfg)
    a2_config.check_seqclf(cfg)
    if a.finalize:
        rec = finalize(cfg)
        json.dump(rec, open(a.record, "w"), indent=1)
        return
    require_final_membership(a.record)
    if a2_config.sha256_file(os.path.join(ROOT, cfg["gates"]["frozen_scorer"]["file"])) != \
            cfg["gates"]["frozen_scorer"]["sha256"]:
        raise SystemExit("gates_a2: frozen scorer hash mismatch")
    res = {}
    if a.gates:
        R, P = collect(cfg)
        res["gates"] = evaluate(R, P, cfg)
    if a.g3:
        res["g3"] = g3(cfg)
    if a.gates:
        res["overall"] = overall(res["gates"], res.get("g3"))
    json.dump(res, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
