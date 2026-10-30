#!/usr/bin/env python3
"""DNIT steps 4-5 (NULL_RULES §9-§12): mixtures on the 29 D10cf hosts, IT-CID production search, L3 features (pinned
prediction store), final A1.3 (c), and the ONE gate run for the frozen itcid_linear_pred.

Statistics are A2's, imported unchanged: scripts/a2/gates_a2.evaluate / overall (tail-G2, G1, G2-original, G5, FS
handling), vulcan.gates.a2_stats (g3_gate_a2), scripts/a2/gates_a2.search_exhaustive_hyp, vulcan.gates.a2_fs
(certified_precursor, a2_hypotheses).  The pipelines are named as A2's evaluate() expects: "fixed" = hyperscore,
"gbm" = the scorer under test (here itcid_linear_pred), "gbm_noRT" = the same with RT features neutralized (reported).
IT-CID adaptations (declared in NULL_RULES before data): hosts = D10cf, search = 3c IT-CID parameters, no FS, scorer
features from vulcan.rescore.itcid + pinned store, G1 classifier rt_pred input = AlphaPeptDeep RT.

    python scripts/a2_itcid/gates_dnit.py build      # mixtures (no label read)
    python scripts/a2_itcid/gates_dnit.py search     # production IT-CID search of every mixture
    python scripts/a2_itcid/gates_dnit.py features   # L3 features, pinned store (fill: label-free Koina queries)
    python scripts/a2_itcid/gates_dnit.py finalize --record results/.../dnit_final_membership.json
    python scripts/a2_itcid/gates_dnit.py gates --record ... --out results/.../dnit_gates.json   # ONCE, after commit
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import pickle
import sys
import time

W = os.environ.setdefault("VULCAN_A2_WORK", "/tmp/vulcan_dnit")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
for p in (os.path.join(ROOT, "scripts", "m2"), os.path.join(ROOT, "scripts", "m3e"), os.path.join(ROOT, "scripts", "m3"),
          os.path.join(ROOT, "scripts", "a2"), os.path.join(ROOT, "scripts", "m3c"), HERE, ROOT):
    sys.path.insert(0, p)

import numpy as np  # noqa: E402
import yaml  # noqa: E402

import walk_dnit as WD  # noqa: E402

MIX = os.path.join(W, "mix")
PSM = os.path.join(W, "psms", "DNITmix")
P_HOST, P_OFF, P_A, P_B = 0, 1, 20, 21
PIPES = ("fixed", "gbm", "gbm_noRT")
POOL_CODE = {20: "A", 21: "B", 1: "OFF"}


def admitted():
    return json.load(open(os.path.join(ROOT, "data", "manifests", "dnit_select.log.json")))["admitted"]


def build(cfg):
    import common as C3
    from g2_nulls import concat_spectra, subset
    from vulcan.gates.a2_fs import certified_precursor
    from vulcan.gates.nulls import quantile_map
    from vulcan.gpu import engine as E
    seed = cfg["seeds"]["mixture_permutation"]
    per = cfg["mixtures"]["n_per_mixture"]
    hosts = [E.load_store(s) for s in C3.d10_stores()]
    H = len(hosts)
    parts, prov, srcrun, single, rt_map, info = [], [], [], [], [[] for _ in range(H)], {"sources": []}
    k = 0
    for acc in admitted():
        for f in sorted(glob.glob(os.path.join(W, "admission", acc, "*.npz"))):
            z = np.load(f)
            sp = E.load_store(str(z["store"]))
            codes = z["codes"]
            a_idx = z["A_spec"]
            b_idx = np.flatnonzero(codes != 1)
            for pv, part, idx, sgl in ((P_A, certified_precursor(sp, a_idx, z["A_mz_mono"], z["A_z"]), a_idx, True),
                                       (P_B, subset(sp, b_idx), b_idx, False)):
                parts.append(part)
                prov.append(np.full(len(part), pv, np.int8))
                srcrun.append(np.full(len(part), k, np.int16))
                single.append(np.full(len(part), sgl, bool))
                for h, host in enumerate(hosts):
                    rt_map[h].append(quantile_map(sp.rt, host.rt)[idx])
            info["sources"].append({"id": k, "accession": acc, "run": sp.run, "n_spectra": len(sp),
                                    "n_A": int(len(a_idx)), "n_B_candidates": int(len(b_idx))})
            k += 1
    pool = concat_spectra(parts, "pool")
    prov, srcrun, single = np.concatenate(prov), np.concatenate(srcrun), np.concatenate(single)
    rt_map = [np.concatenate(x) for x in rt_map]
    perm = np.random.default_rng(seed).permutation(len(pool))
    rng = np.random.default_rng(cfg["seeds"]["offsets_a1_recipe"])
    os.makedirs(MIX, exist_ok=True)
    import build_mixtures as BM                                 # A1 cut rule (scripts/a2/build_mixtures.cut)
    info["mixtures"] = {}
    for m, (a, b) in enumerate(BM.cut(len(pool), per)):
        h = m % H
        host = hosts[h]
        idx = np.sort(perm[a:b])
        ps, pv, sr, sg = [host], [np.zeros(len(host), np.int8)], [np.full(len(host), -1, np.int16)], \
            [np.zeros(len(host), bool)]
        if m < 2:
            known = np.flatnonzero(host.charge > 0)
            for _ in range(2):
                nn = rng.choice([-2, -1, 1, 2], size=len(known))
                delta = nn + 0.5 + rng.uniform(-0.1, 0.1, size=len(known))
                ps.append(subset(host, known, host.precursor_mz[known] + delta / host.charge[known].astype(float)))
                pv.append(np.full(len(known), P_OFF, np.int8))
                sr.append(np.full(len(known), -1, np.int16))
                sg.append(np.zeros(len(known), bool))
        nul = subset(pool, idx)
        nul.rt = rt_map[h][idx].astype(nul.rt.dtype)
        ps.append(nul)
        pv.append(prov[idx])
        sr.append(srcrun[idx])
        sg.append(single[idx])
        name = f"ditmix{m:03d}_h{h:02d}"
        mix = concat_spectra(ps, name)
        with open(os.path.join(MIX, name + ".spectra.pkl"), "wb") as fh:
            pickle.dump(mix, fh, protocol=5)
        np.save(os.path.join(MIX, name + ".prov.npy"), np.concatenate(pv))
        np.save(os.path.join(MIX, name + ".srcrun.npy"), np.concatenate(sr))
        np.save(os.path.join(MIX, name + ".single.npy"), np.concatenate(sg))
        info["mixtures"][name] = {"host": host.run, "n_spectra": len(mix)}
    info["n_pool"] = int(len(pool))
    info["n_by_code"] = {str(c): int((prov == c).sum()) for c in (P_A, P_B)}
    json.dump(info, open(os.path.join(MIX, "mix_info.json"), "w"), indent=1)
    print(json.dumps({k: info[k] for k in ("n_pool", "n_by_code")}), len(info["mixtures"]), "mixtures")


def _params():
    from vulcan.gpu import engine as E
    return E.SearchParams(**WD.SEARCH_IT)


def search(cfg):
    import common as C3
    import pyarrow.parquet as pq
    from vulcan.gates.a2_fs import a2_hypotheses
    from vulcan.gpu import engine as E
    from vulcan.gpu.host import pre_hypotheses
    from vulcan.index.build import L1Index
    os.makedirs(PSM, exist_ok=True)
    p = _params()
    l1 = L1Index(WD.GATE_INDEX)
    st, var = C3.mods("human")
    gi = E.GpuIndex(l1, p, st, var)
    for f in sorted(glob.glob(os.path.join(MIX, "*.spectra.pkl"))):
        name = os.path.basename(f)[:-len(".spectra.pkl")]
        out = os.path.join(PSM, name + ".psms.parquet")
        if os.path.exists(out):
            continue
        mix = pickle.load(open(f, "rb"))
        single = np.load(os.path.join(MIX, name + ".single.npy"))
        valid, _ = pre_hypotheses(mix, gi, p)
        hyp = a2_hypotheses(mix, valid, gi.mass32_np, p, single)
        t0 = time.perf_counter()
        tbl, stt, _ = E.search(mix, gi, l1, p, pre_hyp=(valid, hyp))
        if stt.counts.get("pre_hyp_used") != 1:
            raise RuntimeError(f"{name}: A2 hypotheses were not used")
        pq.write_table(tbl, out)
        print(name, len(mix), tbl.num_rows, round(time.perf_counter() - t0, 1), "s", flush=True)


def features(cfg):
    import pinned_cache as PC
    import pyarrow as pa
    import pyarrow.parquet as pq
    from vulcan.index.build import L1Index
    fb = PC.builder("human", "fill")
    l1 = L1Index(WD.GATE_INDEX)
    for f in sorted(glob.glob(os.path.join(PSM, "*.psms.parquet"))):
        name = os.path.basename(f)[:-len(".psms.parquet")]
        outp = os.path.join(PSM, name + ".feat.parquet")
        if os.path.exists(outp):
            continue
        mix = pickle.load(open(os.path.join(MIX, name + ".spectra.pkl"), "rb"))
        t = pq.read_table(f)
        psm = {c: t.column(c).to_numpy() for c in ("pf_row", "charge", "spec_idx", "rank", "hyperscore", "rt")}
        feats, info = fb.run(l1, psm, mix)
        pq.write_table(pa.table(feats), outp + ".partial")
        os.replace(outp + ".partial", outp)
        print(name, "features", round(info["timing_s"]["total"], 1), "s", flush=True)


# ------------------------------------------------------------------------------------------------ step 5 helpers
def _pipelines():
    import m3_data as MD
    from gates_a1 import _load, neutral_rt
    from vulcan.index.build import L1Index
    from vulcan.rescore.static import StaticModel
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import psm_competition
    it = yaml.safe_load(open(WD.CFG_IT))
    model = StaticModel.load(os.path.join(ROOT, it["scorer"]["file"]))
    l1 = L1Index(WD.GATE_INDEX)
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


def finalize(cfg):
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
    for run, psm, feat, prov, srun, scores, wins, l1, _ in _pipelines():
        rows = psm["pf_row"]
        sp_w = psm["spec_idx"][wins["fixed"]]
        pv = prov[sp_w]
        null = np.isin(pv, [P_A, P_B])
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
    from vulcan.gates import a2_config
    np.savez(os.path.join(W, "final_membership.npz"), **{k: v for k, v in zip(keys, keep_all)})
    rec = {"sha256": a2_config.sha256_file(os.path.join(W, "final_membership.npz")), "mixtures": keys,
           "n_dropped_final_c": int(sum((~k).sum() for k in keep_all)),
           "pinned_store_snapshot": __import__("pinned_cache").store("frozen").snapshot_id()}
    print(json.dumps({k: v for k, v in rec.items() if k != "mixtures"}))
    return rec


def apd_rt(l1, rows):
    """AlphaPeptDeep RT for the G1 classifier input (its training definition; label-blind sequence property)."""
    from vulcan.predict.apd import APD
    from vulcan.predict.cache import PredCache
    global _APD
    if "_APD" not in globals():
        _APD = PredCache(os.path.join(W, "apd_cache"), APD(), 29.0, "QE", {"C": 57.021464}, {"M": 15.9949})
    return _APD.lookup_rt(l1, rows)


def collect(cfg):
    import seqclf
    from vulcan.gates import stats as GS
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import double_competition
    clf = seqclf.load()
    fm = np.load(os.path.join(W, "final_membership.npz"))
    info = json.load(open(os.path.join(MIX, "mix_info.json")))
    src_proj = {s["id"]: s["accession"] for s in info["sources"]}
    R = {"pool": [], "project": []}
    P = {p: {k: [] for k in ("win_t", "score", "tail1", "logit")} for p in PIPES}
    lab = None
    for run, psm, feat, prov, srun, scores, wins, l1, MD in _pipelines():
        lab = lab or MD.labels(l1.path)
        rows = psm["pf_row"]
        dec = lab["is_decoy"][rows]
        tb = hash_tiebreak(psm["pf_id"], 11)
        sp_w = psm["spec_idx"][wins["fixed"]]
        keep = fm[run][sp_w]
        pv = prov[sp_w]
        m = keep & np.isin(pv, [P_OFF, P_A, P_B])
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
            rr = rows[w][m]
            P[p]["logit"].append(clf.logit_rows(l1, rr, apd_rt(l1, rr)))
        print(run, "collected", flush=True)
    return ({k: np.concatenate(v) for k, v in R.items()},
            {p: {k: np.concatenate(v) for k, v in d.items()} for p, d in P.items()})


def g3(cfg):
    import common as C3
    import m3_data as MD
    import torch
    from e2_correctness import _best_per_label, _rescored
    from gates_a2 import search_exhaustive_hyp
    from vulcan.gates import a2_stats as A2S
    from vulcan.gates.a2_fs import a2_hypotheses
    from vulcan.gpu import engine as E
    from vulcan.gpu.host import pre_hypotheses
    from vulcan.index.build import L1Index
    l1 = L1Index(WD.GATE_INDEX)
    isdec = MD.labels(WD.GATE_INDEX)["is_decoy"]
    p0 = _params()
    st, var = C3.mods("human")
    gi = E.GpuIndex(l1, p0, st, var)
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
        sel = np.flatnonzero(np.isin(prov, [P_A, P_B]) & fm[name])
        _, _, dx = search_exhaustive_hyp(sp, gi, l1, p0, hyp, sel)
        sx, hx, rx = _rescored(dx, l1, p0)
        btx, bdx = _best_per_label(sx, hx, isdec[rx], n)
        del dx
        _, _, dk = E.search(sp, gi, l1, p0, pre_hyp=(valid, hyp))
        sk, hk, rk = _rescored(dk, l1, p0)
        btk, bdk = _best_per_label(sk, hk, isdec[rk], n)
        del dk
        for k, v in (("bt_x", btx), ("bd_x", bdx), ("bt_k", btk), ("bd_k", bdk)):
            acc[k].append(v[sel])
        acc["pool"].append(np.where(prov[sel] == P_A, "A", "B"))
        torch.cuda.empty_cache()
        print(name, "g3", flush=True)
    A = {k: np.concatenate(v) for k, v in acc.items()}
    mA = A["pool"] == "A"
    gg = cfg["gates"]["g3"]
    return {"A": A2S.g3_gate_a2(*(np.where(mA, A[k], -np.inf) for k in ("bt_x", "bd_x", "bt_k", "bd_k")),
                                n_min=gg["n_min_per_label_per_bin"], tol=gg["tol"]),
            "AB": A2S.g3_gate_a2(A["bt_x"], A["bd_x"], A["bt_k"], A["bd_k"], n_min=gg["n_min_per_label_per_bin"],
                                 tol=gg["tol"])}


def main():
    import gates_a2 as GA
    from vulcan.gates import a2_config
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["build", "search", "features", "finalize", "gates"])
    ap.add_argument("--record")
    ap.add_argument("--out")
    a = ap.parse_args()
    base, it, cfg = WD.load_cfg()
    if a.what in ("build", "search", "features"):
        return globals()[a.what](cfg)
    a2_config.check_panel_files(base)
    a2_config.check_seqclf(base)
    if a.what == "finalize":
        json.dump(finalize(cfg), open(a.record, "w"), indent=1)
        return
    GA.require_final_membership(a.record)                     # committed + unmodified + same sha256 on disk
    if a2_config.sha256_file(os.path.join(ROOT, it["scorer"]["file"])) != it["scorer"]["sha256"]:
        raise SystemExit("frozen IT-CID scorer hash mismatch")
    t0 = time.time()
    R, P = collect(cfg)
    res = {"gates": GA.evaluate(R, P, cfg)}
    res["g3"] = g3(cfg)
    res["overall"] = GA.overall(res["gates"], res["g3"])
    res["naming"] = {"fixed": "hyperscore", "gbm": "itcid_linear_pred (frozen)", "gbm_noRT": "itcid_linear_pred, RT "
                     "features neutralized (reported only)"}
    res["wall_s"] = time.time() - t0
    json.dump(res, open(a.out, "w"), indent=1, default=float)
    print(json.dumps(res["overall"]))


if __name__ == "__main__":
    main()
