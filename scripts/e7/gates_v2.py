#!/usr/bin/env python3
"""Gates G0, G2 (original, M3 design) and G5 for scorer v2 + the chimeric search (DESIGN §6.4; config.yaml gates_v2).
Dev only; A2 DN2 data are not read.  Tail-G2 / G1 / G3 are pending the A2 DN2 nulls.

G0  D1 (+ E2001) searched in chimeric mode against the E2 index and against the E2 role-swapped, row-shuffled index
    (/tmp/vulcan_e2/g0/index_swapped); chimeric PSM tables, all v2 features (fresh prediction cache per side), the v1
    and v2 scores must be bit-identical.
G2  The M3 G2 mixtures (scripts/m3/g2_nulls.py: each D0 run + 2 offset copies of its charge-known spectra + half of
    the P. furiosus (D1) and D7 yeast spectra; the mixing is reproduced from its seed and checked against the saved
    mixtures), searched in CHIMERIC mode against the human prod index.  Precursor groups come from each source run's
    own MS1 hypothesis table; offset copies shift EVERY group's precursor mass by the copy's offset (so secondary
    hypotheses are non-physical too).  Null query unit = (null spectrum, precursor group); fixed winner = hyperscore,
    learned = v2 (same candidates).  M3 homolog filter on species nulls.  vulcan.stats.g2.g2_differential (unchanged).
G5  leak_auc(v2 score of the null winners, target win) <= 0.55.

    /tmp/vulcan_env_<uid>/pred/bin/python scripts/e7/gates_v2.py --g0 --g2
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))

import numpy as np  # noqa: E402

W = os.environ.get("VULCAN_E7_WORK", "/tmp/vulcan_e7")
RES = os.path.join(ROOT, "results", "2026-10-28_E7_chimeric_v2")
T1W = f"/tmp/vulcan_t1_{os.getuid()}"
STATIC, VAR = {"C": 57.021464}, {"M": 15.9949}
SEED = 11
D0_RUNS = ["b1906_293T_proteinID_01A_QE3_122212", "b1937_293T_proteinID_01B_QE3_122212"]
YEAST = ["LM2480_DDA", "LM2481_DDA", "LM2486_DDA"]
NULL_DIR = "/tmp/vulcan_m3/store_null"
V1 = os.path.join(ROOT, "results", "2026-10-27_M3_scoring", "models", "static_gbm_pred.pkl")
V2 = os.path.join(RES, "models", "static_v2.pkl")


def ms1(store, out):
    run = os.path.basename(store)[:-len(".vspec.arrow")]
    if not os.path.exists(os.path.join(out, run + ".hypotheses.arrow")):
        os.makedirs(out, exist_ok=True)
        subprocess.run([sys.executable, "-m", "vulcan.ms1", store, out, "--threads", "8"], check=True, cwd=ROOT,
                       env={**os.environ, "PYTHONPATH": ROOT}, stdout=subprocess.DEVNULL)
    return os.path.join(out, run + ".hypotheses.arrow"), os.path.join(out, run + ".features.arrow")


def scorer(cache_root):
    from vulcan.predict.apd import APD
    from vulcan.predict.cache import PredCache
    from vulcan.rescore.features_v2 import FeatureBuilderV2
    pred = APD()
    caches = {}

    def fac(nce):
        if nce not in caches:
            caches[nce] = PredCache(cache_root, pred, nce, "QE", STATIC, VAR)
        return caches[nce]
    return FeatureBuilderV2(fac, STATIC, VAR)


def models():
    from vulcan.rescore.static import StaticModel
    return {"v1": StaticModel.load(V1), "v2": StaticModel.load(V2)}


# ------------------------------------------------------------------------------------------------ G0
def g0():
    import shutil

    import torch
    from vulcan.chimeric import search as CS
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    from vulcan.rescore.features_v2 import model_psm, read_ms1_features
    store = os.path.join(T1W, "D1", "stage", "Velos005137.vspec.arrow")
    hyp, feat = ms1(store, os.path.join(W, "ms1", "D1"))
    sp = E.load_store(store)
    groups = CS.precursor_groups(CS.read_hypotheses(hyp), sp)
    mf = read_ms1_features(feat)
    ms = models()
    p = E.SearchParams(report_psms=5)
    outs = []
    for side, idx in (("base", "/tmp/vulcan_e2/index/pfu_contam_E2001"), ("swapped", "/tmp/vulcan_e2/g0/index_swapped")):
        l1 = L1Index(idx)
        gi = E.GpuIndex(l1, p, STATIC, VAR)
        tbl, _ = CS.search_chimeric(sp, gi, l1, p, groups)
        del gi
        torch.cuda.empty_cache()
        psm = {c: tbl.column(c).to_numpy(zero_copy_only=False) for c in tbl.column_names}
        root = os.path.join(W, "g0_cache", side)
        shutil.rmtree(root, ignore_errors=True)
        fb = scorer(root)
        f, info = fb.run_v2(l1, psm, sp, mf)
        mp = model_psm(psm)
        outs.append((tbl, f, {k: m.apply(mp, f) for k, m in ms.items()}, info["nce"]))
    (t0, f0, s0, n0), (t1, f1, s1, n1) = outs
    def col_eq(c):
        a, b = t0.column(c).to_numpy(zero_copy_only=False), t1.column(c).to_numpy(zero_copy_only=False)
        if a.dtype.kind == "f":      # NaN == NaN (ms1_feature_rt_apex is NaN when no MS1 feature): bitwise compare
            return a.shape == b.shape and bool(np.array_equal(a.view(np.uint8), b.view(np.uint8)))
        return bool(np.array_equal(a, b))
    cols = {c: col_eq(c) for c in t0.column_names} if t0.column_names == t1.column_names else {}
    res = {"chimeric_psm_tables_bit_identical": bool(cols) and all(cols.values()),
           "columns_not_identical": [c for c, v in cols.items() if not v], "arrow_equals_raw": bool(t0.equals(t1)),
           "nce_same": n0 == n1,
           "features_bit_identical": {k: bool(np.array_equal(f0[k], f1[k])) for k in f0},
           "scores_bit_identical": {k: bool(np.array_equal(s0[k], s1[k])) for k in s0}, "n_psms": t0.num_rows,
           "n_secondary_psms": int((t0.column("group").to_numpy() > 0).sum())}
    res["pass"] = res["chimeric_psm_tables_bit_identical"] and res["nce_same"] and \
        all(res["features_bit_identical"].values()) and all(res["scores_bit_identical"].values())
    return res


# ------------------------------------------------------------------------------------------------ G2 / G5
def concat_spectra(parts, run):
    from vulcan.gpu.engine import Spectra
    offs, base = [], 0
    for p in parts:
        offs.append(p.offsets[:-1] + base)
        base += p.offsets[-1]
    return Spectra(run, np.concatenate([p.scan for p in parts]), np.concatenate([p.index for p in parts]),
                   np.concatenate([p.rt for p in parts]), np.concatenate([p.precursor_mz for p in parts]),
                   np.concatenate([p.charge for p in parts]), np.concatenate([p.mz for p in parts]),
                   np.concatenate([p.intensity for p in parts]), np.r_[np.concatenate(offs), base].astype(np.int64))


def mixtures():
    """Reproduce scripts/m3/g2_nulls.py's mixing (seed 11) -> per mixture: (mix Spectra, provenance, source run id,
    source spectrum index, mass offset Da) and the source runs."""
    from vulcan.gpu import engine as E
    src_names = D0_RUNS + ["Velos005137"] + YEAST
    stores = [os.path.join(T1W, "D0", "stage", r + ".vspec.arrow") for r in D0_RUNS] + \
        [os.path.join(T1W, "D1", "stage", "Velos005137.vspec.arrow")] + [os.path.join(NULL_DIR, y + ".vspec.arrow") for y in YEAST]
    sps = [E.load_store(s) for s in stores]
    rng = np.random.default_rng(SEED)
    d0, foreign_parts = sps[:2], sps[2:]
    f_src = np.concatenate([np.full(len(p), 2 + i, np.int64) for i, p in enumerate(foreign_parts)])
    f_idx = np.concatenate([np.arange(len(p)) for p in foreign_parts])
    nf = len(f_src)
    halves = np.array_split(rng.permutation(nf), len(d0))
    out = []
    for i, (sp, half) in enumerate(zip(d0, halves)):
        src = [np.full(len(sp), i)]
        idx = [np.arange(len(sp))]
        dl = [np.zeros(len(sp))]
        known = np.flatnonzero(sp.charge > 0)
        for _ in range(2):
            n = rng.choice([-2, -1, 1, 2], size=len(known))
            delta = n + 0.5 + rng.uniform(-0.1, 0.1, size=len(known))
            src.append(np.full(len(known), i))
            idx.append(known)
            dl.append(delta)
        half = np.sort(half)
        src.append(f_src[half])
        idx.append(f_idx[half])
        dl.append(np.zeros(len(half)))
        mix = pickle.load(open(os.path.join(NULL_DIR, sp.run + "_g2mix.spectra.pkl"), "rb"))
        prov = np.load(os.path.join(NULL_DIR, sp.run + "_g2mix.provenance.npy"))
        src, idx, dl = np.concatenate(src), np.concatenate(idx), np.concatenate(dl)
        exp_pmz = np.array([sps[s].precursor_mz[k] for s, k in zip(src, idx)]) + dl / np.maximum(
            np.array([sps[s].charge[k] for s, k in zip(src, idx)], np.float64), 1)
        assert len(mix) == len(src) and np.allclose(mix.precursor_mz, exp_pmz, rtol=0, atol=1e-9), "mixing not reproduced"
        assert np.array_equal(np.diff(mix.offsets), np.array([np.diff(sps[s].offsets)[k] for s, k in zip(src, idx)]))
        out.append({"mix": mix, "prov": prov, "src": src, "idx": idx, "delta": dl})
    return out, sps, stores, src_names


def mixture_groups(m, sps, stores):
    """Precursor groups of a mixture from each source run's own hypothesis table; offset copies shift every group."""
    from vulcan.chimeric import search as CS
    from vulcan.rescore.features_v2 import read_ms1_features
    G, F, foff = [], [], [0]
    for s, st in zip(sps, stores):
        hyp, feat = ms1(st, os.path.join(W, "ms1", "g2src"))
        G.append(CS.precursor_groups(CS.read_hypotheses(hyp), s))
        F.append(read_ms1_features(feat))
        foff.append(foff[-1] + len(F[-1]["trace_cos"]))
    src, idx, dl = m["src"], m["idx"], m["delta"]
    n = len(src)
    prim = {k: np.zeros(n, G[0]["prim"][k].dtype) for k in G[0]["prim"]}
    for k in prim:
        prim[k] = np.concatenate([G[s]["prim"][k][[i]] for s, i in zip(src, idx)]) if n else prim[k]
    prim["feature_id"] = np.where(prim["feature_id"] >= 0, prim["feature_id"] + np.array([foff[s] for s in src]), -1)
    secs = []
    for j, (s, i, d) in enumerate(zip(src, idx, dl)):
        sec = G[s]["sec"]
        a, b = np.searchsorted(sec["spec"], i, "left"), np.searchsorted(sec["spec"], i, "right")
        if b > a:
            secs.append((j, s, d, {k: v[a:b] for k, v in sec.items()}))
    sec = {k: np.concatenate([x[3][k] for x in secs]) for k in G[0]["sec"]}
    sec["spec"] = np.concatenate([np.full(len(x[3]["spec"]), x[0], np.int64) for x in secs])
    dd = np.concatenate([np.full(len(x[3]["spec"]), x[2]) for x in secs])
    sec["mz"] = sec["mz"] + dd / sec["charge"].astype(np.float64)
    sec["feature_id"] = sec["feature_id"] + np.concatenate([np.full(len(x[3]["spec"]), foff[x[1]], np.int64) for x in secs])
    ng = np.ones(n, np.uint16)
    np.add.at(ng, sec["spec"], 1)
    mf = {k: np.concatenate([f[k] for f in F]) for k in F[0]}
    return {"prim": prim, "sec": sec, "n_groups": ng}, mf


def g2():
    import pyarrow.parquet as pq
    import torch

    import gates as M3G
    import m3_common as C
    import m3_data as MD
    from vulcan.chimeric import search as CS
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    from vulcan.rescore.features_v2 import model_psm
    from vulcan.stats import g2 as G2
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import psm_competition
    mixes, sps, stores, names = mixtures()
    idx_dir = C.index_dir("human_contam", "prod")
    lab = MD.labels(idx_dir)
    l1 = L1Index(idx_dir)
    p = E.SearchParams(report_psms=5)
    gi = E.GpuIndex(l1, p, STATIC, VAR)
    fb = scorer(os.path.join(W, "predcache_g2"))
    ms = models()
    rec = {k: [] for k in ("wf", "sf", "wl", "sl", "w1", "s1", "prov", "grp")}
    excl, info = {}, {}
    os.makedirs(os.path.join(W, "g2"), exist_ok=True)
    for m in mixes:
        mix = m["mix"]
        groups, mf = mixture_groups(m, sps, stores)
        t0 = time.perf_counter()
        tbl, st = CS.search_chimeric(mix, gi, l1, p, groups)
        pq.write_table(tbl, os.path.join(W, "g2", mix.run + ".chim.psms.parquet"))
        psm = {c: tbl.column(c).to_numpy(zero_copy_only=False) for c in tbl.column_names}
        feat, finfo = fb.run_v2(l1, psm, mix, mf)
        mp = model_psm(psm)
        sc = {k: mm.apply(mp, feat) for k, mm in ms.items()}
        info[mix.run] = {"n_spectra": len(mix), "n_psms": tbl.num_rows, "n_secondary_groups": st["n_secondary_groups"],
                         "nce": finfo["nce"], "wall_s": time.perf_counter() - t0}
        rows, spec, grp = psm["pf_row"], psm["spec_idx"].astype(np.int64), psm["group"].astype(np.int64)
        dec = lab["is_decoy"][rows]
        tseq = lab["target_seq"][rows]
        pv = m["prov"][spec]
        bad_spec = np.zeros(len(mix), bool)
        excl[mix.run] = {}
        for src in (2, 3):
            spn = pv == src
            bad = M3G.homolog_excluded_seqs(set(tseq[spn].tolist()), src)
            hit = spn & np.isin(tseq, list(bad))
            bad_spec[spec[hit]] = True
            excl[mix.run][f"source{src}"] = {"spectra_excluded": int(np.unique(spec[hit]).size), "homolog_peptides": len(bad)}
        null_row = (pv > 0) & ~bad_spec[spec]
        unit = spec * 64 + grp
        tb = hash_tiebreak(psm["pf_id"], SEED)
        hs = psm["hyperscore"].astype(float)
        w0 = psm_competition(unit, hs, tiebreak=tb)
        w0 = w0[null_row[w0]]
        ou = unit[w0]

        def aligned(score):
            w = psm_competition(unit, score, tiebreak=tb)
            w = w[null_row[w]]
            o = np.argsort(unit[w])
            ww = w[o][np.searchsorted(unit[w][o], ou)]
            assert np.array_equal(unit[ww], ou)
            return ww
        wl = aligned(sc["v2"])
        w1 = aligned(sc["v1"])
        rec["wf"].append(~dec[w0]); rec["sf"].append(hs[w0])
        rec["wl"].append(~dec[wl]); rec["sl"].append(sc["v2"][wl])
        rec["w1"].append(~dec[w1]); rec["s1"].append(sc["v1"][w1])
        rec["prov"].append(m["prov"][spec[w0]]); rec["grp"].append(grp[w0])
        print(mix.run, info[mix.run], flush=True)
    del gi
    torch.cuda.empty_cache()
    cat = {k: np.concatenate(v) for k, v in rec.items()}
    wf, sf, prov, grp = cat["wf"], cat["sf"], cat["prov"], cat["grp"]
    out = {"n_null_competitions": int(wf.size), "n_offset": int((prov == 1).sum()), "n_species_pfu": int((prov == 2).sum()),
           "n_species_yeast": int((prov == 3).sum()), "n_secondary_units": int((grp > 0).sum()),
           "homolog_filter": excl, "search": info, "p_target_win_fixed": float(wf.mean()),
           "budget_g2_sample_size": G2.g2_sample_size(), "scorers": {}}
    for name, wl, sl in (("v2", cat["wl"], cat["sl"]), ("v1_chimeric", cat["w1"], cat["s1"])):
        res = G2.g2_differential(wf, wl, sf)
        by = {}
        for lab_, m_ in (("offset", prov == 1), ("pfu", prov == 2), ("yeast", prov == 3), ("primary_units", grp == 0),
                         ("secondary_units", grp > 0)):
            if m_.sum() > 100:
                r_ = G2.g2_differential(wf[m_], wl[m_], sf[m_])
                by[lab_] = {"pass": r_["pass"], "diff": r_["strata"][0]["overall"]["diff"],
                            "ucb": r_["strata"][0]["overall"]["ucb"], "top_ucb": r_["strata"][0]["top"]["ucb"],
                            "n": int(m_.sum())}
        out["scorers"][name] = {"g2": res, "g2_by_subset": by, "g5": G2.leak_auc(sl, wl),
                                "g5_fixed_reference": G2.leak_auc(sf, wf), "g2_absolute_diag": G2.g2_absolute(wl, sl),
                                "frac_winner_changed": float(np.mean(wl != wf))}
        print(name, "G2", res["pass"], round(res["strata"][0]["overall"]["diff"], 4),
              round(res["strata"][0]["overall"]["ucb"], 4), "top", round(res["strata"][0]["top"]["ucb"], 4),
              "G5", out["scorers"][name]["g5"], flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--g0", action="store_true")
    ap.add_argument("--g2", action="store_true")
    ap.add_argument("--out", default=os.path.join(RES, "gates_v2.json"))
    a = ap.parse_args()
    out = json.load(open(a.out)) if os.path.exists(a.out) else {}
    if a.g0:
        out["G0"] = g0()
        print("G0", out["G0"]["pass"], flush=True)
        json.dump(out, open(a.out, "w"), indent=1, default=float)
    if a.g2:
        out["G2_original_G5"] = g2()
    out["pending"] = {"tail_G2": "pending A2 DN2 nulls", "G1": "pending A2 DN2 nulls", "G3": "pending A2 DN2 nulls"}
    json.dump(out, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
