#!/usr/bin/env python3
"""3c gates for the frozen IT-CID static scorers that do NOT depend on the pending A2 null set (DESIGN §6.4):

G0  label blindness: one D10cf fraction (label-free IT-CID) searched against the E2 Pfu E2001 index and E2's swap
    index (FASTA order shuffled, pair list permuted, target/decoy roles swapped), and one D4 fraction (TMT) against
    TMT-static versions of the same two indexes (built here with the same seeds); top-5 PSM tables, every L3 feature
    (Prosit via Koina, a FRESH prediction cache per side) and every frozen IT-CID score must be bit-identical.
G2  original PSM-level differential null test (M3 implementation, vulcan.stats.g2.g2_differential, class IT-CID):
    mixture design: every D10cf fraction (IT-CID, label-free human host) + K = 2 offset copies of its known-charge
    spectra (n + 0.5 + U(-0.1, 0.1) Da, n in {-2,-1,1,2}, seed 11) + an equal share (seeded permutation) of the yeast
    IT-CID null spectra (set DN3c, data/manifests/null_itcid.tsv), searched with the production IT-CID settings
    (human prod index, top-5), L3 features on the mixture.  Homolog filter for species nulls as M3 (I/L-identical,
    isobaric, single-substitution neighbours of yeast + contaminant tryptic peptides; label-symmetric).
G5  leak AUC on the G2 null winners (vulcan.gates.stats.g5_gate, alarm 0.55; also vulcan.stats.g2.leak_auc).
G4  decoy score distributions across charge strata {1, 2, 3, 4+} on TRAINING decoys (13 IT-CID training runs),
    vulcan.gates.stats.g4_gate (tau = 95th pct of within-stratum between-run KS), as M3e.

    python scripts/m3c/gates_itcid.py g0 | g2prep | g2 | g4
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m3"))
sys.path.insert(0, HERE)
import common as C  # noqa: E402
import m3_data as MD  # noqa: E402

import numpy as np  # noqa: E402

from vulcan.rescore.static import StaticModel  # noqa: E402

STATIC = ["itcid_linear_base", "itcid_gbm_base", "itcid_linear_pred", "itcid_gbm_pred"]
SEED = 11
K_COPIES = 2
MIX_DIR = os.path.join(C.WORK, "g2mix")
OUT = os.path.join(C.RES, "gates.json")


def models():
    return {s: StaticModel.load(os.path.join(C.RES, "models", s + ".pkl")) for s in STATIC}


def psm_dict(tbl):
    return {c: tbl.column(c).to_numpy(zero_copy_only=False) for c in MD.PSM_COLS}


# ------------------------------------------------------------------------------------------------ G0
def g0(shared=False):
    import torch
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index, build
    from vulcan.predict.cache import PredCache
    from vulcan.predict.koina import KoinaProsit
    from vulcan.rescore.itcid import ItcidFeatureBuilder
    ms = models()
    res_all = {}
    for tag, store, tmt in (("LF_D10cf", C.d10_stores()[10], False), ("TMT_D4", C.d4_stores()[0], True)):
        if tmt:
            base = os.path.join(C.WORK, "index", "g0_pfu_tmt_E2001")
            alt = os.path.join(C.WORK, "index", "g0_pfu_tmt_swapped")
            par = {"static_mods": dict(C.CAM, K=C.TMT), "entrap_seeds": [2001]}
            if not os.path.exists(os.path.join(base, "meta.json")):
                build(["/tmp/vulcan_e2/pfu_contam.fasta"], base, par)
            if not os.path.exists(os.path.join(alt, "meta.json")):
                build(["/tmp/vulcan_e2/g0/shuffled.fasta"], alt, par, shuffle_rows_seed=23, swap_labels=True)
            st, var = C.mods("human_tmt")
        else:
            base = "/tmp/vulcan_e2/index/pfu_contam_E2001"
            alt = "/tmp/vulcan_e2/g0/index_swapped"
            st, var = C.mods("human")
        p = E.SearchParams(**C.SEARCH)
        sp = E.load_store(store)
        outs = []
        for idx in (base, alt):
            l1 = L1Index(idx)
            if tmt:
                l1.mass32 = (l1.mass32.astype(np.float32) + np.float32(C.TMT)).astype(np.float32)
                l1.mask = l1.mask | np.uint64(1 << 63)
            gi = E.GpuIndex(l1, p, st, var)
            tbl, _, _ = E.search(sp, gi, l1, p)
            del gi
            torch.cuda.empty_cache()
            psm = psm_dict(tbl)
            root = os.path.join(C.WORK, "g0_cache", tag + "_" + ("shared" if shared else os.path.basename(idx)))
            pred = KoinaProsit("tmt" if tmt else "lf", st, var)
            caches = {}

            def fac(nce, root=root, pred=pred, caches=caches):
                if nce not in caches:
                    caches[nce] = PredCache(root, pred, nce, "IT", st, var)
                return caches[nce]
            feat, info = ItcidFeatureBuilder(fac, st, var).run(l1, psm, sp)
            scores = {s: m.apply(psm, feat) for s, m in ms.items()}
            outs.append((tbl, feat, scores))
        (t0, f0, s0), (t1, f1, s1) = outs
        r = {"store": store, "n_psms": t0.num_rows, "psm_tables_bit_identical": bool(t0.equals(t1)),
             "features_bit_identical": {k: bool(np.array_equal(f0[k], f1[k], equal_nan=True)) for k in f0},
             "features_max_abs_diff": {k: float(np.nanmax(np.abs(f0[k].astype(float) - f1[k].astype(float))))
                                       if f0[k].size else 0.0 for k in f0},
             "static_scores_bit_identical": {k: bool(np.array_equal(s0[k], s1[k])) for k in s0}}
        r["pass"] = r["psm_tables_bit_identical"] and all(r["features_bit_identical"].values()) and \
            all(r["static_scores_bit_identical"].values())
        res_all[tag] = r
        print(tag, "G0", r["pass"], r["psm_tables_bit_identical"], flush=True)
    res_all["pass"] = all(v["pass"] for k, v in res_all.items() if isinstance(v, dict))
    return res_all


# ------------------------------------------------------------------------------------------------ G2 mixtures
def g2prep():
    import pyarrow.parquet as pq
    sys.path.insert(0, os.path.join(HERE, "..", "m3"))
    import g2_nulls as M3N
    from vulcan.gpu import engine as E
    import csv
    os.makedirs(MIX_DIR, exist_ok=True)
    rng = np.random.default_rng(SEED)
    hosts = [E.load_store(s) for s in C.d10_stores()]
    nulls = [r for r in csv.DictReader(open(os.path.join(C.ROOT, "data", "manifests", "null_itcid.tsv")), delimiter="\t")]
    ys = [E.load_store(os.path.join(C.WORK, "store_null", os.path.splitext(r["filename"])[0] + ".vspec.arrow"))
          for r in nulls]
    foreign = M3N.concat_spectra(ys, "foreign")
    parts_idx = np.array_split(rng.permutation(len(foreign)), len(hosts))
    mixes, info = [], {}
    for sp, half in zip(hosts, parts_idx):
        parts, prov = [sp], [np.zeros(len(sp), np.int8)]
        known = np.flatnonzero(sp.charge > 0)
        for _ in range(K_COPIES):
            n = rng.choice([-2, -1, 1, 2], size=len(known))
            delta = n + 0.5 + rng.uniform(-0.1, 0.1, size=len(known))
            parts.append(M3N.subset(sp, known, sp.precursor_mz[known] + delta / sp.charge[known].astype(np.float64)))
            prov.append(np.ones(len(known), np.int8))
        half = np.sort(half)
        parts.append(M3N.subset(foreign, half))
        prov.append(np.full(len(half), 3, np.int8))
        mix = M3N.concat_spectra(parts, sp.run + "_g2mix")
        mixes.append(mix)
        pv = np.concatenate(prov)
        pickle.dump(mix, open(os.path.join(MIX_DIR, mix.run + ".spectra.pkl"), "wb"), protocol=5)
        np.save(os.path.join(MIX_DIR, mix.run + ".provenance.npy"), pv)
        info[mix.run] = {"n_spectra": len(mix), "n_real": int((pv == 0).sum()), "n_offset": int((pv == 1).sum()),
                         "n_yeast": int((pv == 3).sum())}
    p = E.SearchParams(**C.SEARCH)
    l1 = C.load_l1("human", "prod")
    st, var = C.mods("human")
    gi = E.GpuIndex(l1, p, st, var)
    out = C.psm_dir("G2_mix")
    os.makedirs(out, exist_ok=True)
    for mix in mixes:
        tbl, _, _ = E.search(mix, gi, l1, p)
        pq.write_table(tbl, os.path.join(out, mix.run + ".psms.parquet"))
        info[mix.run]["n_psms"] = tbl.num_rows
    del gi
    # L3 features on the mixtures (the run-level calibration sees the mixture, as in production)
    import pyarrow as pa
    import importlib.util
    spec = importlib.util.spec_from_file_location("m3c_features", os.path.join(HERE, "features.py"))
    F = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(F)
    fb, _ = F.builder("human")
    for mix in mixes:
        t = pq.read_table(os.path.join(out, mix.run + ".psms.parquet"))
        psm = {k: t.column(k).to_numpy() for k in ("pf_row", "charge", "spec_idx", "rank", "hyperscore", "rt")}
        feats, fi = fb.run(l1, psm, mix)
        pq.write_table(pa.table(feats), os.path.join(out, mix.run + ".feat.parquet"))
        info[mix.run]["features"] = {k: fi[k] for k in ("rt_mad_min", "frac_pred_missing")}
        print(mix.run, info[mix.run], flush=True)
    json.dump({"k_copies": K_COPIES, "seed": SEED, "null_runs": [r["filename"] for r in nulls], "runs": info},
              open(os.path.join(out, "g2_info.json"), "w"), indent=1, default=float)


def g2():
    sys.path.insert(0, os.path.join(HERE, "..", "m3"))
    import gates as M3G                     # M3 homolog filter (yeast + contaminants + SNCA = source 3)
    from vulcan.gates import stats as GS
    from vulcan.stats import g2 as G2
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import psm_competition
    ms = models()
    lab = MD.labels(C.index_dir("human", "prod"))
    d = C.psm_dir("G2_mix")
    runs = MD.load_runs(d)
    rec = {"win_fixed": [], "score_fixed": [], "prov": [], "win": {s: [] for s in ms}, "score": {s: [] for s in ms}}
    n_excl = 0
    for run, psm, feat in runs:
        prov = np.load(os.path.join(MIX_DIR, run + ".provenance.npy"))
        rows = psm["pf_row"]
        dec = lab["is_decoy"][rows]
        tseq = lab["target_seq"][rows]
        spec = psm["spec_idx"]
        tb = hash_tiebreak(psm["pf_id"], SEED)
        pv = prov[spec]
        bad_spec = np.zeros(prov.size, bool)
        sp_null = pv == 3
        bad = M3G.homolog_excluded_seqs(set(tseq[sp_null].tolist()), 3)
        hit = sp_null & np.isin(tseq, list(bad))
        bad_spec[spec[hit]] = True
        n_excl += int(np.unique(spec[hit]).size)
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
            o = np.argsort(spec[w])
            ww = w[o][np.searchsorted(spec[w][o], order_spec)]
            assert np.array_equal(spec[ww], order_spec)
            rec["win"][s].append(~dec[ww])
            rec["score"][s].append(sc[ww])
    cat = np.concatenate
    wf, sf, prov = cat(rec["win_fixed"]), cat(rec["score_fixed"]), cat(rec["prov"])
    out = {"class": "IT-CID", "n_null_competitions": int(wf.size), "n_offset": int((prov == 1).sum()),
           "n_species_yeast": int((prov == 3).sum()), "homolog_filter_spectra_excluded": n_excl,
           "p_target_win_fixed": float(wf.mean()), "budget_g2_sample_size": G2.g2_sample_size(), "scorers": {}}
    for s in ms:
        wl, sl = cat(rec["win"][s]), cat(rec["score"][s])
        res = G2.g2_differential(wf, wl, sf)
        by_src = {src: G2.g2_differential(wf[prov == k], wl[prov == k], sf[prov == k])["pass"]
                  for src, k in (("offset", 1), ("yeast", 3)) if (prov == k).sum() > 100}
        out["scorers"][s] = {"g2": res, "g2_pass_by_source": by_src, "g5_leak_auc": G2.leak_auc(sl, wl),
                             "g5_gate": GS.g5_gate(sl, wl), "g5_fixed_reference": G2.leak_auc(sf, wf),
                             "frac_winner_changed": float(np.mean(wl != wf))}
        print(s, "G2", res["pass"], json.dumps(res["strata"][0], default=float)[:300], "G5",
              out["scorers"][s]["g5_gate"].get("pass"), flush=True)
    return out


# ------------------------------------------------------------------------------------------------ G4
def g4():
    from vulcan.gates import stats as GS
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import psm_competition
    ms = models()
    acc = {p: {"score": [], "stratum": [], "run": []} for p in ["hyperscore"] + list(ms)}
    import pyarrow.parquet as pq
    for r in C.train_rows():
        run = os.path.splitext(r["filename"])[0]
        lab = MD.labels(C.index_dir(C.train_family(r), "prod"))
        d = C.psm_dir("train")
        psm = psm_dict(pq.read_table(os.path.join(d, run + ".psms.parquet")))
        ft = pq.read_table(os.path.join(d, run + ".feat.parquet"))
        feat = {c: ft.column(c).to_numpy() for c in ft.column_names}
        dec = lab["is_decoy"][psm["pf_row"]]
        tb = hash_tiebreak(psm["pf_id"], SEED)
        scores = {"hyperscore": psm["hyperscore"].astype(float), **{s: m.apply(psm, feat) for s, m in ms.items()}}
        for p, s in scores.items():
            w = psm_competition(psm["spec_idx"], s, tiebreak=tb)
            w = w[dec[w]]
            acc[p]["score"].append(s[w])
            acc[p]["stratum"].append(np.minimum(psm["charge"][w].astype(int), 4))
            acc[p]["run"].append(np.full(w.size, run))
    out = {"partition": "IT-CID x CID x precursor charge {1, 2, 3, 4+} (as M3e; frozen before estimation)",
           "runs": "the 13 IT-CID training runs (train_itcid.tsv)"}
    for p, d in acc.items():
        out[p] = GS.g4_gate(np.concatenate(d["score"]), np.concatenate(d["stratum"]), np.concatenate(d["run"]))
        print(p, out[p]["pass"], round(out[p]["tau"], 4), flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["g0", "g0shared", "g2prep", "g2", "g4"])
    a = ap.parse_args()
    if a.what == "g2prep":
        return g2prep()
    out = json.load(open(OUT)) if os.path.exists(OUT) else {}
    if a.what == "g0shared":
        # G0 variant: one content-keyed (pf_id) prediction cache for both sides, so the remote predictor is queried
        # once per peptidoform; isolates our pipeline from Koina's run-to-run float nondeterminism.
        t0 = time.perf_counter()
        out["G0_shared_prediction_cache"] = g0(shared=True)
        out["G0_shared_prediction_cache"]["wall_s"] = time.perf_counter() - t0
        json.dump(out, open(OUT, "w"), indent=1, default=float)
        return
    t0 = time.perf_counter()
    out[{"g0": "G0", "g2": "G2_G5", "g4": "G4"}[a.what]] = {"wall_s": None, **globals()[a.what]()}
    out[{"g0": "G0", "g2": "G2_G5", "g4": "G4"}[a.what]]["wall_s"] = time.perf_counter() - t0
    json.dump(out, open(OUT, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
