#!/usr/bin/env python3
"""E5: bi-encoder retrieval (the split stage of the frozen ce_s model) vs the K1 inverted index (DESIGN §10 E5).

A. Open-search recall on D0 (M4 open mode: human_contam_prod, delta [-150, +500] Da, K1 top-K 500, top-5 reported).
   Reference = rank-1 TARGET PSMs accepted at PSM-level TDC q <= 1% (hyperscore), split narrow (|mass_shift| < 0.5)
   / shifted.  Bi-encoder: exact brute-force cosine over every peptidoform (targets and decoys) whose mass lies in
   the spectrum's open window; rank of the reference peptidoform among them -> recall@K, K in {10,100,500,1000}.
   Peptides embedded with charge token 2 and 3 (two banks; spectra of charge >= 3 query bank 3).
B. G3-style retrieval symmetry on the existing dev nulls (M3e G2 mixtures 00/01, null-provenance spectra; as
   scripts/m3e/g3.py): exhaustive closed-window K2 rescoring -> best target / best decoy hyperscore per spectrum;
   pruned = bi-encoder top-K (label-blind, by cosine, ties by row) within the same window, then the best hyperscore
   per label among the kept candidates.  vulcan.gates.stats.g3_gate (|recall_T - recall_D| <= 0.01 per decile),
   K in {10, 100}; all nulls and the A1.3-clean species subset.
    python scripts/m7/e5_retrieval.py --open --g3 --model results/2026-10-28_M7_encoders/models/ce_s.pt
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))
sys.path.insert(0, os.path.join(ROOT, "scripts", "m2"))
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3e"))
import m3_common as C  # noqa: E402
import m3_data as MD  # noqa: E402

import numpy as np  # noqa: E402
import torch  # noqa: E402

W = os.environ.get("VULCAN_M7_WORK", "/tmp/vulcan_m7")
RES = os.path.join(ROOT, "results", "2026-10-28_M7_encoders")
D0 = ["/tmp/vulcan_e2/store/b1906_293T_proteinID_01A_QE3_122212.vspec.arrow",
      "/tmp/vulcan_e2/store/b1937_293T_proteinID_01B_QE3_122212.vspec.arrow"]
KS = (10, 100, 500, 1000)
PROTON = 1.00727646688


def peptide_banks(model, l1):
    from vulcan.encoders.score import embed_peptides
    rows = np.arange(len(l1), dtype=np.int64)
    t0 = time.perf_counter()
    banks = {z: embed_peptides(model, l1, rows, charge=z) for z in (2, 3)}
    torch.cuda.synchronize()
    return banks, time.perf_counter() - t0


def open_search(l1):
    import pyarrow.parquet as pq
    from vulcan.gpu import engine as E
    from vulcan.ptm import search as PS
    out = os.path.join(W, "e5", "open")
    os.makedirs(out, exist_ok=True)
    p = E.SearchParams()
    mp = PS.ModeParams()
    mp.open_da = (-150.0, 500.0)
    stats = {}
    gi = None
    for s in D0:
        sp = E.load_store(s)
        f = os.path.join(out, sp.run + ".psms.parquet")
        if os.path.exists(f):
            continue
        gi = gi or E.GpuIndex(l1, p, C.STATIC, C.VAR)
        tbl, st = PS.search_modes(sp, gi, l1, p, mp, report=5)
        pq.write_table(tbl, f)
        stats[sp.run] = {"wall_s": st.wall_s, "n_spectra": len(sp)}
        print("open", sp.run, tbl.num_rows, st.wall_s.get("total"), flush=True)
    if stats:
        json.dump(stats, open(os.path.join(out, "stats.json"), "w"), indent=1, default=float)
    del gi
    torch.cuda.empty_cache()
    return out


def part_open(model):
    import pyarrow.parquet as pq
    from vulcan.encoders.score import embed_spectra
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    from vulcan.stats.tdc import tdc_qvalues
    idx = C.index_dir("human_contam", "prod")
    l1 = L1Index(idx)
    lab = MD.labels(idx)
    d = open_search(l1)
    banks, t_bank = peptide_banks(model, l1)
    mass = torch.from_numpy(l1.mass32.astype(np.float32)).cuda()
    res = {"n_peptidoforms": len(l1), "peptide_embedding_s_two_banks": t_bank, "per_run": {}}
    ranks_all, shifted_all = [], []
    t_query = 0.0
    for s in D0:
        sp = E.load_store(s)
        t = pq.read_table(os.path.join(d, sp.run + ".psms.parquet"))
        rk = t.column("rank").to_numpy()
        m1 = rk == 1
        rows = t.column("pf_row").to_numpy()[m1]
        spec = t.column("spec_idx").to_numpy()[m1]
        hs = t.column("hyperscore").to_numpy()[m1]
        z = t.column("charge").to_numpy()[m1].astype(np.int64)
        expm = t.column("expmass").to_numpy()[m1].astype(np.float64)
        shift = t.column("mass_shift").to_numpy()[m1]
        dec = lab["is_decoy"][rows]
        q = tdc_qvalues(hs, dec)
        acc = (q <= 0.01) & ~dec
        rows, spec, z, expm, shift = rows[acc], spec[acc], z[acc], expm[acc], shift[acc]
        zq = np.where(sp.charge[spec] > 0, sp.charge[spec], z).astype(np.int64)
        semb = embed_spectra(model, sp, np.where(sp.charge > 0, sp.charge, 2))[torch.from_numpy(spec).cuda()]
        bank_id = torch.from_numpy((zq >= 3).astype(np.int64)).cuda()
        lo = torch.from_numpy((expm - 500.0).astype(np.float32)).cuda()
        hi = torch.from_numpy((expm + 150.0).astype(np.float32)).cuda()
        ref = torch.empty(len(rows), device="cuda")
        rr = torch.from_numpy(rows).cuda()
        for zb, z_ in ((0, 2), (1, 3)):
            m = bank_id == zb
            ref[m] = (semb[m].float() * banks[z_][rr[m]].float()).sum(1)
        rank = torch.zeros(len(rows), dtype=torch.int64, device="cuda")
        torch.cuda.synchronize()
        tq = time.perf_counter()
        CH = 1 << 19
        for zb, z_ in ((0, 2), (1, 3)):
            sel = torch.nonzero(bank_id == zb).squeeze(1)
            for b in range(0, len(sel), 2048):
                si = sel[b:b + 2048]
                q_ = semb[si]
                for c0 in range(0, len(l1), CH):
                    sc = q_ @ banks[z_][c0:c0 + CH].T
                    pm = mass[c0:c0 + CH]
                    inwin = (pm[None] >= lo[si, None]) & (pm[None] <= hi[si, None])
                    rank[si] += ((sc.float() > ref[si, None]) & inwin).sum(1)
        torch.cuda.synchronize()
        t_query += time.perf_counter() - tq
        rk_ = rank.cpu().numpy()
        ranks_all.append(rk_)
        shifted_all.append(np.abs(shift) >= 0.5)
        res["per_run"][sp.run] = {"n_reference": int(len(rows)), "n_shifted": int((np.abs(shift) >= 0.5).sum()),
                                  "recall": {str(k): float((rk_ < k).mean()) for k in KS}}
        print(sp.run, res["per_run"][sp.run], flush=True)
    R = np.concatenate(ranks_all)
    S = np.concatenate(shifted_all)
    res["recall_all"] = {str(k): float((R < k).mean()) for k in KS}
    res["recall_narrow"] = {str(k): float((R[~S] < k).mean()) for k in KS}
    res["recall_shifted"] = {str(k): float((R[S] < k).mean()) for k in KS}
    res["median_rank"] = {"all": float(np.median(R)), "narrow": float(np.median(R[~S])), "shifted": float(np.median(R[S]))}
    res["n_reference"] = int(len(R))
    res["n_shifted"] = int(S.sum())
    res["bruteforce_query_s_reference_spectra"] = t_query
    res["bruteforce_query_s_per_spectrum"] = t_query / max(len(R), 1)
    st = json.load(open(os.path.join(d, "stats.json"))) if os.path.exists(os.path.join(d, "stats.json")) else {}
    res["k1_open_stats"] = st
    print(json.dumps({k: v for k, v in res.items() if k != "per_run"}, default=float), flush=True)
    return res


def part_g3(model, Ks=(10, 100)):
    from e2_correctness import _best_per_label, _rescored
    from gates_a1 import homolog_filters
    from vulcan.encoders.score import embed_spectra
    from vulcan.gates import stats as GS
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    from vulcan.encoders.data import peptide_tokens
    idx = C.index_dir("human_contam", "prod")
    l1 = L1Index(idx)
    lab = MD.labels(idx)
    isdec = lab["is_decoy"]
    p0 = E.SearchParams()
    gi = E.GpuIndex(l1, p0, C.STATIC, C.VAR)
    hfs = homolog_filters()
    MIX = "/tmp/vulcan_m3e/null_mix"
    acc = {k: [] for k in ["bt_x", "bd_x", "prov", "clean"] + [f"bt_{K}" for K in Ks] + [f"bd_{K}" for K in Ks]}
    for name in ("g2mix00_h0", "g2mix01_h1"):
        sp = pickle.load(open(os.path.join(MIX, name + ".spectra.pkl"), "rb"))
        prov = np.load(os.path.join(MIX, name + ".prov.npy"))
        n = len(sp)
        nul = np.flatnonzero(prov >= 1)
        _, _, dx = E.search_exhaustive(sp, gi, l1, p0, spec_subset=nul)
        sx, hx, rx = _rescored(dx, l1, p0)
        btx, bdx = _best_per_label(sx, hx, isdec[rx], n)
        # bi-encoder cosine of every exhaustive candidate (label-blind)
        semb = embed_spectra(model, sp, np.where(sp.charge > 0, sp.charge, 2))
        cos = np.zeros(len(sx), np.float32)
        zs = np.where(sp.charge[sx] >= 3, 3, 2)
        with torch.inference_mode():
            for b in range(0, len(sx), 1 << 15):
                sl = slice(b, b + (1 << 15))
                tok = peptide_tokens(l1, rx[sl])
                L = int((tok > 0).sum(1).max())
                for z_ in (2, 3):
                    m = zs[sl] == z_
                    if not m.any():
                        continue
                    t = torch.from_numpy(tok[m][:, :L].astype(np.int64)).cuda()
                    _, _, pe, _ = model.encode_peptides(t, torch.full((len(t),), z_, device="cuda"))
                    se = semb[torch.from_numpy(sx[sl][m]).cuda()].float()
                    cc = np.zeros(m.sum(), np.float32)
                    cc[:] = (se * pe.float()).sum(1).cpu().numpy()
                    idxs = np.flatnonzero(m) + b
                    cos[idxs] = cc
        o = np.lexsort((rx, -cos, sx))
        s_o = sx[o]
        first = np.searchsorted(s_o, s_o, "left")
        rnk = np.empty(len(sx), np.int64)
        rnk[o] = np.arange(len(sx)) - first
        for K in Ks:
            keep = rnk < K
            btk, bdk = _best_per_label(sx[keep], hx[keep], isdec[rx][keep], n)
            acc[f"bt_{K}"].append(btk[nul]); acc[f"bd_{K}"].append(bdk[nul])
        clean = np.ones(n, bool)
        for lab_dec, best in ((False, btx), (True, bdx)):
            m = isdec[rx] == lab_dec
            hit = m & (hx == best[sx])
            brow = np.full(n, -1, np.int64)
            brow[sx[hit]] = rx[hit]
            for src, hf in hfs.items():
                sel = np.flatnonzero((prov == src) & (brow >= 0))
                if sel.size == 0:
                    continue
                rr = np.unique(np.r_[brow[sel], l1.partner[brow[sel]]])
                rr = rr[rr >= 0]
                bad = rr[hf.near([l1.sequence(int(r)) for r in rr])]
                clean[sel[np.isin(brow[sel], bad) | np.isin(l1.partner[brow[sel]], bad)]] = False
        for k, v in (("bt_x", btx), ("bd_x", bdx), ("prov", prov), ("clean", clean)):
            acc[k].append(v[nul])
        cnt = np.bincount(sx, minlength=n)[nul]
        print(name, "nulls", nul.size, "median window candidates", float(np.median(cnt[cnt > 0])), flush=True)
        torch.cuda.empty_cache()
    A = {k: np.concatenate(v) for k, v in acc.items()}
    out = {"mixtures": ["g2mix00_h0", "g2mix01_h1"], "by_K": {}}
    cs = (A["prov"] >= 2) & A["clean"]
    for K in Ks:
        allg = GS.g3_gate(A["bt_x"], A["bd_x"], A[f"bt_{K}"], A[f"bd_{K}"])
        clg = GS.g3_gate(*(np.where(cs, A[k], -np.inf) for k in ("bt_x", "bd_x", f"bt_{K}", f"bd_{K}")))
        out["by_K"][str(K)] = {"all_nulls": allg, "clean_species_nulls": clg, "pass": allg["pass"] and clg["pass"]}
        print("G3-style K", K, allg["pass"], clg["pass"],
              [(round(b["recall_target"], 3), round(b["recall_decoy"], 3)) for b in allg["bins"]], flush=True)
    del gi
    torch.cuda.empty_cache()
    return out


def part_open_null_symmetry(model, Ks=(10, 100, 1000)):
    """ADDED AFTER the pre-stated G3-style result (diagnostic, not a pre-stated criterion): label symmetry of
    open-window bi-encoder retrieval on the same null spectra (mixtures 00/01, null provenance; all and A1.3-clean
    species nulls are not separable here without per-hit homolog checks, so all null spectra are used).
    Per spectrum: top-K by cosine over every peptidoform (targets + decoys) with mass in [M - 500, M + 150];
    target fraction among the K hits.  Pseudo-reverse decoys are mass-identical to their targets, so a
    label-symmetric retriever gives 0.5 in expectation."""
    from vulcan.encoders.score import embed_spectra
    from vulcan.index.build import L1Index
    idx = C.index_dir("human_contam", "prod")
    l1 = L1Index(idx)
    lab = MD.labels(idx)
    isdec = torch.from_numpy(lab["is_decoy"].astype(np.bool_)).cuda()
    banks, _ = peptide_banks(model, l1)
    mass = torch.from_numpy(l1.mass32.astype(np.float32)).cuda()
    MIX = "/tmp/vulcan_m3e/null_mix"
    Kmax = max(Ks)
    res = {"note": part_open_null_symmetry.__doc__, "by_K": {}}
    fr = {K: [] for K in Ks}
    provs = []
    for name in ("g2mix00_h0", "g2mix01_h1"):
        sp = pickle.load(open(os.path.join(MIX, name + ".spectra.pkl"), "rb"))
        prov = np.load(os.path.join(MIX, name + ".prov.npy"))
        nul = np.flatnonzero((prov >= 1) & (sp.charge > 0) & np.isfinite(sp.precursor_mz))
        z = sp.charge[nul].astype(np.int64)
        M = (sp.precursor_mz[nul] - PROTON) * z
        semb = embed_spectra(model, sp, np.where(sp.charge > 0, sp.charge, 2))[torch.from_numpy(nul).cuda()]
        lo = torch.from_numpy((M - 500).astype(np.float32)).cuda()
        hi = torch.from_numpy((M + 150).astype(np.float32)).cuda()
        bank_id = torch.from_numpy((z >= 3).astype(np.int64)).cuda()
        CH = 1 << 19
        tf = {K: torch.zeros(len(nul), device="cuda") for K in Ks}
        with torch.inference_mode():
            for zb, z_ in ((0, 2), (1, 3)):
                sel = torch.nonzero(bank_id == zb).squeeze(1)
                for b in range(0, len(sel), 1024):
                    si = sel[b:b + 1024]
                    q_ = semb[si]
                    bv = torch.full((len(si), Kmax), -1e4, device="cuda")
                    bi = torch.zeros((len(si), Kmax), dtype=torch.int64, device="cuda")
                    for c0 in range(0, len(l1), CH):
                        sc = (q_ @ banks[z_][c0:c0 + CH].T).float()
                        pm = mass[c0:c0 + CH]
                        sc.masked_fill_(~((pm[None] >= lo[si, None]) & (pm[None] <= hi[si, None])), -1e4)
                        v, i = torch.topk(sc, Kmax, dim=1)
                        cv = torch.cat([bv, v], 1)
                        ci = torch.cat([bi, i + c0], 1)
                        o = torch.topk(cv, Kmax, dim=1).indices
                        bv, bi = torch.gather(cv, 1, o), torch.gather(ci, 1, o)
                    ok = bv > -1e3
                    for K in Ks:
                        t = (~isdec[bi[:, :K]]) & ok[:, :K]
                        tf[K][si] = t.sum(1).float() / ok[:, :K].sum(1).clamp_min(1).float()
        for K in Ks:
            fr[K].append(tf[K].cpu().numpy())
        provs.append(prov[nul])
        print(name, "null spectra", len(nul), {K: float(tf[K].mean()) for K in Ks}, flush=True)
    P = np.concatenate(provs)
    for K in Ks:
        f = np.concatenate(fr[K])
        se = f.std(ddof=1) / np.sqrt(len(f))
        res["by_K"][str(K)] = {"mean_target_fraction": float(f.mean()), "se": float(se), "n_spectra": int(len(f)),
                               "by_source": {str(int(u)): float(f[P == u].mean()) for u in np.unique(P)}}
    print(json.dumps(res["by_K"]), flush=True)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--open", action="store_true")
    ap.add_argument("--g3", action="store_true")
    ap.add_argument("--open-null", action="store_true")
    ap.add_argument("--model", default=os.path.join(RES, "models", "ce_s.pt"))
    ap.add_argument("--out", default=os.path.join(RES, "e5.json"))
    a = ap.parse_args()
    from vulcan.encoders.score import load_model
    model = load_model(a.model)
    out = json.load(open(a.out)) if os.path.exists(a.out) else {}
    out["model_sha256"] = model.sha256
    if a.open:
        out["open_recall"] = part_open(model)
        json.dump(out, open(a.out, "w"), indent=1, default=float)
    if a.open_null:
        out["DIAG_open_window_null_symmetry"] = part_open_null_symmetry(model)
        json.dump(out, open(a.out, "w"), indent=1, default=float)
    if a.g3:
        out["g3_symmetry"] = part_g3(model)
        json.dump(out, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
