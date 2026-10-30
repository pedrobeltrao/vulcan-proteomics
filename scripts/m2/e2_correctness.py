#!/usr/bin/env python3
"""E2 correctness: (A) G0 label-swap / storage-permutation invariance, (B) K1 index counts vs brute force,
(C) top-K recall vs exhaustive K2 scoring for targets and decoys, (D) agreement with Sage 0.14.7 (oracle).

    python scripts/m2/e2_correctness.py [--parts A B C D] [--out results/.../correctness.json]

Labels (is_decoy etc.) are read here only to *evaluate* (L4 role); the engine never sees them.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e2_common as C  # noqa: E402

import torch  # noqa: E402

from vulcan.gpu import engine as E  # noqa: E402
from vulcan.gpu import sage_bridge as SB  # noqa: E402
from vulcan.index.build import L1Index, _read, build  # noqa: E402

P = E.SearchParams()


def labels_by_row(index_dir):
    lab = _read(os.path.join(index_dir, "labels.arrow"))
    return {k: lab.column(k).to_numpy(zero_copy_only=False) for k in ("is_decoy", "is_entrapment", "family")}


# ------------------------------------------------------------------------------------------------- A: G0
def part_a():
    """Shuffle FASTA protein order + permute the pair list + swap target/decoy roles -> index side and every
    L2/L3 output must be bit-identical."""
    ds = "D1"
    base = C.index_dir(ds, f"E{C.ENTRAP_SEED}")
    work = os.path.join(C.WORK, "g0")
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    # shuffled FASTA (entry order)
    txt = open(C.fasta_path(ds)).read().split("\n>")
    ents = [txt[0]] + [">" + e for e in txt[1:]]
    ents = [e if e.endswith("\n") else e + "\n" for e in ents]
    rng = np.random.default_rng(11)
    shuf = os.path.join(work, "shuffled.fasta")
    with open(shuf, "w") as fh:
        fh.write("".join(ents[i] for i in rng.permutation(len(ents))))
    alt = os.path.join(work, "index_swapped")
    build([shuf], alt, {"entrap_seeds": [C.ENTRAP_SEED]}, shuffle_rows_seed=23, swap_labels=True)
    a, b = L1Index(base), L1Index(alt)
    res = {"index_side_identical": {}}
    for k in ("pf_id", "mass32", "seq_off", "length", "mask", "n_var", "missed", "partner", "residues"):
        res["index_side_identical"][k] = bool(np.array_equal(getattr(a, k), getattr(b, k)))
    la, lb = labels_by_row(base), labels_by_row(alt)
    res["labels_identical_after_swap(true labels preserved)"] = all(np.array_equal(la[k], lb[k]) for k in la)
    sp = E.load_store(C.store_path(C.DATASETS[ds]["raw"][0]))
    outs = []
    for l1 in (a, b):
        gi = E.GpuIndex(l1, P, C.STATIC, C.VAR)
        t, _, _ = E.search(sp, gi, l1, P, report=5)
        outs.append(t)
        del gi
        torch.cuda.empty_cache()
    res["psm_table_rows"] = outs[0].num_rows
    res["psm_tables_bit_identical"] = bool(outs[0].equals(outs[1]))
    res["pass"] = res["psm_tables_bit_identical"] and all(res["index_side_identical"].values())
    return res


# ------------------------------------------------------------------------------------------------- B: K1 vs brute force
def part_b(n_sample=3000):
    res = {}
    for ds, kind in (("D1", f"E{C.ENTRAP_SEED}"), ("D0", f"E{C.ENTRAP_SEED}")):
        l1 = L1Index(C.index_dir(ds, kind))
        gi = E.GpuIndex(l1, P, C.STATIC, C.VAR)
        sp = E.load_store(C.store_path(C.DATASETS[ds]["raw"][0]))
        st = E.Stats()
        peaks = E.k0_preprocess(sp, P, gi.dev, st)
        hyp = E.make_hypotheses(sp, peaks["valid"], gi.mass32_np, P)
        rng = np.random.default_rng(47)
        specs = np.unique(hyp["spec"])
        sub = np.sort(rng.choice(specs, min(n_sample, len(specs)), replace=False))
        o_h, o_r, o_c, scored, matched, d_spec, d_maxfz = E.run_k1(gi, hyp, peaks, P, st)
        h, r = E.window_candidates(hyp, gi.dev, sub)
        bc = E.k1_brute(gi, hyp, peaks, h, r, P).long()
        keep = torch.from_numpy(np.isin(hyp["spec"], sub)).to(gi.dev)[o_h.long()]
        ki = (o_h.long()[keep] << 32) | o_r.long()[keep]
        kv = o_c.long()[keep]
        bi = (h << 32) | r
        nz = bc > 0
        o1 = torch.argsort(ki)
        o2 = torch.argsort(bi[nz])
        same_set = torch.equal(ki[o1], bi[nz][o2])
        same_cnt = same_set and torch.equal(kv[o1], bc[nz][o2])
        # per-hypothesis Sage stats (scored candidates, matched peaks)
        hs = torch.zeros(len(hyp["spec"]), dtype=torch.int64, device=gi.dev).index_add_(0, h[nz], torch.ones_like(h[nz]))
        hm = torch.zeros(len(hyp["spec"]), dtype=torch.int64, device=gi.dev).index_add_(0, h, bc)
        hsel = torch.from_numpy(np.isin(hyp["spec"], sub)).to(gi.dev)
        res[ds] = {"index": kind, "n_spectra_sampled": int(len(sub)), "n_window_pairs": int(len(h)),
                   "n_nonzero_bruteforce": int(nz.sum()), "n_appended_k1": int(keep.sum()),
                   "same_candidate_set": bool(same_set), "same_counts": bool(same_cnt),
                   "same_scored_candidates": bool(torch.equal(hs[hsel], scored.long()[hsel])),
                   "same_matched_peaks": bool(torch.equal(hm[hsel], matched.long()[hsel]))}
        if not same_cnt and same_set is False:
            a = set(ki.tolist()); b = set(bi[nz].tolist())
            res[ds]["only_k1"] = len(a - b)
            res[ds]["only_brute"] = len(b - a)
        del gi
        torch.cuda.empty_cache()
    res["pass"] = all(v["same_counts"] for k, v in res.items() if isinstance(v, dict))
    return res


# ------------------------------------------------------------------------------------------------- C: recall
def _best_per_label(spec, score, is_dec, n):
    bt = np.full(n, -np.inf)
    bd = np.full(n, -np.inf)
    np.maximum.at(bt, spec[~is_dec], score[~is_dec])
    np.maximum.at(bd, spec[is_dec], score[is_dec])
    return bt, bd


def _rescored(dbg, l1, P):
    k2 = dbg["k2"]
    ok = ((k2["matched_b"].long() + k2["matched_y"].long()) >= P.min_matched_peaks).cpu().numpy()
    spec = k2["spec"].cpu().numpy()[ok]
    hs = k2["hyperscore"].cpu().numpy()[ok]
    row = dbg["r"].cpu().numpy()[ok]
    return spec, hs, row


def part_c():
    res = {}
    for ds, kind, files in (("D1", "prod", 1), ("D1", f"E{C.ENTRAP_SEED}", 1), ("D0", "prod", 1),
                            ("D0", f"E{C.ENTRAP_SEED}", 1)):
        d = C.index_dir(ds, kind)
        l1 = L1Index(d)
        lab = labels_by_row(d)
        gi = E.GpuIndex(l1, P, C.STATIC, C.VAR)
        for raw in C.DATASETS[ds]["raw"][:files]:
            sp = E.load_store(C.store_path(raw))
            t0 = time.perf_counter()
            tk, _, dk = E.search(sp, gi, l1, P, report=1)
            tx, _, dx = E.search_exhaustive(sp, gi, l1, P, report=1)
            t_ex = time.perf_counter() - t0
            n = len(sp)
            sk, hk, rk = _rescored(dk, l1, P)
            sx, hx, rx = _rescored(dx, l1, P)
            btk, bdk = _best_per_label(sk, hk, lab["is_decoy"][rk], n)
            btx, bdx = _best_per_label(sx, hx, lab["is_decoy"][rx], n)
            has_t, has_d = np.isfinite(btx), np.isfinite(bdx)
            rec_t = float(np.mean(btk[has_t] == btx[has_t]))
            rec_d = float(np.mean(bdk[has_d] == bdx[has_d]))
            # top-1 (any label): same score (ties count as agreement), and same peptidoform
            a = tk.to_pandas().set_index("spec_idx")
            b = tx.to_pandas().set_index("spec_idx")
            j = b.join(a, rsuffix="_k", how="left")
            top_score = float(np.mean(j.hyperscore_k.fillna(-1).values == j.hyperscore.values))
            top_pf = float(np.mean(j.pf_row_k.fillna(-1).values == j.pf_row.values))   # rows < 2^27: exact in f64
            # among confident exhaustive winners (top-1 hyperscore >= 25% best quartile) for context
            q = j.hyperscore >= np.quantile(j.hyperscore, 0.5)
            res[f"{ds}/{kind}/{C.run_name(raw)}"] = {
                "n_spectra_with_psm": int(len(b)),
                "candidates_exhaustive": int(len(dx["r"])), "candidates_k1k2": int(len(dk["r"])),
                "recall_best_target": rec_t, "recall_best_decoy": rec_d, "abs_delta_t_d": abs(rec_t - rec_d),
                "top1_same_score": top_score, "top1_same_peptidoform": top_pf,
                "top1_same_score_upper_half": float(np.mean((j.hyperscore_k.fillna(-1) == j.hyperscore)[q])),
                "wall_s_both_searches": round(t_ex, 2)}
            print(ds, kind, res[f"{ds}/{kind}/{C.run_name(raw)}"], flush=True)
        del gi
        torch.cuda.empty_cache()
    res["pass"] = all(v["recall_best_target"] >= 0.995 and v["recall_best_decoy"] >= 0.995 and v["abs_delta_t_d"] <= 0.01
                      for k, v in res.items() if isinstance(v, dict))
    return res


# ------------------------------------------------------------------------------------------------- D: Sage oracle
def run_sage(ds, kind, outdir, threads, raws=None):
    raws = raws or C.DATASETS[ds]["raw"]
    os.makedirs(outdir, exist_ok=True)
    cfg = SB.sage_config(C.ensure_sage_fasta(ds, kind), [C.ensure_mgf(r) for r in raws], outdir, P, C.STATIC, C.VAR,
                         C.MAX_VAR)
    SB.write_config(cfg, os.path.join(outdir, "config.json"))
    env = {**os.environ, "RAYON_NUM_THREADS": str(threads)}
    with open(os.path.join(outdir, "sage.log"), "w") as fh:
        subprocess.run([os.path.join(C.ROOT, "tools", "bin", "sage"), os.path.join(outdir, "config.json"),
                        "--disable-telemetry-i-dont-want-to-improve-sage"], check=True, env=env, stdout=fh,
                       stderr=subprocess.STDOUT)
    return outdir


def compare_with_sage(ds, kind, sage_dir, raw, gi, l1, p=P):
    sp = E.load_store(C.store_path(raw))
    tbl, _, _ = E.search(sp, gi, l1, p)
    ours = tbl.to_pandas()
    ours["seq"] = [l1.sequence(r) for r in ours.pf_row]
    ours["mask"] = l1.mask[ours.pf_row]
    s = SB.read_results(sage_dir).to_pandas()
    s = s[s.filename == C.run_name(raw) + ".mgf"].copy()
    s["spec_idx"] = s.scannr.astype(int)
    pm = [SB.parse_peptide(x, C.VAR) for x in s.peptide]
    s["seq"] = [a for a, _ in pm]
    s["mask"] = np.array([b for _, b in pm], dtype=np.uint64)
    m = ours.merge(s[["spec_idx", "seq", "mask", "hyperscore", "matched_peaks", "delta_next", "longest_b",
                      "longest_y", "poisson", "scored_candidates", "matched_intensity_pct"]],
                   on="spec_idx", how="outer", suffixes=("", "_s"))
    both = m.dropna(subset=["hyperscore", "hyperscore_s"])
    same = (both.seq == both.seq_s).values & (both["mask"].astype(np.uint64).values == both.mask_s.astype(np.uint64).values)
    d = np.abs(both.hyperscore.values - both.hyperscore_s.values)
    feat = {}
    for f in ("matched_peaks", "longest_b", "longest_y", "scored_candidates"):
        feat[f + "_identical_frac(same pep)"] = float(np.mean(both[f].values[same] == both[f + "_s"].values[same]))
    for f in ("poisson", "matched_intensity_pct", "delta_next"):
        a, b = both[f].values[same].astype(float), both[f + "_s"].values[same].astype(float)
        feat[f + "_max_abs_diff(same pep)"] = float(np.nanmax(np.abs(a - b))) if len(a) else None
    return {
        "n_psms_ours": int(len(ours)), "n_psms_sage": int(len(s)), "n_both": int(len(both)),
        "only_ours": int(m.hyperscore_s.isna().sum()), "only_sage": int(m.hyperscore.isna().sum()),
        "top1_same_peptidoform": float(same.mean()),
        "same_pep_hyperscore_max_abs_diff": float(d[same].max()) if same.any() else None,
        "diff_pep_ours_higher": int(((both.hyperscore.values > both.hyperscore_s.values + 1e-9) & ~same).sum()),
        "diff_pep_sage_higher": int(((both.hyperscore_s.values > both.hyperscore.values + 1e-9) & ~same).sum()),
        "diff_pep_exact_tie": int(((d <= 1e-9) & ~same).sum()),
        **feat,
    }


def part_d():
    res = {}
    for ds in ("D1", "D0"):
        kind = "prod"
        sdir = os.path.join(C.WORK, "sage_runs", f"{ds}_{kind}_oracle")
        if not os.path.exists(os.path.join(sdir, "results.sage.tsv")):
            run_sage(ds, kind, sdir, 32)
        l1 = L1Index(C.index_dir(ds, kind))
        gi = E.GpuIndex(l1, P, C.STATIC, C.VAR)
        for raw in C.DATASETS[ds]["raw"]:
            for mode, p in (("sage_equivalent_K50", E.SearchParams(topk=50)), ("default_K500", P)):
                key = f"{ds}/{kind}/{C.run_name(raw)}/{mode}"
                res[key] = compare_with_sage(ds, kind, sdir, raw, gi, l1, p)
                print(key, res[key], flush=True)
        del gi
        torch.cuda.empty_cache()
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", nargs="+", default=["A", "B", "C", "D"])
    ap.add_argument("--out", default=os.path.join(C.RES, "correctness.json"))
    a = ap.parse_args()
    out = json.load(open(a.out)) if os.path.exists(a.out) else {}
    for part in a.parts:
        t = time.perf_counter()
        r = {"A": part_a, "B": part_b, "C": part_c, "D": part_d}[part]()
        r["wall_s"] = round(time.perf_counter() - t, 1)
        out[{"A": "A_G0_invariance", "B": "B_k1_vs_bruteforce", "C": "C_topk_recall_vs_exhaustive",
             "D": "D_sage_oracle"}[part]] = r
        print(part, json.dumps(r, indent=1, default=str), flush=True)
        os.makedirs(os.path.dirname(a.out), exist_ok=True)
        with open(a.out, "w") as fh:
            json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    main()
