#!/usr/bin/env python3
"""M7 encoder training set (DESIGN §4 L3): positives and INCORRECT-TARGET negatives only.  Decoy-recipe sequences are
never used: every candidate's L1 row is checked against labels.arrow (is_decoy False) before it is admitted.

Per positive spectrum (scripts/m7/prep_train.py --label: PSM winner, target, peptidoform accepted at q <= 0.1%,
score >= that list's threshold):
  kind 0  positive              the confident target peptidoform, charge of its hypothesis
  kind 1  hard (top-K)          other TARGET rows of the same spectrum's top-20 (search report_psms 20) whose
                                I/L-collapsed stripped sequence differs from the positive; up to 6 (seeded sample)
  kind 2  shifted-precursor     2 confident target peptidoforms of OTHER spectra of the same species pool whose mass is
                                m_pos +- U(5, 60) Da (nearest in the pool), I/L-different.  Every confident peptide
                                therefore also appears as a negative ("prior balancing": sequence identity alone
                                carries ~no label information)
  kind 3  wrong-species         human/mouse spectra only: top-3 TARGET rows from a search of the run against an
                                E. coli-only index (no contaminants), minus any I/L-sequence present in the human or
                                mouse target space; up to 2
Dedup: exact spectrum hash (vulcan.encoders.data.spectrum_hash) -- duplicates inside training dropped (first kept in
manifest order) and any training spectrum whose hash occurs in an evaluation store (D0, D1, D7, DN, D8 holdout) dropped;
raw-file sha256 compared against every dev/sealed/null/D8 manifest sha256.
Folds: train_m7_folds.tsv (project-disjoint, 5 folds).
Output: $W/dataset/{train.npz, meta.json}; results dir gets the hash list and the fold map copy.
"""
from __future__ import annotations

import csv
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
import m3_common as C  # noqa: E402
import m3_data as MD  # noqa: E402
import prep_train as PT  # noqa: E402

import numpy as np  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from vulcan.encoders.data import MAX_LEN, N_PEAKS, peptide_tokens, spectra_topk, spectrum_hash, stripped_il  # noqa: E402

W = PT.W
RES = os.path.join(ROOT, "results", "2026-10-28_M7_encoders")
SEED = 11
N_HARD, N_SHIFT, N_WS = 6, 2, 2
ECOLI_ONLY = os.path.join(W, "index", "ecoli_only_prod")
EVAL_STORES = (glob.glob("/tmp/vulcan_e2/store/*.vspec.arrow") + glob.glob("/tmp/vulcan_m3/store_null/*.vspec.arrow")
               + glob.glob("/tmp/vulcan_m3e/store_d8/*.vspec.arrow") + glob.glob("/tmp/vulcan_m3e/store_dn/*.vspec.arrow"))


def ensure_ecoli_only():
    if not os.path.exists(os.path.join(ECOLI_ONLY, "meta.json")):
        from vulcan.index.build import build
        build([os.path.join(ROOT, "data", "fasta", "ecoli_K12_UP000000625.fasta")], ECOLI_ONLY, {"entrap_seeds": []})
    return ECOLI_ONLY


def eval_hashes():
    from vulcan.gpu.engine import load_store
    p = os.path.join(W, "dataset", "eval_hashes.npy")
    if os.path.exists(p):
        return np.load(p)
    hs = []
    for s in sorted(EVAL_STORES):
        sp = load_store(s)
        M, I, _ = spectra_topk(sp)
        hs.append(spectrum_hash(M, I, sp.precursor_mz))
        print("eval hash", os.path.basename(s), len(sp), flush=True)
    h = np.unique(np.concatenate(hs))
    np.save(p, h)
    return h


def raw_sha_overlap():
    if not os.path.exists(os.path.join(ROOT, "data", "manifests", "train_m7.sha256")):
        return "pending: train_m7.sha256 not written yet"
    train = {}
    for line in open(os.path.join(ROOT, "data", "manifests", "train_m7.sha256")):
        if line.strip():
            h, f = line.split(None, 1)
            train[h] = f.strip()
    other = {}
    for f in glob.glob(os.path.join(ROOT, "data", "manifests", "*.sha256")):
        if os.path.basename(f) in ("train_m7.sha256", "train.sha256"):
            continue
        for line in open(f):
            parts = line.split(None, 1)
            if len(parts) == 2 and len(parts[0]) == 64:
                other[parts[0]] = (os.path.basename(f), parts[1].strip())
    dup_in_train = len(train) != sum(1 for line in open(os.path.join(ROOT, "data", "manifests", "train_m7.sha256")) if line.strip())
    return {"n_train_files": len(train), "n_other_manifest_hashes": len(other),
            "overlap": [[train[h], other[h]] for h in train if h in other], "duplicate_sha256_within_train": dup_in_train}


def main():
    t0 = time.time()
    os.makedirs(os.path.join(W, "dataset"), exist_ok=True)
    rows = PT.rows()
    fold_of = {r["accession"]: int(r["fold"]) for r in rows}
    run_info = {PT.run_name(r): r for r in rows}
    ok, frac = PT.class_ok()
    ev_h = eval_hashes()
    import torch
    from vulcan.gpu import engine as E
    from vulcan.gpu.engine import load_store
    from vulcan.index.build import L1Index
    rng = np.random.default_rng(SEED)
    l1_ec = L1Index(ensure_ecoli_only())
    lab_ec = MD.labels(ECOLI_ONLY)
    ec_tgt = np.flatnonzero(~lab_ec["is_decoy"])
    p_ws = E.SearchParams(report_psms=3)
    gi_ec = E.GpuIndex(l1_ec, p_ws, C.STATIC, C.VAR)
    host_il = set()
    for sp_h in ("human", "mouse"):           # I/L target space of both host species, before any wrong-species pick
        _l1 = L1Index(C.index_dir(PT.FAM[sp_h], "prod"))
        _lab = MD.labels(C.index_dir(PT.FAM[sp_h], "prod"))
        host_il |= set(stripped_il(_l1, np.flatnonzero(~_lab["is_decoy"])))
        del _l1, _lab

    S_mz, S_int, S_n, S_ch, S_hash, S_run, S_fold, S_spec = [], [], [], [], [], [], [], []
    C_tok, C_ch, C_sid, C_lab, C_kind = [], [], [], [], []
    n_spec = 0
    stats = {"runs": {}, "dropped_dup_within_train": 0, "dropped_in_eval": 0}
    seen = set()
    order = ["human", "mouse", "ecoli"]
    for sp_name in order:
        fam = PT.FAM[sp_name]
        l1 = L1Index(C.index_dir(fam, "prod"))
        lab = MD.labels(C.index_dir(fam, "prod"))
        pool_rows, pool_mass = [], []
        pending = []      # (spectrum global ids, positive rows) for shifted negatives after the pool is complete
        for r in [x for x in rows if x["species"] == sp_name]:
            run = PT.run_name(r)
            lf = os.path.join(PT.psm_dir(sp_name), run + ".label.parquet")
            if run not in ok or not os.path.exists(lf):
                stats["runs"][run] = {"skipped": "not OT-HCD" if run not in ok else "no labels"}
                continue
            L = pq.read_table(lf).to_pydict()
            t = pq.read_table(os.path.join(PT.psm_dir(sp_name), run + ".psms.parquet"),
                              columns=["spec_idx", "pf_row", "charge", "rank"])
            ps = {k: t.column(k).to_numpy() for k in ("spec_idx", "pf_row", "charge")}
            spc = load_store(PT.store_of(run_info[run]))
            M, I, n = spectra_topk(spc)
            pos_spec = np.asarray(L["spec_idx"], np.int64)
            pos_row = np.asarray(L["row"], np.int64)
            h = spectrum_hash(M[pos_spec], I[pos_spec], spc.precursor_mz[pos_spec])
            keep = np.ones(len(pos_spec), bool)
            for i, hv in enumerate(h):
                if int(hv) in seen:
                    keep[i] = False; stats["dropped_dup_within_train"] += 1
                seen.add(int(hv))
            ine = np.isin(h, ev_h)
            stats["dropped_in_eval"] += int((ine & keep).sum())
            keep &= ~ine
            pos_spec, pos_row, h = pos_spec[keep], pos_row[keep], h[keep]
            # group PSM rows by spectrum
            o = np.argsort(ps["spec_idx"], kind="stable")
            ss = ps["spec_idx"][o]
            starts = np.searchsorted(ss, pos_spec, "left")
            ends = np.searchsorted(ss, pos_spec, "right")
            pos_pf = ps["pf_row"][pos_row]
            pos_il = stripped_il(l1, pos_pf)
            # wrong-species search (human/mouse)
            ws = {}
            if sp_name in ("human", "mouse"):
                tbl, _, _ = E.search(spc, gi_ec, l1_ec, p_ws)
                wsi = tbl.column("spec_idx").to_numpy()
                wsr = tbl.column("pf_row").to_numpy()
                wsc = tbl.column("charge").to_numpy()
                tg = ~lab_ec["is_decoy"][wsr]
                for a, b, c in zip(wsi[tg], wsr[tg], wsc[tg]):
                    ws.setdefault(int(a), []).append((int(b), int(c)))
            n_k = np.zeros(4, np.int64)
            for i in range(len(pos_spec)):
                gid = n_spec
                n_spec += 1
                si = pos_spec[i]
                S_spec.append(si)
                S_run.append(run)
                S_fold.append(fold_of[r["accession"]])
                cands = [(pos_pf[i], int(ps["charge"][pos_row[i]]), 1, 0)]
                rr = o[starts[i]:ends[i]]
                hard = [j for j in rr if j != pos_row[i] and not lab["is_decoy"][ps["pf_row"][j]]]
                if hard:
                    hil = stripped_il(l1, ps["pf_row"][hard])
                    hard = [j for j, s_ in zip(hard, hil) if s_ != pos_il[i]]
                if len(hard) > N_HARD:
                    hard = list(rng.choice(hard, N_HARD, replace=False))
                cands += [(int(ps["pf_row"][j]), int(ps["charge"][j]), 0, 1) for j in sorted(hard)]
                if ws.get(int(si)):
                    wl = ws[int(si)]
                    wil = stripped_il(l1_ec, [b for b, _ in wl])
                    wl = [x for x, s_ in zip(wl, wil) if s_ not in host_il][:N_WS]
                    cands += [(b, c, 0, 3) for b, c in wl]
                pending.append((gid, pos_pf[i], int(ps["charge"][pos_row[i]]), pos_il[i]))
                for pf, ch, y, kd in cands:
                    C_sid.append(gid); C_ch.append(ch); C_lab.append(y); C_kind.append(kd)
                    n_k[kd] += 1
                    C_tok.append((sp_name if kd != 3 else "ecoli_only", pf))
                pool_rows.append(pos_pf[i]); pool_mass.append(float(l1.mass32[pos_pf[i]]))
            S_mz.append(M[pos_spec]); S_int.append(I[pos_spec]); S_n.append(n[pos_spec])
            S_ch.append(np.array([int(ps["charge"][j]) for j in pos_row], np.int8)); S_hash.append(h)
            stats["runs"][run] = {"species": sp_name, "fold": fold_of[r["accession"]], "n_pos": int(len(pos_spec)),
                                  "n_by_kind": n_k.tolist()}
            print(run, stats["runs"][run], round(time.time() - t0), flush=True)
        # shifted-precursor (prior-balanced) negatives from the species pool
        pr = np.asarray(pool_rows, np.int64)
        pm = np.asarray(pool_mass)
        o = np.argsort(pm)
        pr, pm = pr[o], pm[o]
        pil = stripped_il(l1, pr)
        n_sh = 0
        for gid, pf, ch, il in pending:
            m0 = float(l1.mass32[pf])
            got = 0
            for _ in range(6):
                if got >= N_SHIFT:
                    break
                tgt_m = m0 + rng.choice([-1, 1]) * rng.uniform(5, 60)
                j = int(np.clip(np.searchsorted(pm, tgt_m), 0, len(pm) - 1))
                if abs(pm[j] - m0) < 4.5 or pil[j] == il:
                    continue
                C_sid.append(gid); C_ch.append(ch); C_lab.append(0); C_kind.append(2)
                C_tok.append((sp_name, int(pr[j])))
                got += 1; n_sh += 1
        stats.setdefault("n_shifted", {})[sp_name] = n_sh
        # tokenise this species' candidates now (index-specific rows)
        del l1
    # tokens
    l1s = {k: None for k in ("human", "mouse", "ecoli", "ecoli_only")}
    tok = np.zeros((len(C_tok), MAX_LEN), np.uint8)
    keys = np.array([k for k, _ in C_tok])
    rws = np.array([r_ for _, r_ in C_tok], np.int64)
    for k in l1s:
        m = keys == k
        if not m.any():
            continue
        l1k = L1Index(ECOLI_ONLY if k == "ecoli_only" else C.index_dir(PT.FAM[k], "prod"))
        labk = MD.labels(ECOLI_ONLY if k == "ecoli_only" else C.index_dir(PT.FAM[k], "prod"))
        assert not labk["is_decoy"][rws[m]].any(), "decoy-recipe sequence in training candidates"
        tok[m] = peptide_tokens(l1k, rws[m])
    out = dict(spec_mz=np.concatenate(S_mz), spec_int=np.concatenate(S_int), spec_n=np.concatenate(S_n),
               spec_charge=np.concatenate(S_ch), spec_hash=np.concatenate(S_hash), spec_fold=np.array(S_fold, np.int8),
               spec_run=np.array(S_run), cand_tok=tok, cand_charge=np.array(C_ch, np.int8),
               cand_sid=np.array(C_sid, np.int64), cand_label=np.array(C_lab, np.int8), cand_kind=np.array(C_kind, np.int8))
    np.savez(os.path.join(W, "dataset", "train.npz"), **out)
    stats["n_spectra"] = int(len(out["spec_fold"]))
    stats["n_candidates"] = int(len(out["cand_sid"]))
    stats["n_by_kind"] = np.bincount(out["cand_kind"], minlength=4).tolist()
    stats["n_spectra_by_fold"] = np.bincount(out["spec_fold"], minlength=5).tolist()
    stats["raw_sha256"] = raw_sha_overlap()
    stats["n_eval_spectrum_hashes"] = int(len(ev_h))
    stats["eval_stores"] = sorted(os.path.basename(s) for s in EVAL_STORES)
    stats["ot_hcd_fraction"] = frac
    stats["wall_s"] = time.time() - t0
    json.dump(stats, open(os.path.join(W, "dataset", "meta.json"), "w"), indent=1, default=str)
    os.makedirs(RES, exist_ok=True)
    import pyarrow as pa
    pq.write_table(pa.table({"spec_hash": out["spec_hash"], "run": out["spec_run"], "fold": out["spec_fold"]}),
                   os.path.join(RES, "train_spectrum_hashes.parquet"), compression="zstd")
    print(json.dumps({k: v for k, v in stats.items() if k not in ("runs", "ot_hcd_fraction")}, default=str)[:2000])
    del gi_ec
    torch.cuda.empty_cache()


if __name__ == "__main__":
    main()
