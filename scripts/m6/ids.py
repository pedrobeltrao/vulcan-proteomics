#!/usr/bin/env python3
"""Experiment-wide identifications for M6 quant (L4): paired double competition over all runs of a PSM dir
(peptidoform units, pair = tdc_pair, hash tie-break seed 11), score = frozen static_gbm_pred (or hyperscore).

Output <out>/ids.parquet: one row per PSM-competition winner (best candidate of each spectrum by the score) with
  run, spec_idx, scan, rt, pf_row, pf_id, charge, isotope_error, calcmass, expmass, score, is_decoy,
  pf_q (q of the peptidoform's pair winner; 1 if the peptidoform lost its pair competition),
  pass_{alpha} = pf_q <= alpha AND score >= the alpha score threshold (targets: run-level evidence of an accepted
  peptidoform; decoys: the decoy peptides that reach the same threshold, used as PIP-ECHO decoy-peptide donors).
and <out>/peptides.parquet: pf_id -> modified sequence, target sequence, proteins, species (HUMAN/YEAST/ECOLI/..),
  contaminant flag, is_decoy.
    python scripts/m6/ids.py --psms DIR --index INDEX --out OUT [--score static|hyperscore]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

import numpy as np  # noqa: E402
import pyarrow as pa  # noqa: E402
import pyarrow.ipc as ipc  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from vulcan.index.build import L1Index  # noqa: E402
from vulcan.stats._util import hash_tiebreak  # noqa: E402
from vulcan.stats.tdc import double_competition  # noqa: E402

ALPHAS = (0.001, 0.01, 0.05)
COLS = ["run", "spec_idx", "scan", "rt", "rank", "pf_row", "pf_id", "charge", "isotope_error", "calcmass", "expmass",
        "hyperscore", "delta_mass_ppm"]


def labels(index_dir):
    with ipc.open_file(os.path.join(index_dir, "labels.arrow")) as r:
        t = r.read_all()
    return {k: t.column(k).to_numpy(zero_copy_only=False) for k in ("is_decoy", "tdc_pair", "target_seq")}


def protein_map(index_dir):
    with ipc.open_file(os.path.join(index_dir, "proteins.arrow")) as rd:
        t = rd.read_all()
    pep = t.column("peptide").to_pylist()
    off = t.column("prot_off").to_numpy()
    ln = t.column("prot_len").to_numpy()
    pidx = np.load(os.path.join(index_dir, "protein_idx.npy"))
    acc = open(os.path.join(index_dir, "accessions.txt")).read().split("\n")
    return dict(zip(pep, range(len(pep)))), off, ln, pidx, acc


def species_of(acc: str) -> str:
    tail = acc.rsplit("_", 1)[-1] if "_" in acc else "?"
    return tail


def modseq(l1, row, var):
    s = l1.sequence(row)
    m = int(np.uint64(l1.mask[row]))
    out = []
    for k, c in enumerate(s):
        out.append(c)
        if (m >> k) & 1:
            out.append(f"[{var[c]:+.4f}]")
    return "".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--psms", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--score", default="static")
    ap.add_argument("--runs", nargs="*", default=None, help="subset of runs (default all in --psms)")
    a = ap.parse_args()
    t0 = time.perf_counter()
    files = sorted(glob.glob(os.path.join(a.psms, "*.psms.parquet")))
    parts = []
    for f in files:
        run = os.path.basename(f)[:-len(".psms.parquet")]
        if a.runs and run not in a.runs:
            continue
        t = pq.read_table(f, columns=COLS)
        d = {c: t.column(c).to_numpy(zero_copy_only=False) for c in COLS}
        if a.score == "static":
            d["score"] = pq.read_table(os.path.join(a.psms, run + ".score.parquet")).column("score").to_numpy()
        else:
            d["score"] = d["hyperscore"].astype(np.float64)
        parts.append(d)
    P = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    lab = labels(a.index)
    r = P["pf_row"]
    is_dec = lab["is_decoy"][r].astype(bool)
    pair = lab["tdc_pair"][r]
    spec = np.char.add(np.char.add(P["run"].astype(str), ":"), P["spec_idx"].astype(str))
    tb = hash_tiebreak(P["pf_id"], C.SEED)
    res = double_competition(spec, P["pf_id"], pair, P["score"], is_dec, seed=C.SEED, tiebreak=tb)
    pfq = {}
    for row, q in zip(res.row, res.q):
        pfq[int(P["pf_id"][row])] = float(q)
    w1 = res.psm_winners
    out = {k: v[w1] for k, v in P.items() if k != "rank"}
    out["is_decoy"] = is_dec[w1]
    out["pf_q"] = np.array([pfq.get(int(u), 1.0) for u in out["pf_id"]])
    dec_winner = np.isin(out["pf_id"], P["pf_id"][res.row[res.is_decoy]])
    summary = {"n_psm_winners": int(len(w1)), "runs": sorted(set(P["run"].tolist())), "score": a.score}
    for al in ALPHAS:
        acc = (res.q <= al) & ~res.is_decoy
        thr = float(res.score[acc].min()) if acc.any() else np.inf
        tpass = ~out["is_decoy"] & (out["pf_q"] <= al) & (out["score"] >= thr)
        dpass = out["is_decoy"] & dec_winner & (out["score"] >= thr)
        out[f"pass_{al}"] = tpass | dpass
        summary[str(al)] = {"peptidoforms": int(acc.sum()), "score_threshold": thr,
                            "decoy_pair_winners_above": int(((res.score >= thr) & res.is_decoy).sum())}
    tbl = pa.table(out)
    os.makedirs(a.out, exist_ok=True)
    pq.write_table(tbl, os.path.join(a.out, "ids.parquet"))
    # peptide annotation for every peptidoform that appears among the winners
    l1 = L1Index(a.index)
    var = {"M": C.VAR["M"]}
    pos, off, ln, pidx, acc = protein_map(a.index)
    rows_u, first = np.unique(out["pf_row"], return_index=True)
    pep = {"pf_id": [], "pf_row": [], "modseq": [], "target_seq": [], "proteins": [], "species": [],
           "contaminant": [], "is_decoy": []}
    for rr in rows_u:
        rr = int(rr)
        ts = lab["target_seq"][rr]
        i = pos.get(ts)
        prots = [] if i is None else [acc[j] for j in pidx[off[i]:off[i] + ln[i]]]
        cont = any("Cont_" in p for p in prots)
        sp = sorted({species_of(p) for p in prots if "Cont_" not in p})
        pep["pf_id"].append(int(l1.pf_id[rr]))
        pep["pf_row"].append(rr)
        pep["modseq"].append(modseq(l1, rr, var))
        pep["target_seq"].append(ts)
        pep["proteins"].append(";".join(prots))
        pep["species"].append(sp[0] if len(sp) == 1 else ("MULTI" if sp else "NONE"))
        pep["contaminant"].append(cont)
        pep["is_decoy"].append(bool(lab["is_decoy"][rr]))
    pep["pf_id"] = np.array(pep["pf_id"], np.uint64)
    pq.write_table(pa.table(pep), os.path.join(a.out, "peptides.parquet"))
    summary["wall_s"] = time.perf_counter() - t0
    json.dump(summary, open(os.path.join(a.out, "ids_summary.json"), "w"), indent=1)
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
