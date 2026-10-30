#!/usr/bin/env python3
"""M5 identifications (L4): paired double competition (peptidoform units, pair = tdc_pair, hash tie-break seed 11),
score = frozen static_gbm_pred (or hyperscore), q-values by

  --strata lip   Group-walk over the specificity class of the pair (specific / N-semi / C-semi; vulcan.lip.digest),
                 the DESIGN L4 closed-search stratum "specificity"; pooled TDC q is kept as pf_q_pooled.
  --strata none  pooled TDC.
  --per-run      the whole competition is done separately for every run (CF-MS: run-specific per-fraction q-values,
                 DESIGN L5); the series-level (pooled over runs) q is stored as pf_q_series.

Output <out>/ids.parquet (one row per PSM-competition winner; schema of scripts/m6/ids.py plus `stratum`) and
<out>/peptides.parquet (pf_id -> modseq, target_seq, proteins, is_decoy, is_entrapment, cls).  Entrapment peptides
(family E) are annotated with NO protein (they are not in the proteome); their target partner is in `entrap_target`.
pass_{a}: target pair-winner with q <= a and score >= its stratum's acceptance threshold (run-level evidence); decoys
reaching the same threshold are flagged too (PIP-ECHO decoy-peptide donors for MBR).
    python scripts/m5lc/ids.py --psms DIR --index INDEX --out OUT [--runs ...] [--strata lip|none] [--per-run]
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
from vulcan.stats.groupwalk import group_walk, label_blind_pair_counts  # noqa: E402
from vulcan.stats.tdc import double_competition, tdc_qvalues  # noqa: E402

ALPHAS = (0.001, 0.01, 0.05)
COLS = ["run", "spec_idx", "scan", "rt", "rank", "pf_row", "pf_id", "charge", "isotope_error", "calcmass", "expmass",
        "hyperscore", "delta_mass_ppm"]
CLS_NAME = {0: "specific", 1: "n_semi", 2: "c_semi"}


def labels(index_dir):
    with ipc.open_file(os.path.join(index_dir, "labels.arrow")) as r:
        t = r.read_all()
    return {k: t.column(k).to_numpy(zero_copy_only=False)
            for k in ("is_decoy", "tdc_pair", "target_seq", "is_entrapment", "fwd_seq")}


def protein_map(index_dir):
    with ipc.open_file(os.path.join(index_dir, "proteins.arrow")) as rd:
        t = rd.read_all()
    pep = t.column("peptide").to_pylist()
    off = t.column("prot_off").to_numpy()
    ln = t.column("prot_len").to_numpy()
    pidx = np.load(os.path.join(index_dir, "protein_idx.npy"))
    acc = open(os.path.join(index_dir, "accessions.txt")).read().split("\n")
    return dict(zip(pep, range(len(pep)))), off, ln, pidx, acc


def modseq(l1, row, var):
    s = l1.sequence(row)
    m = int(np.uint64(l1.mask[row]))
    out = []
    for k, c in enumerate(s):
        out.append(c)
        if (m >> k) & 1:
            out.append(f"[{var[c]:+.4f}]")
    return "".join(out)


def compete(P, lab, strata_of_row, seed):
    """Double competition + (Group-walk | pooled) q.  Returns dict with winners and per-peptidoform q."""
    r = P["pf_row"]
    is_dec = lab["is_decoy"][r].astype(bool)
    pair = lab["tdc_pair"][r]
    spec = np.char.add(np.char.add(P["run"].astype(str), ":"), P["spec_idx"].astype(str))
    tb = hash_tiebreak(P["pf_id"], seed)
    res = double_competition(spec, P["pf_id"], pair, P["score"], is_dec, seed=seed, tiebreak=tb)
    out = {"res": res, "is_dec": is_dec, "q_pooled": res.q}
    if strata_of_row is not None:
        g = strata_of_row[res.row]
        out["counts"] = {str(k): int(v) for k, v in label_blind_pair_counts(res.pair, g).items()}
        out["q"] = group_walk(res.score, res.is_decoy, g, seed=seed)
        out["stratum"] = g
    else:
        out["q"] = res.q
        out["stratum"] = np.zeros(len(res.row), np.int64)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--psms", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--score", default="static")
    ap.add_argument("--runs", nargs="*", default=None)
    ap.add_argument("--strata", default="none", choices=["none", "lip"])
    ap.add_argument("--per-run", action="store_true")
    a = ap.parse_args()
    t0 = time.perf_counter()
    parts = []
    for f in sorted(glob.glob(os.path.join(a.psms, "*.psms.parquet"))):
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
    cls_row = None
    if a.strata == "lip":
        from vulcan.lip.index import load_class
        cmap = load_class(a.index)
        ts = lab["target_seq"][P["pf_row"]]
        cls_row = np.array([cmap.get(s, 0) for s in ts], np.int64)
    runs = sorted(set(P["run"].tolist()))
    groups = [np.flatnonzero(P["run"] == r) for r in runs] if a.per_run else [np.arange(len(P["run"]))]
    out_parts, summary = [], {"runs": runs, "score": a.score, "strata": a.strata, "per_run": a.per_run, "groups": {}}
    series = compete(P, lab, cls_row, C.SEED) if a.per_run else None
    series_q = {}
    if series is not None:
        for row, q in zip(series["res"].row, series["q"]):
            series_q[int(P["pf_id"][row])] = float(q)
    for gi, idx in enumerate(groups):
        Pg = {k: v[idx] for k, v in P.items()}
        cg = cls_row[idx] if cls_row is not None else None
        c = compete(Pg, lab, cg, C.SEED)
        res = c["res"]
        pfq, pfqp, pfs = {}, {}, {}
        for row, q, qp, s in zip(res.row, c["q"], c["q_pooled"], c["stratum"]):
            u = int(Pg["pf_id"][row])
            pfq[u], pfqp[u], pfs[u] = float(q), float(qp), int(s)
        w1 = res.psm_winners
        o = {k: v[w1] for k, v in Pg.items() if k != "rank"}
        o["is_decoy"] = c["is_dec"][w1]
        o["pf_q"] = np.array([pfq.get(int(u), 1.0) for u in o["pf_id"]])
        o["pf_q_pooled"] = np.array([pfqp.get(int(u), 1.0) for u in o["pf_id"]])
        if series is not None:
            o["pf_q_series"] = np.array([series_q.get(int(u), 1.0) for u in o["pf_id"]])
        st_w = cg[w1] if cg is not None else np.zeros(len(w1), np.int64)
        o["stratum"] = st_w
        dec_winner = np.isin(o["pf_id"], Pg["pf_id"][res.row[res.is_decoy]])
        gsum = {"n_psm_winners": int(len(w1)), "strata_counts": c.get("counts")}
        for al in ALPHAS:
            acc = (c["q"] <= al) & ~res.is_decoy
            tpass = np.zeros(len(w1), bool)
            dpass = np.zeros(len(w1), bool)
            rec = {"peptidoforms": int(acc.sum()), "pooled_tdc_peptidoforms": int(((c["q_pooled"] <= al) & ~res.is_decoy).sum()),
                   "by_stratum": {}}
            for s in np.unique(c["stratum"]):
                ins = c["stratum"] == s
                thr = float(res.score[acc & ins].min()) if (acc & ins).any() else np.inf
                ws = st_w == s
                tpass |= ws & ~o["is_decoy"] & (o["pf_q"] <= al) & (o["score"] >= thr)
                dpass |= ws & o["is_decoy"] & dec_winner & (o["score"] >= thr)
                rec["by_stratum"][CLS_NAME.get(int(s), str(s)) if a.strata == "lip" else "all"] = {
                    "accepted": int((acc & ins).sum()), "score_threshold": thr,
                    "decoy_pair_winners_above": int(((res.score >= thr) & res.is_decoy & ins).sum()),
                    "pooled_tdc_accepted": int(((c["q_pooled"] <= al) & ~res.is_decoy & ins).sum()),
                    "separate_tdc_accepted": int(((tdc_qvalues(res.score[ins], res.is_decoy[ins]) <= al)
                                                  & ~res.is_decoy[ins]).sum())}
            o[f"pass_{al}"] = tpass | dpass
            gsum[str(al)] = rec
        summary["groups"][runs[gi] if a.per_run else "all"] = gsum
        out_parts.append(o)
    out = {k: np.concatenate([p[k] for p in out_parts]) for k in out_parts[0]}
    os.makedirs(a.out, exist_ok=True)
    pq.write_table(pa.table(out), os.path.join(a.out, "ids.parquet"))
    # peptide annotation
    l1 = L1Index(a.index)
    pos, off, ln, pidx, acc = protein_map(a.index)
    cmap = None
    if a.strata == "lip":
        from vulcan.lip.index import load_class
        cmap = load_class(a.index)
    rows_u = np.unique(out["pf_row"])
    pep = {k: [] for k in ("pf_id", "pf_row", "modseq", "seq", "target_seq", "proteins", "is_decoy", "is_entrapment",
                           "entrap_target", "cls")}
    for rr in rows_u:
        rr = int(rr)
        ent = bool(lab["is_entrapment"][rr])
        ts = lab["target_seq"][rr]
        i = None if ent else pos.get(ts)
        prots = [] if i is None else [acc[j] for j in pidx[off[i]:off[i] + ln[i]]]
        pep["pf_id"].append(int(l1.pf_id[rr]))
        pep["pf_row"].append(rr)
        pep["modseq"].append(modseq(l1, rr, C.VAR))
        pep["seq"].append(l1.sequence(rr))
        pep["target_seq"].append(ts)
        pep["proteins"].append(";".join(prots))
        pep["is_decoy"].append(bool(lab["is_decoy"][rr]))
        pep["is_entrapment"].append(ent)
        pep["entrap_target"].append(ts if ent else "")
        pep["cls"].append(int(cmap.get(ts, 0)) if cmap is not None else 0)
    pep["pf_id"] = np.array(pep["pf_id"], np.uint64)
    pq.write_table(pa.table(pep), os.path.join(a.out, "peptides.parquet"))
    summary["wall_s"] = time.perf_counter() - t0
    json.dump(summary, open(os.path.join(a.out, "ids_summary.json"), "w"), indent=1)
    print(json.dumps({k: summary[k] for k in ("runs", "wall_s")}))
    if not a.per_run:
        print(json.dumps(summary["groups"]["all"], indent=1))


if __name__ == "__main__":
    main()
