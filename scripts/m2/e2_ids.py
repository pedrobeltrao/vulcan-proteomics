#!/usr/bin/env python3
"""E2 first end-to-end ID count (dev only): peptide-level paired TDC at q <= 1% (vulcan.stats.double_competition).

Arms per dataset (all files of the dataset pooled, top-1 PSM per spectrum):
  ours_hyperscore      our K1+K2 PSMs, score = hyperscore, paired double competition (pair = target<->own decoy)
  sage_hyperscore      Sage's PSMs (same spectra, same index) through the SAME estimator and score
  sage_own_q           Sage's own peptide_q <= 1% (its LDA discriminant + its peptide-level TDC)  [reference only]
plus, from a search of the entrapment-augmented index (T + E_2001 + decoys), the paired entrapment FDP estimate
of our hyperscore list at 1% (dev diagnostic; DESIGN §6.3 arm P, one seed).

    python scripts/m2/e2_ids.py --ours DIR_WITH_psms.parquet --sage SAGE_OUTDIR --dataset D0 [--entrap-ours DIR]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e2_common as C  # noqa: E402

from vulcan.gpu import sage_bridge as SB  # noqa: E402
from vulcan.index.build import L1Index, _read  # noqa: E402
from vulcan.stats import entrapment as ENT  # noqa: E402
from vulcan.stats._util import hash_tiebreak  # noqa: E402
from vulcan.stats.tdc import double_competition  # noqa: E402

ALPHA = 0.01
TB_SEED = 11


def labels(index_dir):
    lab = _read(os.path.join(index_dir, "labels.arrow"))
    return {k: lab.column(k).to_numpy(zero_copy_only=False) for k in
            ("pf_id", "is_decoy", "is_entrapment", "tdc_pair", "entrap_pair", "family")}


def paired_tdc(spec_key, rows, score, lab, l1):
    """rows: index rows of each PSM. Returns (result, accepted target rows mask over winners)."""
    unit = l1.pf_id[rows]
    pair = lab["tdc_pair"][rows]
    dec = lab["is_decoy"][rows]
    tb = hash_tiebreak(unit, TB_SEED)
    r = double_competition(spec_key, unit, pair, score, dec, seed=TB_SEED, tiebreak=tb)
    return r


def ours_psms(d):
    import pyarrow as pa
    t = pa.concat_tables([pq.read_table(f) for f in sorted(glob.glob(os.path.join(d, "*.psms.parquet")))])
    t = t.filter(pa.compute.equal(t.column("rank"), 1))
    spec_key = np.array([f"{a}:{b}" for a, b in zip(t.column("run").to_pylist(), t.column("spec_idx").to_pylist())])
    return spec_key, t.column("pf_row").to_numpy(), t.column("hyperscore").to_numpy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--ours", required=True)
    ap.add_argument("--sage", required=True)
    ap.add_argument("--entrap-ours", default=None)
    ap.add_argument("--out", default=os.path.join(C.RES, "ids.json"))
    a = ap.parse_args()
    ds = a.dataset
    res = {}
    idx = C.index_dir(ds, "prod")
    l1 = L1Index(idx)
    lab = labels(idx)
    # ours
    sk, rows, hs = ours_psms(a.ours)
    r = paired_tdc(sk, rows, hs, lab, l1)
    acc = r.accepted(ALPHA)
    res["ours_hyperscore"] = {"n_psms": int(len(rows)), "peptides_at_1pct": int(len(acc)),
                              "peptidoforms_unique_seq_at_1pct": int(len({l1.sequence(rows[i]) for i in acc}))}
    # Sage psms -> our rows
    s = SB.read_results(a.sage).to_pandas()
    s = s[s["rank"] == 1]
    key = {}
    pm = [SB.parse_peptide(x, C.VAR) for x in s.peptide]
    need = {(q, m) for q, m in pm}
    for i in range(len(l1)):
        k = (l1.sequence(i), int(l1.mask[i]))
        if k in need:
            key[k] = i
    srows = np.array([key[k] for k in pm])
    sk2 = (s.filename + ":" + s.scannr.astype(str)).values
    assert np.array_equal(lab["is_decoy"][srows], (s.label.values == -1)), "label mapping mismatch"
    r2 = paired_tdc(sk2, srows, s.hyperscore.values, lab, l1)
    acc2 = r2.accepted(ALPHA)
    res["sage_hyperscore"] = {"n_psms": int(len(srows)), "peptides_at_1pct": int(len(acc2))}
    # Sage own peptide-level q (unique peptidoforms among targets with peptide_q <= 1%)
    tg = s[(s.label == 1) & (s.peptide_q <= ALPHA)]
    res["sage_own_q"] = {"peptides_at_1pct": int(tg.peptide.nunique()), "psms_spectrum_q_at_1pct":
                         int(((s.label == 1) & (s.spectrum_q <= ALPHA)).sum()), "score": "Sage LDA discriminant"}
    # entrapment (dev diagnostic)
    if a.entrap_ours:
        ie = C.index_dir(ds, f"E{C.ENTRAP_SEED}")
        l1e = L1Index(ie)
        labe = labels(ie)
        sk3, rows3, hs3 = ours_psms(a.entrap_ours)
        r3 = paired_tdc(sk3, rows3, hs3, labe, l1e)
        accepted = (r3.q <= ALPHA) & ~r3.is_decoy
        tside = ~r3.is_decoy
        wr = r3.row[tside]
        cnt = ENT.paired_counts(labe["entrap_pair"][rows3[wr]], labe["is_entrapment"][rows3[wr]],
                                r3.score[tside], accepted[tside])
        nt, ne, neo, negt = cnt.totals()
        res["ours_hyperscore_entrapment_E2001"] = {
            "targets_plus_entrapments_at_1pct": int(accepted.sum()), "N_T": int(nt), "N_E": int(ne),
            "paired_FDP_hat": float(ENT.paired(nt, ne, neo, negt)), "lower_bound": float(ENT.lower_bound(nt, ne)),
            "note": "single entrapment seed, single decoy seed, no CI: dev diagnostic only (DESIGN §6.3)"}
    out = json.load(open(a.out)) if os.path.exists(a.out) else {}
    out[ds] = res
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({ds: res}, indent=1))


if __name__ == "__main__":
    main()
