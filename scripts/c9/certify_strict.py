#!/usr/bin/env python3
"""C9 pilot: label-free yield cost of a stricter certification (PSM-level q <= 1e-4 instead of 1e-3), per DEV run.

Same certification search as scripts/c9/pilot_yield.py (full prod index; labels read only for certification).
Writes results/2026-10-28_C9_pilot/certify_strict.json: per run, certified counts at q_psm 1e-3 and 1e-4 (both with
peptidoform q <= 0.01, human, not contaminant-near).  No null search, no null candidate, no per-label null quantity.
"""
from __future__ import annotations

import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pilot_yield as PY  # noqa: E402


def main():
    import m3_common as C
    from vulcan.gates import a2_admission as AD
    from vulcan.gates.homolog import HomologFilter, read_fasta_seqs
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    from vulcan.stats._util import hash_tiebreak
    l1 = L1Index(PY.IDX)
    p = E.SearchParams(report_psms=C.REPORT)
    gi = E.GpuIndex(l1, p, C.STATIC, C.VAR)
    lab, src = PY.full_labels(PY.IDX)
    near = HomologFilter(read_fasta_seqs([PY.CONTAM_FASTA])).near
    out = {}
    for path in PY.D0_STORES + PY.D8H_STORES:
        sp = E.load_store(path)
        tbl, _, _ = E.search(sp, gi, l1, p)
        spec = tbl.column("spec_idx").to_numpy().astype(np.int64)
        rows = tbl.column("pf_row").to_numpy().astype(np.int64)
        pf = tbl.column("pf_id").to_numpy()
        hs = tbl.column("hyperscore").to_numpy().astype(float)
        tb = hash_tiebreak(pf, PY.TB_SEED)
        useq = {int(r): l1.sequence(int(r)) for r in np.unique(rows)}
        il = np.array([useq[int(r)].replace("I", "L") for r in rows], dtype=object)
        rec = {}
        for q in (1e-3, 1e-4):
            res = AD.certify_run(spec, pf, lab["tdc_pair"][rows], hs, lab["is_decoy"][rows], il, src, near, tb, q,
                                 PY.Q_PEP)
            rec[f"q_psm_{q:g}"] = res["n_certified"]
        out[os.path.basename(path)[:-len(".vspec.arrow")]] = rec
        print(json.dumps({os.path.basename(path)[:30]: rec}), flush=True)
    tot = {k: sum(v[k] for v in out.values()) for k in ("q_psm_0.001", "q_psm_0.0001")}
    json.dump({"per_run": out, "totals": tot, "retention_1e-4_vs_1e-3": tot["q_psm_0.0001"] / tot["q_psm_0.001"]},
              open(os.path.join(PY.ROOT, "results", "2026-10-28_C9_pilot", "certify_strict.json"), "w"), indent=1)
    print(tot)


if __name__ == "__main__":
    main()
