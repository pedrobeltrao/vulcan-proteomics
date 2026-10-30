#!/usr/bin/env python3
"""A6.5.1 reduced index: the prod index with the rows of R_P removed (row subset of the frozen L1Index arrays).

LABEL-FREE: this module never opens labels.arrow / proteins.arrow.  It takes an L1Index (label-free by construction,
vulcan.index.build.L1Index) and a boolean `remove` mask over its rows, and returns an object with the L1Index
interface (pf_id, mass32, seq_off, length, mask, n_var, missed, partner, partners, residues, meta, sequence(),
__len__) plus `row_map` (reduced row -> full-index row).  `residues` is shared; `partner` is remapped to reduced rows
(R_P is pair-closed, so a kept row's partner is kept; asserted).  vulcan.gpu.engine.GpuIndex, search and
search_exhaustive are used UNCHANGED on it.

Exactness test (A6.5.1; `exactness_test`): with R = {} the PSM tables are bit-identical to the full index (after
mapping pf_row back); with a non-empty R every returned row maps outside R and has the same hyperscore as the same
(spectrum, row) on the full index (exhaustive search).

    conda_envs/pe/bin/python scripts/a6/reduced_index.py --exactness  # dev D0 run; writes results/.../exactness.json
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)


class ReducedL1:
    """Row subset of an L1Index (same interface as vulcan.index.build.L1Index)."""

    def __init__(self, l1, remove: np.ndarray):
        remove = np.asarray(remove, bool)
        assert remove.shape == (len(l1),)
        keep = np.flatnonzero(~remove)
        self.full = l1
        self.path = getattr(l1, "path", None)
        self.row_map = keep.astype(np.int64)
        inv = np.full(len(l1), -1, np.int64)
        inv[keep] = np.arange(keep.size)
        self.pf_id = l1.pf_id[keep]
        self.mass32 = l1.mass32[keep]
        self.seq_off = l1.seq_off[keep]
        self.length = l1.length[keep]
        self.mask = l1.mask[keep]
        self.n_var = l1.n_var[keep]
        self.missed = l1.missed[keep]
        par = l1.partner[keep].astype(np.int64)
        has = par >= 0
        newp = np.full(keep.size, -1, np.int64)
        newp[has] = inv[par[has]]
        if np.any(newp[has] < 0):
            raise ValueError("removal set is not pair-closed: a kept row's partner was removed")
        self.partner = newp.astype(l1.partner.dtype)
        self.partners = [self.partner]
        self.residues = l1.residues
        self.meta = dict(l1.meta, a6_reduced={"n_removed": int(remove.sum()), "n_kept": int(keep.size)})
        if np.any(np.diff(self.mass32.astype(np.float64)) < 0) and not np.any(np.diff(l1.mass32.astype(np.float64)) < 0):
            raise AssertionError("row subset broke mass order")

    def __len__(self):
        return len(self.pf_id)

    def sequence(self, i: int) -> str:
        o = int(self.seq_off[i])
        return self.residues[o:o + int(self.length[i])].tobytes().decode()


def pair_close(l1, rows_flag: np.ndarray) -> np.ndarray:
    """Close a row flag under the label-free pair map (l1.partner): flagged iff the row or its partner is flagged."""
    f = np.asarray(rows_flag, bool).copy()
    par = l1.partner.astype(np.int64)
    has = par >= 0
    f[has] |= rows_flag[par[has]]
    return f


def map_back(tbl, red: ReducedL1):
    """PSM table from a reduced-index search -> numpy dict with pf_row mapped to full-index rows."""
    out = {c: tbl.column(c).to_numpy() for c in tbl.column_names}
    out["pf_row"] = red.row_map[out["pf_row"].astype(np.int64)]
    return out


def exactness_test(store_path: str, n_spec: int = 3000, frac_remove: float = 0.3, seed: int = 5) -> dict:
    """A6.5.1 exactness test on a dev store (D0).  Single-hypothesis sub-spectra as in the null searches."""
    from dataclasses import replace

    import m3_common as C
    from vulcan.gates.a2_fs import certified_precursor
    from vulcan.gpu import engine as E
    from vulcan.gpu.host import NEUTRON  # noqa: F401
    from vulcan.index.build import L1Index
    idx = "/tmp/vulcan_e2/index/human_contam_prod"
    t0 = time.perf_counter()
    l1 = L1Index(idx)
    p = E.SearchParams(report_psms=C.REPORT)
    gi = E.GpuIndex(l1, p, C.STATIC, C.VAR)
    sp = E.load_store(store_path)
    tbl, _, _ = E.search(sp, gi, l1, p)
    r1 = tbl.column("rank").to_numpy() == 1
    s = tbl.column("spec_idx").to_numpy()[r1][:n_spec].astype(np.int64)
    z = tbl.column("charge").to_numpy()[r1][:n_spec].astype(np.int64)
    iso = tbl.column("isotope_error").to_numpy()[r1][:n_spec].astype(np.int64)
    mz = sp.precursor_mz[s] - iso * 1.00335483 / z
    sub = certified_precursor(sp, s, mz, z)
    p0 = replace(p, isotope_errors=(0, 0))
    res = {"store": store_path, "n_spec": int(len(s))}

    def cols(t):
        return {c: t.column(c).to_numpy() for c in t.column_names}

    full_prod = cols(E.search(sub, gi, l1, p0)[0])
    full_exh = cols(E.search_exhaustive(sub, gi, l1, p0, report=1 << 30)[0])
    # (1) R = {} : bit-identical
    red0 = ReducedL1(l1, np.zeros(len(l1), bool))
    gi0 = E.GpuIndex(red0, p, C.STATIC, C.VAR)
    a = map_back(E.search(sub, gi0, red0, p0)[0], red0)
    same = set(a) == set(full_prod) and all(np.array_equal(a[c], full_prod[c]) for c in a if c != "pf_row") \
        and np.array_equal(a["pf_row"], full_prod["pf_row"])
    res["empty_R_prod_bit_identical"] = bool(same)
    del gi0
    # (2) R = random pair-closed subset of the rows in the windows
    rng = np.random.default_rng(seed)
    cand = np.unique(full_exh["pf_row"].astype(np.int64))
    flag = np.zeros(len(l1), bool)
    flag[rng.choice(cand, int(frac_remove * cand.size), replace=False)] = True
    flag = pair_close(l1, flag)
    red = ReducedL1(l1, flag)
    t1 = time.perf_counter()
    gir = E.GpuIndex(red, p, C.STATIC, C.VAR)
    res["reduced_gpuindex_build_s"] = time.perf_counter() - t1
    ex = map_back(E.search_exhaustive(sub, gir, red, p0, report=1 << 30)[0], red)
    pr = map_back(E.search(sub, gir, red, p0)[0], red)
    res["R_rows"] = int(flag.sum())
    res["exh_rows_outside_R"] = bool(not flag[ex["pf_row"].astype(np.int64)].any())
    res["prod_rows_outside_R"] = bool(not flag[pr["pf_row"].astype(np.int64)].any())
    key_f = dict(zip(zip(full_exh["spec_idx"].tolist(), full_exh["pf_row"].tolist()),
                     full_exh["hyperscore"].tolist()))
    hs_ok = [key_f.get((a_, b_)) == h_ for a_, b_, h_ in zip(ex["spec_idx"].tolist(), ex["pf_row"].tolist(),
                                                           ex["hyperscore"].tolist())]
    res["exh_hyperscore_equal_full"] = bool(all(hs_ok))
    # production search: compare with the full-index PRODUCTION search where the same (spectrum, row) was reported,
    # else with the exhaustive table (the K1 shortlist and the >= 4-peak report filter differ between the two paths)
    key_p = dict(zip(zip(full_prod["spec_idx"].tolist(), full_prod["pf_row"].tolist()),
                     full_prod["hyperscore"].tolist()))
    hs_ok2 = [key_p.get((a_, b_), key_f.get((a_, b_), h_)) == h_
              for a_, b_, h_ in zip(pr["spec_idx"].tolist(), pr["pf_row"].tolist(), pr["hyperscore"].tolist())]
    res["prod_keys_found_in_full"] = float(np.mean([(a_, b_) in key_p or (a_, b_) in key_f
                                                    for a_, b_ in zip(pr["spec_idx"].tolist(),
                                                                      pr["pf_row"].tolist())])) if len(pr["spec_idx"]) else 1.0
    res["prod_hyperscore_equal_full"] = bool(all(hs_ok2))
    # the exhaustive reduced table must be exactly the full exhaustive table minus R rows
    keepf = ~flag[full_exh["pf_row"].astype(np.int64)]
    fs = set(zip(full_exh["spec_idx"][keepf].tolist(), full_exh["pf_row"][keepf].tolist()))
    rs = set(zip(ex["spec_idx"].tolist(), ex["pf_row"].tolist()))
    res["exh_equals_full_minus_R"] = fs == rs
    res["pass"] = all(res[k] for k in ("empty_R_prod_bit_identical", "exh_rows_outside_R", "prod_rows_outside_R",
                                       "exh_hyperscore_equal_full", "prod_hyperscore_equal_full",
                                       "exh_equals_full_minus_R"))
    res["wall_s"] = time.perf_counter() - t0
    return res


if __name__ == "__main__":
    import argparse
    import glob
    sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))
    ap = argparse.ArgumentParser()
    ap.add_argument("--exactness", action="store_true")
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "2026-10-28_A6_hn1", "reduced_index_exactness.json"))
    a = ap.parse_args()
    if a.exactness:
        st = sorted(glob.glob("/tmp/vulcan_e2/store/b19*_293T_proteinID_*.vspec.arrow"))[0]
        r = exactness_test(st)
        json.dump(r, open(a.out, "w"), indent=1, default=float)
        print(json.dumps(r, indent=1, default=float))
        sys.exit(0 if r["pass"] else 1)
