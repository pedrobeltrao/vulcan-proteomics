#!/usr/bin/env python3
"""3c follow-up, item 2: the strongest comparator (A1.9) for the IT-CID class: Sage 0.14.7 PSMs + the SAME Prosit CID /
TMT-CID predicted features our frozen scorer uses (vulcan.rescore.itcid on Sage's own PSMs, every prediction read from
the pinned store) appended to Sage's PIN -> seeded mokapot 0.10.0 (bench/t1/mokapot_seeded.py via
bench.t1.comparators.run_mokapot: --seed 11, canonical PIN row order).  Also seeded mokapot on Sage's plain PIN
("sage+mokapot(seeded)"), since the 3c sage+mokapot run used mokapot's unseeded CLI.

Inputs: the 3c Sage searches ($VULCAN_3C_WORK/comp/sage/<ds>/<db>, M1e DBs: prod, P2001-P2003, C; MGFs from our stores,
TITLE = spectrum index of our store), unchanged.  Features: ITCID_PRED_FEATURES (pred_sa, pred_sa_sqrt, pred_pearson,
pred_entropy, pred_frac_obs, pred_top6, pred_obs_tic, rt_abs_resid, rt_resid_z, pred_missing), label-blind, on every PIN
row (Sage report_psms = 1); RT calibration per run on the top-decile Sage hyperscore PSMs of both labels (as ours).

    python scripts/m3c/sage_pred_itcid.py [--datasets D4 D10cf] [--dbs prod P2001 P2002 P2003 C]
Outputs: $VULCAN_3C_FU/comp/{sage_pred,sage_pred_mokapot,sage_mokapot_seeded}/<ds>/<db>/ ; timings in
results/2026-10-27_3c_itcid/followup/comparator_timings.json
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import types

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402
import pinned_cache as PC  # noqa: E402

while C.ROOT in sys.path:
    sys.path.remove(C.ROOT)
sys.path.insert(0, C.ROOT)                        # before scripts/m3c (its bench.py shadows the bench package)
import types as _types  # noqa: E402
_bench = _types.ModuleType("bench")                # scripts/m3c/bench.py would shadow the namespace package
_bench.__path__ = [os.path.join(C.ROOT, "bench")]
sys.modules.setdefault("bench", _bench)
import numpy as np  # noqa: E402

from bench.t1 import comparators as CMP  # noqa: E402

FU_COMP = os.path.join(PC.FU_WORK, "comp")
SAGE_DIR = os.path.join(C.WORK, "comp", "sage")
DBS = ["prod", "P2001", "P2002", "P2003", "C"]
FEATS = ["pred_sa", "pred_sa_sqrt", "pred_pearson", "pred_entropy", "pred_frac_obs", "pred_top6", "pred_obs_tic",
         "rt_abs_resid", "rt_resid_z", "pred_missing"]
_MOD = re.compile(r"\[([+-]?[\d.]+)\]")


def parse_sage_peptide(p: str, tmt: bool):
    """Sage peptide ('[+229.16293]-PEPM[+15.9949]K[+229.16293]') -> (sequence as written, mask: bit k = Met-ox at
    residue k; bit 63 = TMT N-term (TMT family, as our index)).  Static C/K/N-term brackets are dropped."""
    s = p
    nterm = False
    if s.startswith("["):
        j = s.index("]")
        nterm = abs(float(s[1:j]) - C.TMT) < 1e-3
        s = s[j + 1:].lstrip("-")
    seq, mask, i = [], 0, 0
    while i < len(s):
        c = s[i]
        i += 1
        if i < len(s) and s[i] == "[":
            j = s.index("]", i)
            m = float(s[i + 1:j])
            i = j + 1
            if c == "M" and abs(m - 15.9949) < 1e-2:
                mask |= 1 << len(seq)
        seq.append(c)
    if tmt:
        if not nterm:
            raise ValueError(f"TMT family peptide without N-term TMT: {p}")
        mask |= 1 << 63
    return "".join(seq), mask


def mini_l1(peps):
    """Unique (sequence, mask) -> an L1-like object (residues, seq_off, length, mask) for the feature builder."""
    u = sorted(set(peps))
    pos = {k: i for i, k in enumerate(u)}
    length = np.array([len(s) for s, _ in u], np.int64)
    off = np.zeros(len(u), np.int64)
    off[1:] = np.cumsum(length)[:-1]
    res = np.frombuffer("".join(s for s, _ in u).encode(), np.uint8).copy()
    mask = np.array([m for _, m in u], np.uint64)
    return types.SimpleNamespace(residues=res, seq_off=off, length=length, mask=mask), pos


def stores(ds):
    return {os.path.basename(s)[:-len(".vspec.arrow")]: s for s in C.DEV[ds]["stores"]()}


def sage_pred_pin(ds, db, fb, log):
    import pandas as pd
    from vulcan.gpu.engine import load_store
    sd = os.path.join(SAGE_DIR, ds, db)
    out = os.path.join(FU_COMP, "sage_pred", ds, db)
    info_p = os.path.join(out, "pin_info.json")
    if os.path.exists(info_p):
        return json.load(open(info_p))
    os.makedirs(out, exist_ok=True)
    t0 = time.perf_counter()
    hdr = open(os.path.join(sd, "results.sage.pin")).readline().rstrip("\n").split("\t")
    pin = pd.DataFrame(list(CMP.read_ragged(os.path.join(sd, "results.sage.pin"))))
    tsv = pd.read_csv(os.path.join(sd, "results.sage.tsv"), sep="\t",
                      usecols=["psm_id", "filename", "scannr", "charge", "hyperscore", "expmass"])
    tsv["SpecId"] = tsv.psm_id.astype(str)
    pin = pin.merge(tsv[["SpecId", "scannr", "charge", "hyperscore", "expmass"]], on="SpecId", how="left")
    assert pin.scannr.notna().all(), "PIN/TSV join failed"
    tmt = C.is_tmt(C.DEV[ds]["family"])
    parsed = [parse_sage_peptide(p, tmt) for p in pin.Peptide]
    l1, pos = mini_l1(parsed)
    pf_row = np.array([pos[x] for x in parsed], np.int64)
    pin["run"] = pin.FileName.str.replace(r"\.mgf$", "", regex=True)
    feats = {k: np.zeros(len(pin), np.float32) for k in FEATS}
    info = {"n_psms": int(len(pin)), "n_peptidoforms": int(len(l1.length)), "runs": {}}
    scan = pin.scannr.values.astype(np.int64)
    for run, path in sorted(stores(ds).items()):
        m = np.flatnonzero(pin.run.values == run)
        if not len(m):
            continue
        sp = load_store(path)
        si = scan[m]                                    # MGF TITLE / Sage scannr = store spectrum index (checked)
        assert si.max() < len(sp)
        z = pin.charge.values[m].astype(np.int64)
        dm = np.abs((sp.precursor_mz[si].astype(np.float64) - 1.007276) * z - pin.expmass.values[m].astype(float))
        assert np.median(dm) < 0.02, f"scan mapping check failed for {run}: median |dM| {np.median(dm)}"
        psm = {"pf_row": pf_row[m], "charge": z.astype(np.int8), "spec_idx": si,
               "rank": pin["rank"].values[m].astype(np.int64).astype(np.int8), "hyperscore": pin.hyperscore.values[m].astype(np.float64),
               "rt": sp.rt[si].astype(np.float64)}
        f, finfo = fb.run(l1, psm, sp)
        for k in FEATS:
            feats[k][m] = f[k]
        info["runs"][run] = {"n_psms": int(len(m)), "rt_mad_min": finfo["rt_mad_min"],
                             "frac_pred_missing": finfo["frac_pred_missing"],
                             "median_abs_expmass_diff_da": float(np.median(dm))}
    for k in FEATS:
        pin[k] = feats[k]
    cols = hdr[:-2] + FEATS + hdr[-2:]
    with open(os.path.join(out, "sage_pred.pin.partial"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in pin[cols].itertuples(index=False):
            fh.write("\t".join(str(x) for x in row) + "\n")
    os.replace(os.path.join(out, "sage_pred.pin.partial"), os.path.join(out, "sage_pred.pin"))
    info["wall_s"] = time.perf_counter() - t0
    json.dump(info, open(info_p, "w"), indent=1, default=float)
    log(f"{ds} {db}: sage+pred PIN {info['n_psms']} PSMs in {info['wall_s']:.0f}s")
    return info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=["D10cf", "D4"])
    ap.add_argument("--dbs", nargs="+", default=DBS)
    ap.add_argument("--mode", default="fill", choices=["fill", "frozen"])
    a = ap.parse_args()
    CMP.ensure_mokapot()
    tp = os.path.join(PC.FU_RES, "comparator_timings.json")
    tim = json.load(open(tp)) if os.path.exists(tp) else {}
    for ds in a.datasets:
        fb = PC.builder(C.DEV[ds]["family"], a.mode)
        for db in a.dbs:
            info = sage_pred_pin(ds, db, fb, print)
            mo = os.path.join(FU_COMP, "sage_pred_mokapot", ds, db)
            r1 = CMP.run_mokapot([os.path.join(FU_COMP, "sage_pred", ds, db, "sage_pred.pin")], mo)
            ms = os.path.join(FU_COMP, "sage_mokapot_seeded", ds, db)
            r2 = CMP.run_mokapot([os.path.join(SAGE_DIR, ds, db, "results.sage.pin")], ms)
            tim[f"{ds}/{db}"] = {"features_wall_s": info.get("wall_s"), "mokapot_pred": r1, "mokapot_plain": r2}
            json.dump(tim, open(tp, "w"), indent=1, default=float)
            print(ds, db, "mokapot done", flush=True)


if __name__ == "__main__":
    main()
