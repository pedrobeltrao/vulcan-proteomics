#!/usr/bin/env python3
"""M6d fix selection: evaluate quant variants (isotope set, window, charge handling, apex re-centring) on a run subset
of a dev set, with the unchanged §6.8 metrics.  The selection subset of D9 is A_REP1, A_REP2, B_REP1 (D9 links,
alignment and window widths are rebuilt from those 3 runs only); the other three runs are held out.
    scripts/m2b/vpy.sh scripts/m6d/grid.py DATASET SUBSET OUT.tsv      (DATASET d9|d11; SUBSET sel|hold|all)
"""
from __future__ import annotations

import itertools
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import datasets as DS  # noqa: E402
import evalq as E  # noqa: E402
import harness as H  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from vulcan.ms1 import averagine as AV  # noqa: E402

ISO = {"M0-2": [0, 1, 2], "M0-1": [0, 1], "M0": [0]}


def intensity(a, calcmass, iso):
    if iso == "avgproj":
        env = np.stack([AV.envelope(m, 3) for m in calcmass])
        X = a[:, 1:4]
        return (X * env).sum(1) / (env * env).sum(1) * env.sum(1)
    j = [list(H.ISO_OFFS).index(c) for c in ISO[iso]]
    return a[:, j].sum(1)


def recentred(obj, k_sd=5.0, min_min=0.03):
    """Rows whose feature apex (on the consensus scale) is off the precursor's median over the OTHER runs by more than
    max(k_sd x run alignment SD, min_min) get their window centred on the aligned prediction instead."""
    al = obj["align"]
    cons = {}
    for r in obj["runs"]:
        t = obj["tabs"][r]
        v = np.full(len(t), np.nan)
        ok = np.isfinite(t["rt_apex"].to_numpy())
        v[ok] = al.to_consensus(r, t["rt_apex"].to_numpy()[ok])
        cons[r] = pd.Series(v, index=pd.MultiIndex.from_arrays([t["pf_id"], t["charge"]]))
    C = pd.concat(cons, axis=1)
    out, n = {}, 0
    for r in obj["runs"]:
        t = obj["tabs"][r].copy()
        idx = pd.MultiIndex.from_arrays([t["pf_id"], t["charge"]])
        oth = C.drop(columns=r).reindex(idx)
        med = oth.median(axis=1).to_numpy()
        nn = oth.notna().sum(axis=1).to_numpy()
        me = C[r].reindex(idx).to_numpy()
        off = np.abs(me - med) > max(k_sd * al.maps[r][3], min_min)
        off &= nn >= 2
        t.loc[off, "rt_apex"] = al.to_run(r, med[off])
        n += int(off.sum())
        out[r] = t
    return out, n


def run_grid(ds, subset):
    cfg = DS.get(ds)
    runs, A, B = cfg.subset(subset)
    obj = H.prepare(cfg.ids, runs, cfg.ms1, cfg.store, cache=f"/tmp/vulcan_m6d/{ds}_{subset}_obj.pkl")
    rows = []
    rc_tabs, n_rc = recentred(obj)
    sib = H.sibling_rows(obj)
    for wm, rec in itertools.product((1.0, 1.25, 1.5), (False, True)):
        tabs = rc_tabs if rec else obj["tabs"]
        base = H.requant(obj, cfg.store, tabs, width_mult=wm)
        sres = H.requant(obj, cfg.store, sib, width_mult=wm, mz_source="theory")
        for iso in ("M0-2", "M0-1", "avgproj"):
            vals = {r: intensity(base[r][0], obj["tabs"][r]["calcmass"].to_numpy(), iso) for r in runs}
            ext = []
            for r in runs:
                s = sib[r]
                v = intensity(sres[r][0], s["calcmass"].to_numpy(), iso)
                okk = sres[r][1] >= 3          # mono seen in >= 3 MS1 scans inside the window
                ext.append(pd.DataFrame({"run": r, "pf_id": s["pf_id"], "charge": s["charge"], "intensity": np.where(okk, v, np.nan)}))
            ext = pd.concat(ext, ignore_index=True)
            for sibling in (False, True):
                L = E.long_table(obj, vals, ext if sibling else None)
                Mp, Pp = E.matrix(L, "precursor")
                op, _ = E.score(Mp, Pp, A, B, cfg.expected)
                rec_ = {"dataset": ds, "subset": subset, "width_mult": wm, "recentre": rec, "iso": iso, "sibling": sibling,
                        "prec_n": op["n"], "prec_eps": op["eps"], "prec_cv": op.get("cv", op.get("cv2")),
                        "prec_completeness": op["completeness"]}
                for s_, v_ in op["eps_sp"].items():
                    rec_[f"prec_eps_{s_}"] = v_
                for rule in ("sum", "dominant"):
                    Mf, Pf = E.matrix(L, "peptidoform", rule)
                    of, _ = E.score(Mf, Pf, A, B, cfg.expected)
                    rec_[f"pf_{rule}_n"] = of["n"]
                    rec_[f"pf_{rule}_eps"] = of["eps"]
                    rec_[f"pf_{rule}_cv"] = of.get("cv", of.get("cv2"))
                rows.append(rec_)
                print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in rec_.items()}), flush=True)
    T = pd.DataFrame(rows)
    T["n_recentred"] = n_rc
    return T


if __name__ == "__main__":
    T = run_grid(sys.argv[1], sys.argv[2])
    T.to_csv(sys.argv[3], sep="\t", index=False)
