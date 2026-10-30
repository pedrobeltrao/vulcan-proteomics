#!/usr/bin/env python3
"""G2 null-spectra generator + mixture search (DESIGN §6.4 G2, mixture design).

Null provenance (two sources, §6.4):
  offset   every D0 MS2 spectrum with a known precursor charge, copied K=2 times with the precursor MASS shifted by
           n + 0.5 + U(-0.1, 0.1) Da, n uniform in {-2,-1,1,2} (seed 11): non-physical mass defects, no true match.
           (D0 has no pilot recalibration yet; the shift is applied to the store's precursor m/z.)
  species  wrong-species spectra: D1 (P. furiosus, prov 2) and 3 runs of D7 (PXD031616 LM2480/2481/2486, yeast
           lysate, Fusion Lumos OT-HCD, prov 3; human alpha-synuclein spike excluded by the filter) MS2 searched
           against the human+contaminant DB.  Homolog filter
           (applied in g2_analyze.py, label-symmetric on the pair's target sequence): drop a null spectrum if any
           candidate's target peptide is I/L-identical, isobaric (N<->GG, Q<->AG) or a single-substitution
           neighbour of a P. furiosus or contaminant tryptic peptide.
Mixture: each D0 run + its offset copies + half of the Pfu spectra form one run (<run>_g2mix), searched with the
production settings (prod index, top-5), then L3 features (scripts/m3/features.py --spectra) -> run-level steps
(hs_runz, NCE choice, RT calibration) see the mixture, as in production.

    conda_envs/pred/bin/python scripts/m3/g2_nulls.py
"""
from __future__ import annotations

import json
import os
import pickle
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402

import numpy as np  # noqa: E402

K_COPIES = 2
SEED = 11
NULL_DIR = os.path.join(C.WORK, "store_null")
YEAST_RUNS = ["LM2480_DDA", "LM2481_DDA", "LM2486_DDA"]   # ingested by vulcan.io.convert into NULL_DIR


def concat_spectra(parts, run):
    from vulcan.gpu.engine import Spectra
    offs, base = [], 0
    for p in parts:
        offs.append(p.offsets[:-1] + base)
        base += p.offsets[-1]
    return Spectra(run, np.concatenate([p.scan for p in parts]), np.concatenate([p.index for p in parts]),
                   np.concatenate([p.rt for p in parts]), np.concatenate([p.precursor_mz for p in parts]),
                   np.concatenate([p.charge for p in parts]), np.concatenate([p.mz for p in parts]),
                   np.concatenate([p.intensity for p in parts]), np.r_[np.concatenate(offs), base].astype(np.int64))


def subset(sp, idx, pmz=None):
    from vulcan.gpu.engine import Spectra
    lens = np.diff(sp.offsets)[idx]
    starts = sp.offsets[:-1][idx]
    g = np.repeat(starts, lens) + (np.arange(lens.sum()) - np.repeat(np.cumsum(lens) - lens, lens))
    off = np.r_[0, np.cumsum(lens)].astype(np.int64)
    return Spectra(sp.run, sp.scan[idx], sp.index[idx], sp.rt[idx],
                   sp.precursor_mz[idx] if pmz is None else pmz, sp.charge[idx], sp.mz[g], sp.intensity[g], off)


def main():
    import pyarrow.parquet as pq
    import torch
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    os.makedirs(NULL_DIR, exist_ok=True)
    rng = np.random.default_rng(SEED)
    d0 = [E.load_store(C.E2.store_path(r)) for r in C.DEV["D0"]["raw"]]
    pfu = E.load_store(C.E2.store_path(C.DEV["D1"]["raw"][0]))
    yeast = [E.load_store(os.path.join(NULL_DIR, f + ".vspec.arrow")) for f in YEAST_RUNS]
    foreign = concat_spectra([pfu] + yeast, "foreign")
    fprov = np.concatenate([np.full(len(pfu), 2, np.int8)] + [np.full(len(y), 3, np.int8) for y in yeast])
    halves = np.array_split(rng.permutation(len(foreign)), len(d0))
    mixes, prov_all = [], {}
    for sp, half in zip(d0, halves):
        parts, prov = [sp], [np.zeros(len(sp), np.int8)]       # 0 real, 1 offset null, 2 species null
        known = np.flatnonzero(sp.charge > 0)
        for _ in range(K_COPIES):
            n = rng.choice([-2, -1, 1, 2], size=len(known))
            delta = n + 0.5 + rng.uniform(-0.1, 0.1, size=len(known))
            z = sp.charge[known].astype(np.float64)
            parts.append(subset(sp, known, sp.precursor_mz[known] + delta / z))
            prov.append(np.ones(len(known), np.int8))
        half = np.sort(half)
        parts.append(subset(foreign, half))
        prov.append(fprov[half])
        mix = concat_spectra(parts, sp.run + "_g2mix")
        mixes.append(mix)
        prov_all[mix.run] = np.concatenate(prov)
        with open(os.path.join(NULL_DIR, mix.run + ".spectra.pkl"), "wb") as fh:
            pickle.dump(mix, fh, protocol=5)
        np.save(os.path.join(NULL_DIR, mix.run + ".provenance.npy"), prov_all[mix.run])
    # production search of the mixtures
    p = E.SearchParams(report_psms=C.REPORT)
    l1 = L1Index(C.index_dir("human_contam", "prod"))
    gi = E.GpuIndex(l1, p, C.STATIC, C.VAR)
    out = C.psm_dir("G2_mix")
    os.makedirs(out, exist_ok=True)
    info = {}
    for mix in mixes:
        t0 = time.perf_counter()
        tbl, st, _ = E.search(mix, gi, l1, p)
        pq.write_table(tbl, os.path.join(out, mix.run + ".psms.parquet"))
        pv = prov_all[mix.run]
        info[mix.run] = {"n_spectra": len(mix), "n_real": int((pv == 0).sum()), "n_offset_null": int((pv == 1).sum()),
                         "n_species_null_pfu": int((pv == 2).sum()), "n_species_null_yeast": int((pv == 3).sum()), "search_s": time.perf_counter() - t0,
                         "n_psms": tbl.num_rows}
    json.dump({"k_copies": K_COPIES, "seed": SEED, "runs": info}, open(os.path.join(out, "g2_info.json"), "w"), indent=1)
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
