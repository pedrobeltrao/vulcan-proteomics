#!/usr/bin/env python3
"""A1.3 G2 null mixtures (DESIGN §6.4 mixture design + §6.0 A1.3), production search.  Dev only.

Null sources (instrument class OT-HCD; >= 3 projects required, 6 used):
  prov 1  offset   D0 MS2 copied K = 2 times with precursor mass + n + 0.5 + U(-0.1, 0.1) Da, n in {-2,-1,1,2}
                   (M3 recipe, seed 11), kept as specified; only in the first two mixtures (one per D0 host) so that
                   no host spectrum is duplicated more than K times
  prov 2  P. furiosus  D1 PXD001077 (1 run)
  prov 3  S. cerevisiae D7 PXD031616 (6 runs; human alpha-synuclein spike handled by the homolog filter)
  prov 4  T. brucei PXD015648, 5  S. aureus PXD023811, 6  Synechococcus PCC 7002 PXD021787, 7  B. cereus PXD055392
          (dev set DN, data/manifests/g2null.tsv, 3 runs each)
Recipe (fixed before searching):
  * every species-null MS2 spectrum is used once: pooled over sources, permuted (numpy default_rng(11)), cut into
    mixtures of N_PER_MIX = 20,000 (a final remainder >= 10,000 forms the last mixture, else is appended to it);
  * mixture k: host = D0 run k % 2 (all its spectra, prov 0) + its species nulls (+ offsets for k < 2);
  * A1.3 RT mapping: each species-null RT is quantile-mapped from its source run's MS2 RT distribution onto the
    host run's MS2 RT distribution (vulcan.gates.nulls.quantile_map) before mixing;
  * production search: prod index (T + pseudo-reverse decoys), E2 SearchParams with report_psms = 5.
Outputs ($VULCAN_M3E_WORK/null_mix/): <mix>.spectra.pkl, <mix>.prov.npy (code per spectrum), <mix>.srcrun.npy
(source run id), mix_info.json; PSMs to psms/G2mix/<mix>.psms.parquet.  Features: scripts/m3/features.py --spectra.

    conda_envs/pred/bin/python scripts/m3e/g2mix.py [--build] [--search]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import pickle
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m3"))
sys.path.insert(0, HERE)
import m3_common as C  # noqa: E402

import numpy as np  # noqa: E402

W = os.environ.get("VULCAN_M3E_WORK", "/tmp/vulcan_m3e")
MIX = os.path.join(W, "null_mix")
PSM = os.path.join(W, "psms", "G2mix")
N_PER_MIX = 20000
K_COPIES = 2
SEED = 11
D7 = ["LM2480_DDA", "LM2481_DDA", "LM2486_DDA", "LM2487_DDA", "LM2490_DDA", "LM2492_DDA"]
DN_PROV = {"PXD015648": 4, "PXD023811": 5, "PXD021787": 6, "PXD055392": 7}
PROV_NAME = {0: "real", 1: "offset", 2: "pfu", 3: "yeast", 4: "tbrucei", 5: "saureus", 6: "synechococcus",
             7: "bcereus"}


def sources():
    """-> list of (prov, run_name, store_path)."""
    import csv
    out = [(2, "Velos005137", C.E2.store_path(C.DEV["D1"]["raw"][0]))]
    out += [(3, r, os.path.join(C.WORK, "store_null", r + ".vspec.arrow")) for r in D7]
    for r in csv.DictReader(open(os.path.join(C.ROOT, "data", "manifests", "g2null.tsv")), delimiter="\t"):
        b = r["filename"][:-4]
        out.append((DN_PROV[r["accession"]], b, os.path.join(W, "store_dn", b + ".vspec.arrow")))
    return out


def build():
    from g2_nulls import concat_spectra, subset
    from vulcan.gates.nulls import quantile_map
    from vulcan.gpu import engine as E
    os.makedirs(MIX, exist_ok=True)
    hosts = [E.load_store(C.E2.store_path(r)) for r in C.DEV["D0"]["raw"]]
    src = sources()
    parts, prov, srcrun, rt_map = [], [], [], [[], []]
    info = {"sources": []}
    for k, (pv, name, path) in enumerate(src):
        sp = E.load_store(path)
        parts.append(sp)
        prov.append(np.full(len(sp), pv, np.int8))
        srcrun.append(np.full(len(sp), k, np.int16))
        for h, host in enumerate(hosts):
            rt_map[h].append(quantile_map(sp.rt, host.rt))
        info["sources"].append({"id": k, "prov": pv, "prov_name": PROV_NAME[pv], "run": name, "n_spectra": len(sp),
                                "rt_q01_50_99": np.quantile(sp.rt, [0.01, 0.5, 0.99]).tolist()})
        print(name, len(sp), flush=True)
    pool = concat_spectra(parts, "pool")
    prov, srcrun = np.concatenate(prov), np.concatenate(srcrun)
    rt_map = [np.concatenate(x) for x in rt_map]
    n = len(pool)
    perm = np.random.default_rng(SEED).permutation(n)
    nmix = max(1, n // N_PER_MIX)
    cuts = [i * N_PER_MIX for i in range(nmix)] + [n]
    if n - nmix * N_PER_MIX >= N_PER_MIX // 2:
        cuts = [i * N_PER_MIX for i in range(nmix + 1)] + [n]
    chunks = [perm[cuts[i]:cuts[i + 1]] for i in range(len(cuts) - 1)]
    rng = np.random.default_rng(SEED)
    info.update(n_species_null=int(n), n_mixtures=len(chunks), n_per_mix=N_PER_MIX, mixtures={})
    for k, idx in enumerate(chunks):
        h = k % 2
        host = hosts[h]
        idx = np.sort(idx)
        ps, pv, sr = [host], [np.zeros(len(host), np.int8)], [np.full(len(host), -1, np.int16)]
        if k < 2:
            known = np.flatnonzero(host.charge > 0)
            for _ in range(K_COPIES):
                nn = rng.choice([-2, -1, 1, 2], size=len(known))
                delta = nn + 0.5 + rng.uniform(-0.1, 0.1, size=len(known))
                z = host.charge[known].astype(np.float64)
                ps.append(subset(host, known, host.precursor_mz[known] + delta / z))
                pv.append(np.ones(len(known), np.int8))
                sr.append(np.full(len(known), -1, np.int16))
        nul = subset(pool, idx)
        nul.rt = rt_map[h][idx].astype(nul.rt.dtype)
        ps.append(nul)
        pv.append(prov[idx])
        sr.append(srcrun[idx])
        name = f"g2mix{k:02d}_h{h}"
        mix = concat_spectra(ps, name)
        with open(os.path.join(MIX, name + ".spectra.pkl"), "wb") as fh:
            pickle.dump(mix, fh, protocol=5)
        np.save(os.path.join(MIX, name + ".prov.npy"), np.concatenate(pv))
        np.save(os.path.join(MIX, name + ".srcrun.npy"), np.concatenate(sr))
        info["mixtures"][name] = {"host": host.run, "n_spectra": len(mix),
                                  "n_by_prov": {PROV_NAME[int(u)]: int(c) for u, c in
                                                zip(*np.unique(np.concatenate(pv), return_counts=True))}}
        print(name, info["mixtures"][name], flush=True)
    json.dump(info, open(os.path.join(MIX, "mix_info.json"), "w"), indent=1)


def search():
    import pyarrow.parquet as pq
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    os.makedirs(PSM, exist_ok=True)
    p = E.SearchParams(report_psms=C.REPORT)
    l1 = L1Index(C.index_dir("human_contam", "prod"))
    gi = E.GpuIndex(l1, p, C.STATIC, C.VAR)
    for f in sorted(glob.glob(os.path.join(MIX, "*.spectra.pkl"))):
        name = os.path.basename(f)[:-len(".spectra.pkl")]
        out = os.path.join(PSM, name + ".psms.parquet")
        if os.path.exists(out):
            continue
        mix = pickle.load(open(f, "rb"))
        t0 = time.perf_counter()
        tbl, st, _ = E.search(mix, gi, l1, p)
        pq.write_table(tbl, out)
        print(name, len(mix), tbl.num_rows, round(time.perf_counter() - t0, 1), "s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--search", action="store_true")
    a = ap.parse_args()
    if a.build:
        build()
    if a.search:
        search()


if __name__ == "__main__":
    main()
