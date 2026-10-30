#!/usr/bin/env python3
"""Amendment A2.9.3 / A2.9.4 step 4: DN2 null mixtures (as A1, scripts/m3e/g2mix.py) + FS nulls + single-precursor
hypotheses for A and FS spectra, and the production search of every mixture.  No label is read.

Provenance codes (<mix>.prov.npy): 0 host (D0), 1 A1-recipe offset (reported only), 20 stratum-A candidate,
21 B candidate (not A; final (c) decides B), 22 FS (fragment-shifted copy of an A spectrum).
Recipe (config/a2.yaml "mixtures", "fs"):
  * sources = every run of every admitted DN2 project (admission/<acc>/<run>.npz from admit.py), all spectra;
  * A spectra carry the certified precursor (mono m/z, z) and ONE hypothesis (iso 0) (A2.2(e));
  * FS = every A spectrum, all MS2 peaks + 2/7 Th, certified precursor, one hypothesis (A2.4);
  * pooled nulls (A, B candidates, FS) permuted with numpy default_rng(23), cut into mixtures of 20,000 (a remainder
    >= 10,000 is its own mixture, else appended to the last), host = D0 run k % 2; RT quantile-mapped from the source
    run onto the host (FS: from the donor's run);
  * A1-recipe offsets of the host (2 copies, n + 0.5 + U(-0.1, 0.1) Da, n in {-2,-1,1,2}, rng default_rng(23)) in
    the first two mixtures only (reported only, A2.4 rule 1);
  * search: production E2 SearchParams (report_psms 5), prod gate index, hypotheses = vulcan.gates.a2_fs.a2_hypotheses
    (production for host / offsets / B; single for A / FS).
    python scripts/a2/build_mixtures.py --build --search
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
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))

import numpy as np  # noqa: E402

from vulcan.gates import a2_config  # noqa: E402

W = os.environ.get("VULCAN_A2_WORK", "/tmp/vulcan_a2")
MIX = os.path.join(W, "mix")
PSM = os.path.join(W, "psms", "A2mix")
P_HOST, P_OFF, P_A, P_B, P_FS = 0, 1, 20, 21, 22


def cut(n: int, per: int) -> list[tuple[int, int]]:
    """A1 g2mix cut rule."""
    nmix = max(1, n // per)
    cuts = [i * per for i in range(nmix)] + [n]
    if n - nmix * per >= per // 2:
        cuts = [i * per for i in range(nmix + 1)] + [n]
    return [(cuts[i], cuts[i + 1]) for i in range(len(cuts) - 1)]


def build(cfg):  # pragma: no cover (data)
    import m3_common as C
    from g2_nulls import concat_spectra, subset
    from vulcan.gates.a2_fs import certified_precursor, fragment_shift
    from vulcan.gates.nulls import quantile_map
    from vulcan.gpu import engine as E
    seed = cfg["seeds"]["mixture_permutation"]
    per = cfg["mixtures"]["n_per_mixture"]
    log = json.load(open(os.path.join(ROOT, "data", "manifests", "dn2_select.log.json")))
    hosts = [E.load_store(C.E2.store_path(r)) for r in C.DEV["D0"]["raw"]]
    parts, prov, srcrun, single, rt_map, info = [], [], [], [], [[], []], {"sources": []}
    k = 0
    for acc in log["admitted"]:
        for f in sorted(glob.glob(os.path.join(W, "admission", acc, "*.npz"))):
            z = np.load(f)
            sp = E.load_store(str(z["store"]))
            codes = z["codes"]
            a_idx = z["A_spec"]
            b_idx = np.flatnonzero(codes != 1)
            blocks = [(P_A, certified_precursor(sp, a_idx, z["A_mz_mono"], z["A_z"]), a_idx, True),
                      (P_B, subset(sp, b_idx), b_idx, False),
                      (P_FS, fragment_shift(sp, a_idx, z["A_mz_mono"], z["A_z"], cfg["fs"]["delta_th"]), a_idx, True)]
            for pv, part, idx, sgl in blocks:
                parts.append(part)
                prov.append(np.full(len(part), pv, np.int8))
                srcrun.append(np.full(len(part), k, np.int16))
                single.append(np.full(len(part), sgl, bool))
                for h, host in enumerate(hosts):
                    rt_map[h].append(quantile_map(sp.rt, host.rt)[idx])
            info["sources"].append({"id": k, "accession": acc, "run": sp.run, "n_spectra": len(sp)})
            k += 1
    pool = concat_spectra(parts, "pool")
    prov, srcrun, single = np.concatenate(prov), np.concatenate(srcrun), np.concatenate(single)
    rt_map = [np.concatenate(x) for x in rt_map]
    perm = np.random.default_rng(seed).permutation(len(pool))
    rng = np.random.default_rng(cfg["seeds"]["offsets_a1_recipe"])
    os.makedirs(MIX, exist_ok=True)
    info["mixtures"] = {}
    for m, (a, b) in enumerate(cut(len(pool), per)):
        h = m % 2
        host = hosts[h]
        idx = np.sort(perm[a:b])
        ps, pv, sr, sg = [host], [np.zeros(len(host), np.int8)], [np.full(len(host), -1, np.int16)], \
            [np.zeros(len(host), bool)]
        if m < 2:
            known = np.flatnonzero(host.charge > 0)
            for _ in range(2):
                nn = rng.choice([-2, -1, 1, 2], size=len(known))
                delta = nn + 0.5 + rng.uniform(-0.1, 0.1, size=len(known))
                ps.append(subset(host, known, host.precursor_mz[known] + delta / host.charge[known].astype(float)))
                pv.append(np.full(len(known), P_OFF, np.int8))
                sr.append(np.full(len(known), -1, np.int16))
                sg.append(np.zeros(len(known), bool))
        nul = subset(pool, idx)
        nul.rt = rt_map[h][idx].astype(nul.rt.dtype)
        ps.append(nul)
        pv.append(prov[idx])
        sr.append(srcrun[idx])
        sg.append(single[idx])
        name = f"a2mix{m:03d}_h{h}"
        mix = concat_spectra(ps, name)
        with open(os.path.join(MIX, name + ".spectra.pkl"), "wb") as fh:
            pickle.dump(mix, fh, protocol=5)
        np.save(os.path.join(MIX, name + ".prov.npy"), np.concatenate(pv))
        np.save(os.path.join(MIX, name + ".srcrun.npy"), np.concatenate(sr))
        np.save(os.path.join(MIX, name + ".single.npy"), np.concatenate(sg))
        info["mixtures"][name] = {"host": host.run, "n_spectra": len(mix)}
    json.dump(info, open(os.path.join(MIX, "mix_info.json"), "w"), indent=1)


def search(cfg):  # pragma: no cover (GPU)
    import m3_common as C
    import pyarrow.parquet as pq
    from vulcan.gates.a2_fs import a2_hypotheses
    from vulcan.gpu import engine as E
    from vulcan.gpu.host import pre_hypotheses
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
        single = np.load(os.path.join(MIX, name + ".single.npy"))
        valid, _ = pre_hypotheses(mix, gi, p)
        hyp = a2_hypotheses(mix, valid, gi.mass32_np, p, single)
        t0 = time.perf_counter()
        tbl, st, _ = E.search(mix, gi, l1, p, pre_hyp=(valid, hyp))
        if st.counts.get("pre_hyp_used") != 1:
            raise RuntimeError(f"{name}: A2 hypotheses were not used")
        pq.write_table(tbl, out)
        print(name, len(mix), tbl.num_rows, round(time.perf_counter() - t0, 1), "s", flush=True)


def main():  # pragma: no cover
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--search", action="store_true")
    a = ap.parse_args()
    cfg = a2_config.load()
    if a.build:
        build(cfg)
    if a.search:
        search(cfg)


if __name__ == "__main__":
    main()
