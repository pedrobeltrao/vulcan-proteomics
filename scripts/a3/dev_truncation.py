#!/usr/bin/env python3
"""A3.1.9 input for V12's n_min: the DEV top-5 truncation rate of the exhaustive v2-best candidate.  Dev only.

V12 (A3.2.4) asks, per A spectrum, whether the exhaustive v2-best target / decoy of the single-hypothesis window W is
inside the production top-5 (hyperscore-ranked); n_min "from A3.1.9 by the V6 rule" = the smallest n per label per bin
at which a symmetric process at the DEV rate passes the 90% Newcombe TOST in [-0.01, +0.01] with probability >= 0.90.
This script measures that dev rate on A-like dev nulls, with the same scoring path step 4b uses for V12
(scripts/a3/step4b_chimeric.exhaustive_v12_scores):
  * a dev MIXTURE built like a DN2 mixture (A2.9.3 / A1): the D0 host run b1906 (production hypotheses) + a seeded
    subsample of 12,000 charge-known MS2 spectra of the D7 yeast null run LM2480_DDA (numpy default_rng(47); an M3
    G2 species-null source; yeast is not in DN2 and not in any training set), RT quantile-mapped onto the host
    (vulcan.gates.nulls.quantile_map); the null spectra carry ONE hypothesis (instrument precursor m/z, charge,
    isotope error 0), as A2.2(e) for A spectra;
  * index: the human prod gate index (E2 human_contam_prod), as DN2;
  * exhaustive (null spectra only): every candidate (report = all), v2 features calibrated on the mixture's
    production primary PSMs with the production NCE (independent review M1; step4b.exhaustive_v12_scores), frozen v2;
  * production: vulcan.gpu.engine.search with the same hypotheses, report_psms 5 (K1 top-K 500 + top-5);
  * per spectrum and label: v2-best candidate of the exhaustive set; "kept" iff its peptidoform row is among the
    production top-5 rows of that spectrum.  Rates overall and per decile of the exhaustive v2-best score
    max(best T, best D) (A2.8 binning).  Labels are read (dev data).
Output: results/2026-10-28_A3_preconditions/truncation_dev.json
    /tmp/vulcan_env_<uid>/pred/bin/python scripts/a3/dev_truncation.py
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))

import numpy as np  # noqa: E402

STORE = os.path.join(os.environ.get("VULCAN_M3_STORE_NULL", "/tmp/vulcan_m3/store_null"), "LM2480_DDA.vspec.arrow")
INDEX = "/tmp/vulcan_e2/index/human_contam_prod"
N_SPEC = 12_000
SEED = 47


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "2026-10-28_A3_preconditions", "truncation_dev.json"))
    a = ap.parse_args()
    import m3_data as MD
    import step4b_chimeric as S4
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    t0 = time.perf_counter()
    from vulcan.gates.nulls import quantile_map
    import m3_common as C
    host = E.load_store(C.E2.store_path(C.DEV["D0"]["raw"][0]))
    sp = E.load_store(STORE)
    known = np.flatnonzero(sp.charge > 0)
    pick = np.sort(np.random.default_rng(SEED).choice(known, size=min(N_SPEC, known.size), replace=False))
    sub = S4.subset_spectra(sp, pick)
    sub.rt = quantile_map(sp.rt, host.rt)[pick].astype(sub.rt.dtype)
    from g2_nulls import concat_spectra
    mix = concat_spectra([host, sub], "dev_v12_mix")
    null_idx = np.arange(len(host), len(mix))
    single = np.zeros(len(mix), bool)
    single[null_idx] = True
    l1 = L1Index(INDEX)
    lab = MD.labels(INDEX)
    r = S4.exhaustive_v12_scores(mix, l1, single, null_idx, want_ce_s=False,
                                 predcache=os.path.join(os.environ.get("VULCAN_A3_WORK", "/tmp/vulcan_a3"), "predcache_trunc"))
    r["spec_idx"] = r["spec_idx"] - len(host)
    r["top5"] = {k - len(host): v for k, v in r["top5"].items() if k >= len(host)}
    dec = lab["is_decoy"][r["pf_row"]]
    spec = r["spec_idx"]
    v2 = r["v2"]
    n = len(sub)
    best = {}
    for name, m in (("T", ~dec), ("D", dec)):
        bs = np.full(n, -np.inf)
        brow = np.full(n, -1, np.int64)
        idx = np.flatnonzero(m)
        o = idx[np.lexsort((-v2[idx], spec[idx]))]
        first = np.r_[True, spec[o][1:] != spec[o][:-1]] if o.size else np.zeros(0, bool)
        w = o[first]
        bs[spec[w]] = v2[w]
        brow[spec[w]] = r["pf_row"][w]
        kept = np.array([(int(brow[s]) in r["top5"].get(int(s), ())) if brow[s] >= 0 else False for s in range(n)])
        best[name] = (bs, kept)
    bmax = np.maximum(best["T"][0], best["D"][0])
    fin = np.isfinite(bmax)
    edges = np.quantile(bmax[fin], np.linspace(0, 1, 11))
    b = np.clip(np.searchsorted(edges, bmax, side="right") - 1, 0, 9)
    out = {"design": __doc__.split("Output:")[0].strip(), "n_spectra": n, "n_exhaustive_psms": int(len(v2)),
           "per_label": {}, "per_decile": [], "wall_s": None}
    for name in ("T", "D"):
        bs, kept = best[name]
        h = np.isfinite(bs)
        out["per_label"][name] = {"n": int(h.sum()), "kept": int(kept[h].sum()), "miss_rate": float(1 - kept[h].mean())}
    for i in range(10):
        rec = {"decile": i}
        for name in ("T", "D"):
            bs, kept = best[name]
            h = np.isfinite(bs) & fin & (b == i)
            rec[name] = {"n": int(h.sum()), "miss_rate": float(1 - kept[h].mean()) if h.any() else float("nan")}
        pooled = [best[nm][1][np.isfinite(best[nm][0]) & fin & (b == i)] for nm in ("T", "D")]
        allk = np.concatenate(pooled)
        rec["pooled_miss_rate"] = float(1 - allk.mean()) if allk.size else float("nan")
        out["per_decile"].append(rec)
    pm = [d["pooled_miss_rate"] for d in out["per_decile"] if np.isfinite(d["pooled_miss_rate"])]
    allk = np.concatenate([best[nm][1][np.isfinite(best[nm][0])] for nm in ("T", "D")])
    out["pooled_miss_rate"] = float(1 - allk.mean())
    out["max_variance_decile_miss_rate"] = float(max(pm, key=lambda p: p * (1 - p)))
    out["production_nce"] = r["nce"]
    out["n_host_spectra"] = int(len(host))
    out["wall_s"] = time.perf_counter() - t0
    json.dump(out, open(a.out, "w"), indent=1, default=float)
    print(json.dumps({k: out[k] for k in ("n_spectra", "n_exhaustive_psms", "per_label", "pooled_miss_rate",
                                          "max_variance_decile_miss_rate")}, indent=1))
    print([round(d["pooled_miss_rate"], 4) for d in out["per_decile"]])


if __name__ == "__main__":
    main()
