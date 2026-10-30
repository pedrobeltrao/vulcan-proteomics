#!/usr/bin/env python3
"""A6.5.2-A6.5.3 self-host mixtures (step 4; LABEL-FREE).  One mixture per admitted run.

  * Host: the null's own run, UNMODIFIED: all MS2 spectra, production hypotheses, searched against the FULL prod index.
    This is the certification search (scripts/a6/admit_hn1.py), reused: its raw PSM table is passed in.
  * Nulls: provenance-flagged COPIES of the run's A spectra (prov 20) and of their FS versions (prov 22, all MS2 peaks
    + 2/7 Th, vulcan.gates.a2_fs.fragment_shift), each with the single certified hypothesis (mono m/z, z, isotope 0;
    vulcan.gates.a2_fs.a2_hypotheses), searched against the project's REDUCED index (scripts/a6/reduced_index.py).
  * The two PSM tables are concatenated (null spec_idx offset by the host size; null pf_row mapped back to full-index
    rows through row_map) BEFORE any feature, calibration or scoring step.  Then, label-free: L3 features
    (vulcan.rescore.features.FeatureBuilder, as scripts/m3/features.py for the A2 mixtures), StaticModel.apply for
    hyperscore, v1 and v1 noRT (scripts/m3e/gates_a1.neutral_rt), psm_competition winners (hash tie-break seed 11).
  * RT is the run's real RT (no quantile mapping).
  * G3 inputs: exhaustive and production top-K of every A copy on the reduced index (scripts/a2/gates_a2.
    search_exhaustive_hyp + scripts/m2/e2_correctness._rescored, as gates_a2.g3), stored as (spectrum, hyperscore,
    full-index row) only.
NOTHING here reads labels.arrow / proteins.arrow, computes a TDC threshold, a tail membership or any target/decoy
split: those are per-label quantities on null-copy rows and exist only at step 5 (scripts/a6/gates_hn1.py; A6 B2).

A6.11 R12 (closed rule; listed before the tag, reported with this direction, never a reason to drop spectra or change
the pipeline): v1 features that depend on the host or the null copy's twin.  NO v1 feature counts host PSMs,
proteins or charges of a candidate peptide (the v1 feature list has no peptide-/protein-count feature; checked
against the frozen model's `names`).  The run-level features below see the host spectra and the null twins:
"""
from __future__ import annotations

import os
import sys
import time
from dataclasses import replace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
for _p in (ROOT, HERE, os.path.join(ROOT, "scripts", "m3"), os.path.join(ROOT, "scripts", "m3e"),
           os.path.join(ROOT, "scripts", "a2"), os.path.join(ROOT, "scripts", "m2")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

P_HOST, P_A, P_FS = 0, 20, 22
TB_SEED = 11
FS_DELTA_TH = 2.0 / 7.0
PIPES = ("fixed", "v1", "v1_noRT")

# A6.11 R12: v1 features that count host PSMs, proteins or charges of a candidate peptide.  Empty by inspection of the
# frozen v1 feature list (static_gbm_pred.names, 31 features).  Direction if any existed: pro-target.
R12_FEATURES: dict = {}
# Run-level v1 features computed over the whole concatenated mixture (host + null copies), hence seeing the host and
# the null twin.  Direction (A6.11 R12): pro-target or neutral (real but unidentified peptides are supported by host
# context; decoys rarely are) -> fail or invalid, never a false pass.  Reported, never used to drop spectra.
R12_RUN_LEVEL_FEATURES = {
    "hs_runz": "median and 90th percentile of rank-1 hyperscore over the mixture (host + copies) standardise each PSM",
    "rt_abs_resid": "RT calibration (rt_calibrate) fitted on the top 10% PSMs of the mixture, incl. host twins",
    "rt_resid_z": "same RT calibration; divided by the mixture's calibration MAD",
    "pred_*": "NCE chosen per mixture (FeatureBuilder.choose_nce) on calibration PSMs incl. host twins; all seven "
              "pred_* fragment features use that NCE",
}


def concat_tables(host: dict, null: dict, n_host: int) -> dict:
    """Concatenate the host PSM table (full index) and the null PSM table (reduced index, pf_row ALREADY mapped back to
    full-index rows): null spec_idx += n_host.  Pure, label-free; columns present in both are kept."""
    cols = [c for c in host if c in null]
    out = {}
    for c in cols:
        a, b = np.asarray(host[c]), np.asarray(null[c])
        if c == "spec_idx":
            b = b.astype(np.int64) + int(n_host)
            a = a.astype(np.int64)
        out[c] = np.concatenate([a, b])
    return out


def winners(spec, pf_id, scores: dict) -> dict:
    """PSM-level competition winners (row indices into the mixture PSM table, sorted by spectrum) per pipeline.
    Label-free: one winner per spectrum over all hypotheses and candidates (psm_competition)."""
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import psm_competition
    tb = hash_tiebreak(pf_id, TB_SEED)
    out = {}
    for p, sc in scores.items():
        w = psm_competition(spec, sc, tiebreak=tb)
        out[p] = w[np.argsort(spec[w], kind="stable")]
    return out


def null_spectra(host_sp, a_spec, a_mz, a_z):
    """A copies (prov 20) then FS copies (prov 22), single certified hypothesis.  -> (Spectra, prov, src store idx)."""
    from g2_nulls import concat_spectra
    from vulcan.gates.a2_fs import certified_precursor, fragment_shift
    a = certified_precursor(host_sp, a_spec, a_mz, a_z)
    fs = fragment_shift(host_sp, a_spec, a_mz, a_z, FS_DELTA_TH)
    sp = concat_spectra([a, fs], host_sp.run + "_hn1null")
    prov = np.r_[np.full(len(a), P_A, np.int8), np.full(len(fs), P_FS, np.int8)]
    src = np.r_[np.asarray(a_spec, np.int64), np.asarray(a_spec, np.int64)]
    return sp, prov, src


def search_nulls(nsp, gir, red, p):  # pragma: no cover (GPU)
    """Production search of the null copies on the reduced index with exactly one hypothesis each."""
    from reduced_index import map_back
    from vulcan.gates.a2_fs import a2_hypotheses
    from vulcan.gpu import engine as E
    from vulcan.gpu.host import pre_hypotheses
    valid, _ = pre_hypotheses(nsp, gir, p)
    hyp = a2_hypotheses(nsp, valid, gir.mass32_np, p, np.ones(len(nsp), bool))
    tbl, st, _ = E.search(nsp, gir, red, p, pre_hyp=(valid, hyp))
    if st.counts.get("pre_hyp_used") != 1:
        raise RuntimeError("A6 null search: single hypotheses were not used")
    return map_back(tbl, red)


def g3_tables(asp, gir, red):  # pragma: no cover (GPU)
    """Exhaustive and production top-K of the A copies on the reduced index (as scripts/a2/gates_a2.g3), label-free:
    (spectrum, hyperscore, FULL-index row) for each."""
    import torch
    from e2_correctness import _rescored
    from gates_a2 import search_exhaustive_hyp
    from vulcan.gates.a2_fs import a2_hypotheses
    from vulcan.gpu import engine as E
    from vulcan.gpu.host import pre_hypotheses
    p0 = E.SearchParams()
    valid, _ = pre_hypotheses(asp, gir, p0)
    hyp = a2_hypotheses(asp, valid, gir.mass32_np, p0, np.ones(len(asp), bool))
    _, _, dx = search_exhaustive_hyp(asp, gir, red, p0, hyp, np.arange(len(asp)))
    sx, hx, rx = _rescored(dx, red, p0)
    _, _, dk = E.search(asp, gir, red, p0, pre_hyp=(valid, hyp))
    sk, hk, rk = _rescored(dk, red, p0)
    torch.cuda.empty_cache()
    m = red.row_map
    return {"exh_spec": np.asarray(sx, np.int64), "exh_hs": np.asarray(hx, float), "exh_row": m[np.asarray(rx, np.int64)],
            "prod_spec": np.asarray(sk, np.int64), "prod_hs": np.asarray(hk, float),
            "prod_row": m[np.asarray(rk, np.int64)]}


class Scorer:  # pragma: no cover (GPU / predictor)
    """L3 features (FeatureBuilder, as scripts/m3/features.py) + frozen v1 StaticModel; label-free."""

    def __init__(self, dataset: str = ""):
        import m3_common as C
        from vulcan.predict.apd import APD, check_allowed
        from vulcan.predict.cache import PredCache
        from vulcan.rescore.features import FeatureBuilder
        from vulcan.rescore.static import StaticModel
        check_allowed(dataset)
        pred = APD()
        caches = {}

        def factory(nce):
            if nce not in caches:
                caches[nce] = PredCache(os.path.join(C.WORK, "predcache"), pred, nce, "QE", C.STATIC, C.VAR)
            return caches[nce]
        self.fb = FeatureBuilder(factory, C.STATIC, C.VAR)
        self.model = StaticModel.load(os.path.join(C.RES, "models", "static_gbm_pred.pkl"))

    def score(self, l1_full, psm: dict, mix_sp):
        from gates_a1 import neutral_rt
        feat, info = self.fb.run(l1_full, psm, mix_sp, nce=None)
        scores = {"fixed": np.asarray(psm["hyperscore"], float), "v1": self.model.apply(psm, feat),
                  "v1_noRT": self.model.apply(psm, neutral_rt(feat))}
        return feat, scores, info


def build_mixture(run, host_sp, host_tbl, a_spec, a_mz, a_z, gir, red, p, l1_full, scorer, out_dir):  # pragma: no cover
    """One self-host mixture: search the null copies, concatenate, featurise, score, competition winners, G3 tables.
    Writes <out_dir>/<run>.hn1mix.npz (label-free) and returns label-free counts + timing."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    from g2_nulls import concat_spectra
    T = {}
    t0 = time.perf_counter()
    nsp, prov_n, src_n = null_spectra(host_sp, a_spec, a_mz, a_z)
    null = search_nulls(nsp, gir, red, p)
    T["null_search"] = time.perf_counter() - t0
    host = {c: host_tbl.column(c).to_numpy() for c in host_tbl.column_names} if hasattr(host_tbl, "column") \
        else dict(host_tbl)
    n_host = len(host_sp)
    psm = concat_tables(host, null, n_host)
    mix_sp = concat_spectra([host_sp, nsp], run + "_hn1mix")
    prov = np.r_[np.full(n_host, P_HOST, np.int8), prov_n]
    src = np.r_[np.arange(n_host, dtype=np.int64), src_n]
    t0 = time.perf_counter()
    feat, scores, finfo = scorer.score(l1_full, psm, mix_sp)
    T["features_scores"] = time.perf_counter() - t0
    win = winners(np.asarray(psm["spec_idx"], np.int64), psm["pf_id"], scores)
    t0 = time.perf_counter()
    na = len(a_spec)
    from g2_nulls import subset
    g3 = g3_tables(subset(nsp, np.arange(na)), gir, red)
    T["g3_searches"] = time.perf_counter() - t0
    os.makedirs(out_dir, exist_ok=True)
    pq.write_table(pa.table({k: np.asarray(v) for k, v in psm.items()}), os.path.join(out_dir, run + ".hn1mix.psms.parquet"))
    pq.write_table(pa.table({k: np.asarray(v) for k, v in feat.items()}), os.path.join(out_dir, run + ".hn1mix.feat.parquet"))
    np.savez(os.path.join(out_dir, run + ".hn1mix.npz"), prov=prov, src_store_idx=src,
             **{f"score_{k}": v for k, v in scores.items()}, **{f"win_{k}": v for k, v in win.items()},
             **{f"g3_{k}": v for k, v in g3.items()})
    n_null_spec_with_psm = int(np.unique(np.asarray(null["spec_idx"])).size) if len(null.get("spec_idx", [])) else 0
    return {"run": run, "n_host": n_host, "n_A_copies": na, "n_FS_copies": int((prov_n == P_FS).sum()),
            "n_null_spectra_with_psm": n_null_spec_with_psm, "nce": finfo.get("nce"), "timing_s": T}
