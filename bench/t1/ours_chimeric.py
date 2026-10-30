"""T1 driver, chimeric methods (E7/v2; results/2026-10-28_E7_chimeric_v2/config.yaml).  Dev only.

  vulcan_chim_v1  chimeric search (vulcan.chimeric.search) + the FROZEN static_gbm_pred (v1, sha checked)
  vulcan_v2       chimeric search + the v2 model (results/2026-10-28_E7_chimeric_v2/models/static_v2.pkl)
Both use the chimeric competition (vulcan.chimeric.compete): per precursor group winners, CONGA neighbour filter,
paired double competition, Group-walk over chimeric-rank strata {1, 2, >=3}.

Stages (each resumable): ms1 (vulcan.ms1 features + hypotheses per kept store) -> search (every replicate index; the
primary group reuses the closed `ours/<rep>/<run>.psms.parquet`, identical by construction and checked on prod) ->
score (v2 features incl. the v1 features; v1 and v2 scores) -> <work>/ours_chim/<rep>/<run>.chim.scored.parquet.
Labels are never read before the evaluate loader.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time

import numpy as np

from .ours import MODEL, MODEL_SHA256_PREFIX, REPORT, STATIC, VAR, core_s

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
MODEL_V2 = os.path.join(ROOT, "results", "2026-10-28_E7_chimeric_v2", "models", "static_v2.pkl")
PSM_COLS = ["run", "spec_idx", "rank", "pf_row", "pf_id", "charge", "isotope_error", "hyperscore", "delta_next",
            "delta_best", "matched_b", "matched_y", "matched_peaks", "matched_intensity_pct", "average_ppm",
            "delta_mass_ppm", "longest_b", "longest_y", "longest_y_pct", "peptide_len", "missed_cleavages", "n_var",
            "scored_candidates", "poisson", "rt", "expmass", "calcmass", "group", "chim_rank", "n_groups",
            "ms1_feature_id", "ms1_isolation_fraction", "ms1_feature_intensity", "ms1_iso_cos", "ms1_feature_rt_apex",
            "n_features_in_window", "inst_charge"]


def run_ms1(stores, out_dir, log, threads=8) -> dict:
    os.makedirs(out_dir, exist_ok=True)
    info = {}
    for s in stores:
        run = os.path.basename(s)[:-len(".vspec.arrow")]
        js = os.path.join(out_dir, run + ".ms1.json")
        if not os.path.exists(js):
            subprocess.run([sys.executable, "-m", "vulcan.ms1", s, out_dir, "--threads", str(threads)], check=True,
                           cwd=ROOT, env={**os.environ, "PYTHONPATH": ROOT}, stdout=subprocess.DEVNULL)
        j = json.load(open(js))
        info[run] = {"wall_s": j["total_wall_s"], "process_cpu_s": j["total_process_cpu_s"]}
        log(f"ms1 {run}: {info[run]}")
    return info


def search_replicates(stores, reps, ours_root, out_root, ms1_dir, log) -> dict:
    import pyarrow.parquet as pq
    import torch

    from vulcan.chimeric import search as CS
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    p = E.SearchParams(report_psms=REPORT)
    todo = [r for r in reps if not os.path.exists(os.path.join(out_root, r, "search_done.json"))]
    if todo:
        sps, groups = {}, {}
        for s in stores:
            sp = E.load_store(s)
            sps[sp.run] = sp
            groups[sp.run] = CS.precursor_groups(CS.read_hypotheses(os.path.join(ms1_dir, sp.run + ".hypotheses.arrow")), sp)
    for r in todo:
        d = os.path.join(out_root, r)
        os.makedirs(d, exist_ok=True)
        t0, c0 = time.perf_counter(), core_s()
        l1 = L1Index(reps[r]["index"])
        gi = E.GpuIndex(l1, p, STATIC, VAR)
        rec = {"index_build_s": time.perf_counter() - t0, "runs": {}}
        for run, sp in sps.items():
            prim = pq.read_table(os.path.join(ours_root, r, f"{run}.psms.parquet"))
            if r.startswith("prod_") and r.endswith("1001"):
                ref, _, _ = E.search(sp, gi, l1, p)
                rec.setdefault("primary_reuse_identical", {})[run] = bool(ref.equals(prim))
            t1 = time.perf_counter()
            tbl, st = CS.search_chimeric(sp, gi, l1, p, groups[run], primary=prim)
            pq.write_table(tbl, os.path.join(d, f"{run}.chim.psms.parquet"))
            rec["runs"][run] = {"secondary_search_s": time.perf_counter() - t1, "n_psms": tbl.num_rows, **st}
        rec["wall_s"] = time.perf_counter() - t0
        rec["core_s"] = core_s() - c0
        del gi, l1
        torch.cuda.empty_cache()
        json.dump(rec, open(os.path.join(d, "search_done.json"), "w"), indent=1, default=float)
        log(f"chimeric search {r}: {rec['wall_s']:.1f}s")
    return {r: json.load(open(os.path.join(out_root, r, "search_done.json"))) for r in reps}


class ScorerV2:
    def __init__(self, predcache: str, instrument: str, dataset: str, model_v2: str | None = MODEL_V2):
        import hashlib

        from vulcan.predict.apd import APD, check_allowed
        from vulcan.predict.cache import PredCache
        from vulcan.rescore.features_v2 import FeatureBuilderV2
        from vulcan.rescore.static import StaticModel
        check_allowed(dataset)
        sha = hashlib.sha256(open(MODEL, "rb").read()).hexdigest()
        assert sha.startswith(MODEL_SHA256_PREFIX), f"frozen model hash mismatch: {sha}"
        self.v1 = StaticModel.load(MODEL)
        self.v2 = None
        self.v2_sha = None
        if model_v2 and os.path.exists(model_v2):
            self.v2 = StaticModel.load(model_v2)
            self.v2_sha = hashlib.sha256(open(model_v2, "rb").read()).hexdigest()
        self.pred = APD()
        self.caches = {}

        def factory(nce):
            if nce not in self.caches:
                self.caches[nce] = PredCache(predcache, self.pred, nce, instrument, STATIC, VAR)
            return self.caches[nce]
        self.fb = FeatureBuilderV2(factory, STATIC, VAR)


def score_replicates(scorer: ScorerV2, stores, reps, out_root, ms1_dir, log, force_v2=False) -> dict:
    import pyarrow as pa
    import pyarrow.parquet as pq
    import torch

    from vulcan.gpu.engine import load_store
    from vulcan.index.build import L1Index
    from vulcan.rescore.features_v2 import model_psm, read_ms1_features
    done = lambda r: os.path.join(out_root, r, "score_done.json")  # noqa: E731

    def need(r):
        if not os.path.exists(done(r)):
            return True
        return force_v2 and scorer.v2 is not None and json.load(open(done(r))).get("v2_sha256") != scorer.v2_sha
    todo = [r for r in reps if need(r)]
    if todo:
        sps, mf = {}, {}
        for s in stores:
            sp = load_store(s)
            sps[sp.run] = sp
            mf[sp.run] = read_ms1_features(os.path.join(ms1_dir, sp.run + ".features.arrow"))
    for r in todo:
        d = os.path.join(out_root, r)
        t0, c0 = time.perf_counter(), core_s()
        l1 = L1Index(reps[r]["index"])
        rec = {"runs": {}}
        for run, sp in sps.items():
            fp = os.path.join(d, f"{run}.chim.feat.parquet")
            t = pq.read_table(os.path.join(d, f"{run}.chim.psms.parquet"))
            psm = {c: t.column(c).to_numpy(zero_copy_only=False) for c in PSM_COLS}
            t1 = time.perf_counter()
            if os.path.exists(fp):
                ft = pq.read_table(fp)
                feats = {c: ft.column(c).to_numpy() for c in ft.column_names}
                finfo = json.load(open(fp + ".json"))
            else:
                feats, finfo = scorer.fb.run_v2(l1, psm, sp, mf[run])
                pq.write_table(pa.table({k: np.asarray(v, np.float32) for k, v in feats.items()}), fp)
                json.dump(finfo, open(fp + ".json", "w"), default=float)
            t2 = time.perf_counter()
            mp = model_psm(psm)
            s1 = scorer.v1.apply(mp, feats).astype(np.float64)
            cols = {"run": psm["run"], "spec_idx": psm["spec_idx"], "group": psm["group"], "chim_rank": psm["chim_rank"],
                    "rank": psm["rank"], "pf_row": psm["pf_row"], "pf_id": psm["pf_id"], "charge": psm["charge"],
                    "hyperscore": psm["hyperscore"], "score_v1": s1}
            if scorer.v2 is not None:
                cols["score_v2"] = scorer.v2.apply(mp, feats).astype(np.float64)
            pq.write_table(pa.table(cols), os.path.join(d, f"{run}.chim.scored.parquet"))
            rec["runs"][run] = {"features_s": t2 - t1, "apply_s": time.perf_counter() - t2, "nce": finfo["nce"]}
        rec["wall_s"] = time.perf_counter() - t0
        rec["core_s"] = core_s() - c0
        rec["v2_sha256"] = scorer.v2_sha
        json.dump(rec, open(done(r), "w"), indent=1, default=float)
        torch.cuda.empty_cache()
        log(f"chimeric scored {r}: {rec['wall_s']:.1f}s")
    return {r: json.load(open(done(r))) for r in reps}


# ------------------------------------------------------------------------------------------- evaluate loader
class _L1View:
    """IonSets needs seq_off/length/residues/mask; RepLabels has them (label-free columns)."""

    def __init__(self, rl):
        self.seq_off, self.length, self.residues, self.mask = rl.off, rl.ln, rl.res, rl.mask


ION_CACHE = {}


def chim_units(score_dir, rl, runs, score_col, return_detail=False, neighbour=True, qcol="q_gw"):
    """-> peptidoform DataFrame for bench.t1.evaluate (q = Group-walk over chimeric-rank strata) + 'stratum'."""
    import pandas as pd
    import pyarrow.parquet as pq

    from vulcan.chimeric.compete import IonSets, chimeric_competition

    from .evaluate import SEED, _finish
    parts = [pq.read_table(os.path.join(score_dir, f"{r}.chim.scored.parquet")).to_pandas() for r in runs]
    psm = pd.concat(parts, ignore_index=True)
    rix = {r: i for i, r in enumerate(runs)}
    row = psm.pf_row.values
    dec = rl.is_decoy_row[row]
    unit = psm.pf_id.values
    spec = np.array([rix[r] for r in psm.run.values], np.int64) * (1 << 26) + psm.spec_idx.values.astype(np.int64)
    if ION_CACHE.get("rl") is not rl:          # identity (holds a reference, so an id can never be recycled)
        ION_CACHE.clear()
        ION_CACHE.update(rl=rl, ions=IonSets(_L1View(rl)))
    cc = chimeric_competition(spec, psm.group.values, psm.chim_rank.values, unit, rl.tdc_pair[row],
                              psm[score_col].values, dec, row, ION_CACHE["ions"], seed=SEED, neighbour=neighbour)
    res, win = cc["res"], cc["win"]
    rows = row[win]
    sv = cc["survivors"]          # run membership: every surviving group winner (as ours_units: PSM winners)
    bits = pd.DataFrame({"unit": unit[sv], "bit": [1 << rix[r] for r in psm.run.values[sv]]})
    df = pd.DataFrame({"seq": rl.row_seq(rows), "mask": rl.mask[rows].astype(np.uint64), "is_decoy": res.is_decoy,
                       "q": cc[qcol], "score": res.score})
    b = bits.groupby("unit").bit.agg(lambda x: int(np.bitwise_or.reduce(x.values)))
    df["runbits"] = b.reindex(unit[win]).fillna(0).astype(np.int64).values
    df["stratum"] = cc["stratum"]
    df["q_pooled"] = cc["q_pooled"]
    df["q_gw"] = cc["q_gw"]
    df["q_sep"] = cc["q_sep"]
    out = _finish(df, rl, runs)
    if return_detail:
        return out, cc
    return out
