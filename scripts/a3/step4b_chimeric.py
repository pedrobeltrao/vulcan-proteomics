#!/usr/bin/env python3
"""Amendment A3.2.2, step 4b: label-free chimeric search of the DN2 mixtures (between A2 step 4 and step 5).

Run by an agent that does NOT run step 5, from the A3 freeze commit on a clean tree (refuses otherwise; config/a3.yaml
`freeze`).  No label of the gate DB is read here, and no per-label count is printed or stored (A2.9.4 extended to
scripts/a3/ by tests/a3/test_no_label_counts_a3.py).  Inputs: the committed A2 step-4 products under $VULCAN_A2_WORK
(mix/<mix>.{spectra.pkl,prov.npy,srcrun.npy,single.npy}, mix/mix_info.json, psms/A2mix/<mix>.psms.parquet, admission/
<acc>/<run>.npz) and the frozen offsets table.

Phases (each writes under $VULCAN_A3_WORK/step4b and is resumable; --record writes the committed record):
  --ms1      vulcan.ms1 features + hypotheses (E7 defaults) of every admitted DN2 source run (store path from the
             admission npz) and of the two D0 hosts.  A source run whose store has no MS1 hills gets no secondary groups
             (its spectra keep the closed search; counted).
  --asyn     A-syn hypotheses (vulcan.gates.a3_asyn = specification): per A spectrum, ranks 2 and 3, offsets table,
             seed 47, exclusions env / (b) / (d) at W'.  Needs the certified peptide p (sequence, mask, hyperscore) of
             every A spectrum: recomputed with the unchanged A2 own-species search (scripts/a2/own_search.py
             build_own_index -> search_run -> certify_and_markers, as admit.py --swap-test does); the recomputed
             certified precursor of every A spectrum must equal the committed admission npz (A_mz_mono, A_z) exactly,
             else the phase stops.  (d) uses vulcan.gpu.engine.search_exhaustive on a virtual copy of the A spectrum
             carrying (mz', z'), isotope errors 0..2, max hyperscore over every candidate (label-free).
  --search   per mixture: E7 chimeric search (vulcan.chimeric.search.search_chimeric, E2 production SearchParams,
             report_psms 5, A2 prod gate index).  Group 0 = the committed step-4 PSM table (A2 production hypothesis
             set incl. A2.2(e) single hypotheses for A and FS), reused as `primary` (identical by construction; the
             reuse is checked on the first mixture by a fresh search).  Secondary groups: host (prov 0) and B (prov 21)
             spectra get the secondary groups of their SOURCE run's hypothesis table (spectrum matched by store index);
             A1-recipe offset copies (prov 1) get their host spectrum's groups shifted by the copy's mass offset
             (recovered as (mz_mix - mz_host) * z); A (20) and FS (22) get none (A2.2(e) kept).  MS1 context columns
             and feature ids follow scripts/e7/gates_v2.mixture_groups (per-source feature-id offsets).
             A-syn: every A spectrum's synthetic hypotheses as virtual spectra, searched SEPARATELY with isotope errors
             0..1 (E7 sec_isotope_errors), report_psms 5; group = rank - 1, chim_rank = rank; never part of the
             mixture's competition (<mix>.asyn.psms.parquet).
  --score    v2 features (FeatureBuilderV2.run_v2, primary-only calibration on the mixture's REAL primary PSMs) for the
             chimeric table; A-syn rows are featurised by run_v2 on the concatenation (real + A-syn) and the real rows
             must be bit-identical to the real-only run (asserted), so A-syn never enters calibration.  Scores: v2
             (frozen static_v2.pkl), hyperscore, v2 with neutral_rt (noRT); rt_pred (AlphaPeptDeep RT of the candidate,
             the frozen sequence classifier's second input, as gates_a2).  -> <mix>.chim.scored.parquet,
             <mix>.asyn.scored.parquet.
  --v12      A3.2.2 "for V12 only": every A spectrum of every mixture, single hypothesis W (certified mono m/z, z,
             isotope 0 = A2.2(e)), EXHAUSTIVE candidates (search_exhaustive with the given hypotheses, every candidate
             with >= min_matched_peaks, report = all), v2 features calibrated on the MIXTURE's production primary PSMs
             with the production NCE (exhaustive rows appended as non-calibration rows; asserted bit-identical
             production features), v2 score and ce_s logit of every candidate, hyperscore, and whether the
             candidate row is in the production top-5 of that spectrum (step-4 table).  -> <mix>.v12.parquet.
  --record   per-mixture sha256 of every chimeric / A-syn / V12 table, of every source run's MS1 hypothesis table, of
             the A-syn hypothesis tables; label-free unit counts by chimeric rank and pool; A-syn construction counts
             (accepted / env / b / d / none); the P1 swap-test result.  -> --out (committed before step 5).
  --swap-test  P1 on the first mixture: rebuild the gate index with swap_labels=True (the G0 knob), check its
             label-free files byte-identical, re-run the chimeric search of the first mixture and the A-syn
             construction of its A spectra on it; PASS iff every chimeric PSM table column and the A-syn table hash
             are identical.
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import pickle
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "scripts", "a2"))
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))

import numpy as np  # noqa: E402

from vulcan.gates import a3_asyn as AS  # noqa: E402

A2W = os.environ.get("VULCAN_A2_WORK", "/tmp/vulcan_a2")
A3W = os.environ.get("VULCAN_A3_WORK", "/tmp/vulcan_a3")
MIX = os.path.join(A2W, "mix")
PSM2 = os.path.join(A2W, "psms", "A2mix")
OUT = os.path.join(A3W, "step4b")
STATIC, VAR = {"C": 57.021464}, {"M": 15.9949}
P_HOST, P_OFF, P_A, P_B, P_FS = 0, 1, 20, 21, 22
CHIM_COLS = ["group", "chim_rank", "n_groups", "ms1_feature_id", "ms1_isolation_fraction", "ms1_feature_intensity",
             "ms1_iso_cos", "ms1_feature_rt_apex", "n_features_in_window", "inst_charge"]


def sha256_file(p: str) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------------------------------------ pure helpers
def subset_spectra(sp, idx, pmz=None, charge=None):
    """Copy of spectra idx (optionally with new precursor m/z / charge)."""
    from vulcan.gpu.host import Spectra
    idx = np.asarray(idx, np.int64)
    lens = np.diff(sp.offsets)[idx]
    starts = sp.offsets[:-1][idx]
    g = np.repeat(starts, lens) + (np.arange(lens.sum()) - np.repeat(np.cumsum(lens) - lens, lens))
    off = np.r_[0, np.cumsum(lens)].astype(np.int64)
    return Spectra(sp.run, sp.scan[idx].copy(), sp.index[idx].copy(), sp.rt[idx].copy(),
                   np.asarray(sp.precursor_mz[idx] if pmz is None else pmz, np.float64).copy(),
                   np.asarray(sp.charge[idx] if charge is None else charge).astype(sp.charge.dtype),
                   sp.mz[g].copy(), sp.intensity[g].copy(), off)


def a_spec_keys(mix_runs_by_src: dict, srcrun, index):
    """Label-free A-spectrum key '<source run>|<store spectrum index>' (independent of the mixture layout)."""
    return np.array([f"{mix_runs_by_src[int(s)]}|{int(i)}" for s, i in zip(srcrun, index)], dtype=object)


def offset_delta_mass(mz_mix, mz_host, z):
    """A1-recipe offset copies: mass offset (Da) recovered from the mixture precursor (mz_mix) and the host's."""
    return (np.asarray(mz_mix, np.float64) - np.asarray(mz_host, np.float64)) * np.asarray(z, np.float64)


RT_FEATS = ("rt_abs_resid", "rt_resid_z")


def neutral_rt(feat):
    """noRT ablation: RT features replaced by their median (identical to scripts/m3e/gates_a1.neutral_rt, which is
    not imported here because gates_a1 also holds gate statistics; equality pinned by tests/a3/test_a3_step4b_pure)."""
    f = dict(feat)
    for k in RT_FEATS:
        f[k] = np.full_like(feat[k], np.median(feat[k]))
    return f


def unit_counts(prov, chim_rank, spec_idx, group):
    """Label-free counts of query units (spectrum, group) by pool and chimeric rank {1, 2, >=3}."""
    key = np.asarray(spec_idx, np.int64) * 64 + np.asarray(group, np.int64)
    u, first = np.unique(key, return_index=True)
    pv = np.asarray(prov)[np.asarray(spec_idx)[first]]
    rk = np.minimum(np.asarray(chim_rank)[first], 3)
    out = {}
    for code, name in ((P_HOST, "host"), (P_OFF, "offset"), (P_A, "A"), (P_B, "B"), (P_FS, "FS")):
        m = pv == code
        out[name] = {str(k): int((m & (rk == k)).sum()) for k in (1, 2, 3)}
    return out


# ------------------------------------------------------------------------------------------------ GPU helpers
def search_exhaustive_table(sp, gi, l1, p, hyp, report: int = 1 << 30, batch: int = 20_000_000,
                            spec_subset=None):  # pragma: no cover
    """vulcan.gpu.engine.search_exhaustive with GIVEN hypotheses and a full PSM table (every candidate with
    >= min_matched_peaks; report = all), optionally restricted to the spectra spec_subset."""
    import torch
    from vulcan.gpu import engine as E
    st = E.Stats()
    peaks = E.k0_preprocess(sp, p, gi.dev, st)
    o_h, o_r, o_c, scored, matched, d_spec, d_maxfz = E.run_k1(gi, hyp, peaks, p, st)
    h, r = E.window_candidates(hyp, gi.dev, spec_subset)
    spec = d_spec[h].long()
    parts, i, n = [], 0, len(h)
    while i < n:
        j = min(n, i + batch)
        if j < n:
            j = int(torch.searchsorted(spec, spec[j - 1], side="right"))
        parts.append(E.run_k2(gi, h[i:j], r[i:j], d_spec, d_maxfz, peaks, p))
        i = j
    k2 = {k: torch.cat([q[k] for q in parts]) for k in parts[0]}
    c = torch.full((len(h),), -1, dtype=torch.int64, device=gi.dev)
    return E.psm_table(sp, gi, l1, hyp, h, r, c, k2, scored, matched, peaks, p, report)


def _psm_dict(tbl, chim_defaults=True):
    """PSM columns as numpy; chim_defaults adds every missing chimeric / MS1-context column with the no-feature
    values of vulcan.chimeric.search (group 0, rank 1, one group, feature id -1, 0 / NaN context)."""
    psm = {c: tbl.column(c).to_numpy(zero_copy_only=False) for c in tbl.column_names}
    n = len(psm["pf_row"])
    if chim_defaults:
        dflt = {"group": np.zeros(n, np.uint16), "chim_rank": np.ones(n, np.uint16), "n_groups": np.ones(n, np.uint16),
                "ms1_feature_id": np.full(n, -1, np.int64), "ms1_isolation_fraction": np.zeros(n, np.float32),
                "ms1_feature_intensity": np.zeros(n, np.float32), "ms1_iso_cos": np.zeros(n, np.float32),
                "ms1_feature_rt_apex": np.full(n, np.nan, np.float32), "n_features_in_window": np.zeros(n, np.uint16),
                "inst_charge": np.zeros(n, np.int8)}
        for k, v in dflt.items():
            psm.setdefault(k, v)
    return psm


def models(cfg):  # pragma: no cover
    from vulcan.gates import a3_config
    from vulcan.rescore.static import StaticModel
    a3_config.check_models(cfg)
    return StaticModel.load(os.path.join(ROOT, cfg["scorers"]["v2_sep"]["model"]))


def feature_builder(predcache):  # pragma: no cover
    from vulcan.predict.apd import APD
    from vulcan.predict.cache import PredCache
    from vulcan.rescore.features_v2 import FeatureBuilderV2
    pred = APD()
    caches = {}

    def fac(nce):
        if nce not in caches:
            caches[nce] = PredCache(predcache, pred, nce, "QE", STATIC, VAR)
        return caches[nce]
    return FeatureBuilderV2(fac, STATIC, VAR)


def exhaustive_v12_scores(sp, l1, single, null_idx, want_ce_s: bool, predcache: str, v2_model=None,
                          ce_model=None):  # pragma: no cover
    """V12 scoring of every candidate in the single-hypothesis window W of the spectra null_idx of the MIXTURE `sp`.

    Calibration (independent review M1): v2 features are calibrated exactly as production v2_sep calibrates the
    mixture -- on the mixture's production primary PSMs (A2 hypotheses: production for host / B / offsets, single for
    A / FS; K1 top-K + top-5; group 0), with the production NCE.  The exhaustive rows are appended to the production
    table as non-calibration rows (group 1, so they never enter calibration and model_psm moves them out of hs_runz),
    and the production rows' features must be bit-identical to the production-only run (asserted).
    Returns per exhaustive candidate: spec_idx (mixture index), pf_row, hyperscore, v2 (and ce_s), in_top5 (row is
    among the production top-5 rows of that spectrum), plus top5 and the production nce.  Label-free."""
    import torch
    from vulcan.gates.a2_fs import a2_hypotheses
    from vulcan.gpu import engine as E
    from vulcan.gpu.host import pre_hypotheses
    from vulcan.rescore.features_v2 import model_psm
    from vulcan.rescore.static import StaticModel
    null_idx = np.asarray(null_idx, np.int64)
    if not np.asarray(single, bool)[null_idx].all():
        raise ValueError("V12 spectra must carry the single (A2.2(e)) hypothesis")
    p = E.SearchParams(report_psms=5)
    gi = E.GpuIndex(l1, p, STATIC, VAR)
    valid, _ = pre_hypotheses(sp, gi, p)
    hyp = a2_hypotheses(sp, valid, gi.mass32_np, p, single)
    prod, _, _ = E.search(sp, gi, l1, p, pre_hyp=(valid, hyp))
    ps, pr = prod.column("spec_idx").to_numpy(), prod.column("pf_row").to_numpy()
    top5 = {}
    for s_, r_ in zip(ps, pr):
        top5.setdefault(int(s_), set()).add(int(r_))
    st = E.Stats()
    peaks = E.k0_preprocess(sp, p, gi.dev, st)
    hyp_x = a2_hypotheses(sp, peaks["valid"], gi.mass32_np, p, single)
    tbl = search_exhaustive_table(sp, gi, l1, p, hyp_x, spec_subset=null_idx)
    del gi
    torch.cuda.empty_cache()
    prod_psm = _psm_dict(prod)
    exh = _psm_dict(tbl)
    exh["group"] = np.ones(len(exh["pf_row"]), np.uint16)          # never a calibration row
    fb = feature_builder(predcache)
    f0, info = fb.run_v2(l1, prod_psm, sp, None)
    both = {k: np.concatenate([np.asarray(prod_psm[k]), np.asarray(exh[k]).astype(np.asarray(prod_psm[k]).dtype)])
            for k in prod_psm}
    fb2, _ = fb.run_v2(l1, both, sp, None, nce=info["nce"])
    n0 = len(prod_psm["pf_row"])
    if not all(np.array_equal(np.asarray(fb2[k][:n0]).view(np.uint8), np.asarray(f0[k]).view(np.uint8)) for k in f0):
        raise RuntimeError("V12: exhaustive rows changed the production rows' features (calibration leak)")
    if v2_model is None:
        v2_model = StaticModel.load(os.path.join(ROOT, "results", "2026-10-28_E7_chimeric_v2", "models", "static_v2.pkl"))
    v2_all = v2_model.apply(model_psm(both), fb2).astype(np.float64)
    out = {"spec_idx": exh["spec_idx"].astype(np.int64), "pf_row": exh["pf_row"].astype(np.int64),
           "hyperscore": exh["hyperscore"].astype(np.float64), "v2": v2_all[n0:], "top5": top5, "nce": info["nce"]}
    if want_ce_s:
        from vulcan.encoders.score import score_table
        ce_raw, _ = score_table(ce_model, l1, {"pf_row": exh["pf_row"], "charge": exh["charge"],
                                              "spec_idx": exh["spec_idx"]}, sp)
        out["ce_s"] = np.asarray(ce_raw, np.float64)
    out["in_top5"] = np.array([int(r_) in top5.get(int(s_), ()) for s_, r_ in zip(out["spec_idx"], out["pf_row"])])
    return out


# ------------------------------------------------------------------------------------------------ guards
def require_freeze(cfg):  # pragma: no cover
    from vulcan.gates import a3_config
    a3_config.require_freeze_commit(cfg)


# ------------------------------------------------------------------------------------------------ phases (data)
def admitted_sources():  # pragma: no cover
    info = json.load(open(os.path.join(MIX, "mix_info.json")))
    return info, {s["id"]: s for s in info["sources"]}


def phase_ms1():  # pragma: no cover
    info, src = admitted_sources()
    d = os.path.join(OUT, "ms1")
    os.makedirs(d, exist_ok=True)
    stores = {}
    for s in src.values():
        z = np.load(os.path.join(A2W, "admission", s["accession"], s["run"] + ".npz"))
        stores[s["run"]] = str(z["store"])
    import m3_common as C
    for r in C.DEV["D0"]["raw"]:
        stores[os.path.basename(C.E2.store_path(r))[:-len(".vspec.arrow")]] = C.E2.store_path(r)
    rec = {}
    for run, st in sorted(stores.items()):
        if not os.path.exists(st):
            raise SystemExit(f"step4b: store of {run} is not retained ({st}); A3.1.5 not met")
        if not os.path.exists(os.path.join(d, run + ".hypotheses.arrow")):
            subprocess.run([sys.executable, "-m", "vulcan.ms1", st, d, "--threads", "8"], check=True, cwd=ROOT,
                           env={**os.environ, "PYTHONPATH": ROOT}, stdout=subprocess.DEVNULL)
        rec[run] = {"store": st, "hypotheses_sha256": sha256_file(os.path.join(d, run + ".hypotheses.arrow")),
                    "features_sha256": sha256_file(os.path.join(d, run + ".features.arrow"))}
    json.dump(rec, open(os.path.join(OUT, "ms1_record.json"), "w"), indent=1)
    return rec


def certified_A(acc: str, cfg2: dict):  # pragma: no cover
    """Recompute the A2 certification of one admitted project (unchanged own_search code) and return, per run, the
    certified peptide of every A spectrum, checked against the committed admission npz."""
    import own_search as OS
    from vulcan.gpu import engine as E
    ref = json.load(open(os.path.join(A2W, "admission", acc, "counts.json")))
    runs = sorted(ref["runs"])
    npz = {r: np.load(os.path.join(A2W, "admission", acc, r + ".npz")) for r in runs}
    idx, src = OS.build_own_index(ref["upid"], cfg2)
    rr = {r: OS.search_run(str(npz[r]["store"]), idx) for r in runs}
    cm = OS.certify_and_markers(rr, idx, src, cfg2)
    out = {}
    for r in runs:
        c = cm["certified"][r]
        a = npz[r]["A_spec"]
        pos = {int(s): i for i, s in enumerate(c["spec"])}
        if any(int(s) not in pos for s in a):
            raise SystemExit(f"step4b: recomputed certification of {acc}/{r} does not contain every A spectrum")
        ii = np.array([pos[int(s)] for s in a], np.int64)
        if not (np.array_equal(np.asarray(c["mz_mono"])[ii], npz[r]["A_mz_mono"]) and
                np.array_equal(np.asarray(c["z"])[ii], npz[r]["A_z"])):
            raise SystemExit(f"step4b: recomputed certified precursors of {acc}/{r} differ from the admission record")
        sp = E.load_store(str(npz[r]["store"]))
        out[r] = {"A_spec": a, "index": np.asarray(sp.index)[a], "mz": npz[r]["A_mz_mono"], "z": npz[r]["A_z"],
                  "p_seq": np.asarray(c["p_seq"], dtype=object)[ii], "p_mask": np.asarray(c["p_mask"])[ii],
                  "hyperscore_p": np.asarray(c["hyperscore_p"], np.float64)[ii], "sp": sp}
    return out


def family_screen_gpu(gi, l1, sp, a_local, mz_s, z_s):  # pragma: no cover
    """Max exhaustive hyperscore over (T∪D) ∩ W' for virtual copies of A spectra a_local with precursors (mz_s, z_s),
    isotope errors 0..2 (label-free)."""
    from dataclasses import replace

    from vulcan.gpu import engine as E
    if len(a_local) == 0:
        return np.zeros(0)
    v = subset_spectra(sp, a_local, pmz=mz_s, charge=z_s)
    p = replace(E.SearchParams(), isotope_errors=(0, 2))
    tbl, _, _ = E.search_exhaustive(v, gi, l1, p, report=1)
    s = tbl.column("spec_idx").to_numpy()
    hs = tbl.column("hyperscore").to_numpy().astype(np.float64)
    mx = np.full(len(v), -np.inf)
    np.maximum.at(mx, s, hs)
    return mx


def asyn_for_run(cert: dict, db, gi, l1, offs, seed: int, step: int = 4):  # pragma: no cover
    """Batched precomputation of the screens, then vulcan.gates.a3_asyn.construct on lookups (single specification)."""
    keys = np.array([f"{cert['run']}|{int(i)}" for i in cert["index"]], dtype=object)
    mz_a, z_a = np.asarray(cert["mz"], np.float64), np.asarray(cert["z"], np.int64)
    idf = AS.identity_screen_fn(db, cert["p_seq"], cert["p_mask"])
    look = {}
    pending = [(i, r) for i in range(len(keys)) for r in AS.RANKS]
    t0 = 0
    while pending and t0 < AS.MAX_TRIES:
        t1 = min(AS.MAX_TRIES, t0 + step)
        q = []
        for i, r in pending:
            for t in range(t0, t1):
                row = AS.draw_row(str(keys[i]), r, t, len(offs), seed)
                mz_s, z_s = float(mz_a[i]) + float(offs.delta_mz[row]), int(offs.charge[row])
                if not AS.envelope_excluded(float(mz_a[i]), int(z_a[i]), mz_s, z_s):
                    q.append((i, r, t, mz_s, z_s))
        eb = {(i, mz_s, z_s): idf(i, mz_s, z_s) for i, r, t, mz_s, z_s in q}
        need = [(i, mz_s, z_s) for (i, mz_s, z_s), b in eb.items() if not b]
        mx = family_screen_gpu(gi, l1, cert["sp"], cert["A_spec"][[x[0] for x in need]],
                               np.array([x[1] for x in need]), np.array([x[2] for x in need]))
        for (i, mz_s, z_s), m in zip(need, mx):
            look[(i, mz_s, z_s)] = (False, bool(m >= cert["hyperscore_p"][i]))
        for k, b in eb.items():
            if b:
                look[k] = (True, False)
        nxt = []
        for i, r in pending:
            ok = False
            for t in range(t0, t1):
                row = AS.draw_row(str(keys[i]), r, t, len(offs), seed)
                mz_s, z_s = float(mz_a[i]) + float(offs.delta_mz[row]), int(offs.charge[row])
                k = (i, mz_s, z_s)
                if k in look and not any(look[k]):
                    ok = True
                    break
            if not ok:
                nxt.append((i, r))
        pending, t0 = nxt, t1

    def screen(i, mz_s, z_s):
        return look[(i, mz_s, z_s)]
    tab = AS.construct(keys, mz_a, z_a, offs, screen, seed=seed)
    tab["key"] = keys
    return tab


def phase_asyn(cfg):  # pragma: no cover
    from vulcan.gates import a2_config
    from vulcan.gates.a2_identity import GateDB
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    import m3_common as C
    cfg2 = a2_config.load()
    offs = AS.load_offsets(os.path.join(ROOT, cfg["asyn"]["offsets_table"]), cfg["asyn"]["offsets_sha256"])
    idx = C.index_dir("human_contam", "prod")
    l1 = L1Index(idx)
    db = GateDB.from_l1(idx)
    gi = E.GpuIndex(l1, E.SearchParams(), STATIC, VAR)
    info, src = admitted_sources()
    d = os.path.join(OUT, "asyn")
    os.makedirs(d, exist_ok=True)
    rec = {}
    for acc in sorted({s["accession"] for s in src.values()}):
        cert = certified_A(acc, cfg2)
        for run, c in cert.items():
            f = os.path.join(d, run + ".asyn.npz")
            if not os.path.exists(f):
                c["run"] = run
                tab = asyn_for_run(c, db, gi, l1, offs, cfg["asyn"]["seed"])
                np.savez(f, **{k: v for k, v in tab.items() if k not in ("counts",)})
                json.dump({"counts": tab["counts"], "table_sha256": AS.table_hash(tab)}, open(f + ".json", "w"))
            rec[run] = json.load(open(f + ".json"))
    json.dump(rec, open(os.path.join(OUT, "asyn_record.json"), "w"), indent=1)
    return rec


def mixture_groups(mix, prov, srcrun, src, hosts):  # pragma: no cover
    """Precursor groups of one mixture: host and B spectra from their source run's hypothesis table, offsets from
    their host spectrum's shifted by the copy's offset, A / FS none.  Feature ids offset per source (as
    scripts/e7/gates_v2.mixture_groups)."""
    from vulcan.chimeric import search as CS
    from vulcan.rescore.features_v2 import read_ms1_features
    d = os.path.join(OUT, "ms1")
    n = len(mix)
    runs = {k: s["run"] for k, s in src.items()}
    host_run = mix.run.split("_h")[-1]
    hrun = hosts[int(host_run)]
    names = sorted(set(runs.values()) | {hrun})
    foff, F = {}, []
    tot = 0
    for r in names:
        mf = read_ms1_features(os.path.join(d, r + ".features.arrow"))
        foff[r] = tot
        tot += len(mf["trace_cos"])
        F.append(mf)
    mfc = {k: np.concatenate([f[k] for f in F]) for k in F[0]}
    prim = {"feature_id": np.full(n, -1, np.int64), "isolation_fraction": np.zeros(n, np.float32),
            "feature_intensity": np.zeros(n, np.float32), "iso_cos": np.zeros(n, np.float32),
            "rt_apex": np.full(n, np.nan, np.float32), "nfw": np.zeros(n, np.uint16)}
    secs = []
    from vulcan.gpu import engine as E
    host_sp = E.load_store(hosts["store"][int(host_run)])

    def take(run, sel, delta_mass=None):
        g = CS.precursor_groups(CS.read_hypotheses(os.path.join(d, run + ".hypotheses.arrow")),
                                _View(np.asarray(mix.index)[sel]))
        for k in prim:
            prim[k][sel] = g["prim"][k]
        prim["feature_id"][sel] = np.where(g["prim"]["feature_id"] >= 0, g["prim"]["feature_id"] + foff[run], -1)
        s = dict(g["sec"])
        if len(s["spec"]):
            s["feature_id"] = s["feature_id"] + foff[run]
            if delta_mass is not None:
                s["mz"] = s["mz"] + delta_mass[s["spec"]] / s["charge"].astype(np.float64)
            s["spec"] = sel[s["spec"]]
            secs.append(s)
    take(hrun, np.flatnonzero(prov == P_HOST))
    off = np.flatnonzero(prov == P_OFF)
    if off.size:
        hpos = {int(i): j for j, i in enumerate(np.asarray(host_sp.index))}
        hz = np.asarray(mix.charge)[off]
        hmz = np.asarray(host_sp.precursor_mz)[[hpos[int(i)] for i in np.asarray(mix.index)[off]]]
        dm = offset_delta_mass(np.asarray(mix.precursor_mz)[off], hmz, hz)
        take(hrun, off, delta_mass=dm)
    for k, run in runs.items():
        sel = np.flatnonzero((prov == P_B) & (srcrun == k))
        if sel.size:
            take(run, sel)
    if secs:
        sec = {k: np.concatenate([s[k] for s in secs]) for k in secs[0]}
        o = np.lexsort((sec["group"], sec["spec"]))
        sec = {k: v[o] for k, v in sec.items()}
    else:
        sec = {"spec": np.zeros(0, np.int64)}
    ng = np.ones(n, np.uint16)
    np.add.at(ng, sec["spec"], 1)
    return {"prim": prim, "sec": sec, "n_groups": ng}, mfc


class _View:
    """Minimal Spectra stand-in for precursor_groups (needs .index and len)."""

    def __init__(self, index):
        self.index = index

    def __len__(self):
        return len(self.index)


def asyn_virtual(mix, prov, srcrun, src):  # pragma: no cover
    """Virtual spectra of every A spectrum's A-syn hypotheses in this mixture; returns (Spectra, a-row table)."""
    d = os.path.join(OUT, "asyn")
    a = np.flatnonzero(prov == P_A)
    keys = a_spec_keys({k: s["run"] for k, s in src.items()}, srcrun[a], np.asarray(mix.index)[a])
    tabs = {}
    rows = {"spec": [], "rank": [], "mz": [], "charge": []}
    for j, key in zip(a, keys):
        run = key.split("|")[0]
        if run not in tabs:
            t = np.load(os.path.join(d, run + ".asyn.npz"), allow_pickle=True)
            kk = t["key"]
            tabs[run] = {"pos": {}, "t": t}
            for r_, ai in enumerate(t["a_idx"]):
                tabs[run]["pos"].setdefault(str(kk[ai]), []).append(r_)
        T = tabs[run]
        for r_ in T["pos"].get(str(key), []):
            rows["spec"].append(int(j))
            rows["rank"].append(int(T["t"]["rank"][r_]))
            rows["mz"].append(float(T["t"]["mz"][r_]))
            rows["charge"].append(int(T["t"]["charge"][r_]))
    rows = {k: np.asarray(v) for k, v in rows.items()}
    v = subset_spectra(mix, rows["spec"].astype(np.int64), pmz=rows["mz"], charge=rows["charge"]) if rows["spec"].size \
        else None
    return v, rows


def phase_search(cfg):  # pragma: no cover
    import pyarrow as pa
    import pyarrow.parquet as pq
    from dataclasses import replace
    import m3_common as C
    from vulcan.chimeric import search as CS
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    p = E.SearchParams(report_psms=5)
    l1 = L1Index(C.index_dir("human_contam", "prod"))
    gi = E.GpuIndex(l1, p, STATIC, VAR)
    info, src = admitted_sources()
    hosts = {"store": [C.E2.store_path(r) for r in C.DEV["D0"]["raw"]]}
    hosts.update({h: os.path.basename(C.E2.store_path(r))[:-len(".vspec.arrow")] for h, r in enumerate(C.DEV["D0"]["raw"])})
    d = os.path.join(OUT, "psms")
    os.makedirs(d, exist_ok=True)
    first = True
    for f in sorted(glob.glob(os.path.join(MIX, "*.spectra.pkl"))):
        name = os.path.basename(f)[:-len(".spectra.pkl")]
        out = os.path.join(d, name + ".chim.psms.parquet")
        if os.path.exists(out):
            first = False
            continue
        mix = pickle.load(open(f, "rb"))
        prov = np.load(os.path.join(MIX, name + ".prov.npy"))
        srcrun = np.load(os.path.join(MIX, name + ".srcrun.npy"))
        prim = pq.read_table(os.path.join(PSM2, name + ".psms.parquet"))
        if first:
            from vulcan.gates.a2_fs import a2_hypotheses
            from vulcan.gpu.host import pre_hypotheses
            single = np.load(os.path.join(MIX, name + ".single.npy"))
            valid, _ = pre_hypotheses(mix, gi, p)
            ref, _, _ = E.search(mix, gi, l1, p, pre_hyp=(valid, a2_hypotheses(mix, valid, gi.mass32_np, p, single)))
            if not ref.equals(prim):
                raise SystemExit(f"step4b: committed step-4 table of {name} is not reproduced by a fresh search")
            first = False
        groups, _ = mixture_groups(mix, prov, srcrun, src, hosts)
        tbl, st = CS.search_chimeric(mix, gi, l1, p, groups, primary=prim)
        pq.write_table(tbl, out)
        v, rows = asyn_virtual(mix, prov, srcrun, src)
        if v is not None:
            p2 = replace(p, isotope_errors=(0, 1))
            t2, _, _ = E.search(v, gi, l1, p2)
            vi = t2.column("spec_idx").to_numpy()
            t2 = t2.set_column(t2.column_names.index("spec_idx"), "spec_idx", pa.array(rows["spec"][vi].astype(np.int64)))
            t2 = t2.append_column("group", pa.array((rows["rank"][vi] - 1).astype(np.uint16)))
            t2 = t2.append_column("chim_rank", pa.array(rows["rank"][vi].astype(np.uint16)))
            pq.write_table(t2, os.path.join(d, name + ".asyn.psms.parquet"))
        json.dump(st, open(out + ".json", "w"), default=float)


def phase_score(cfg):  # pragma: no cover
    import pyarrow as pa
    import pyarrow.parquet as pq
    import m3_common as C
    from vulcan.index.build import L1Index
    from vulcan.rescore.features_v2 import model_psm
    v2 = models(cfg)
    fb = feature_builder(os.path.join(A3W, "predcache_dn2"))
    l1 = L1Index(C.index_dir("human_contam", "prod"))
    info, src = admitted_sources()
    hosts = {"store": [C.E2.store_path(r) for r in C.DEV["D0"]["raw"]]}
    hosts.update({h: os.path.basename(C.E2.store_path(r))[:-len(".vspec.arrow")] for h, r in enumerate(C.DEV["D0"]["raw"])})
    d = os.path.join(OUT, "psms")
    for f in sorted(glob.glob(os.path.join(d, "*.chim.psms.parquet"))):
        name = os.path.basename(f)[:-len(".chim.psms.parquet")]
        out = os.path.join(d, name + ".chim.scored.parquet")
        if os.path.exists(out):
            continue
        mix = pickle.load(open(os.path.join(MIX, name + ".spectra.pkl"), "rb"))
        prov = np.load(os.path.join(MIX, name + ".prov.npy"))
        srcrun = np.load(os.path.join(MIX, name + ".srcrun.npy"))
        _, mf = mixture_groups(mix, prov, srcrun, src, hosts)
        psm = _psm_dict(pq.read_table(f), chim_defaults=False)
        feats, finfo = fb.run_v2(l1, psm, mix, mf)
        mp = model_psm(psm)
        rt_of = fb.cache_factory(finfo["nce"]).lookup_rt     # rt_pred: input of the frozen sequence classifier (G1)
        cols = {"spec_idx": psm["spec_idx"], "group": psm["group"], "chim_rank": psm["chim_rank"], "rank": psm["rank"],
                "pf_row": psm["pf_row"], "pf_id": psm["pf_id"], "charge": psm["charge"],
                "hyperscore": psm["hyperscore"].astype(np.float64), "v2": v2.apply(mp, feats).astype(np.float64),
                "v2_noRT": v2.apply(mp, neutral_rt(feats)).astype(np.float64),
                "rt_pred": np.asarray(rt_of(l1, psm["pf_row"]), np.float32)}
        pq.write_table(pa.table(cols), out)
        fa = os.path.join(d, name + ".asyn.psms.parquet")
        if os.path.exists(fa):
            ta = _psm_dict(pq.read_table(fa))
            na = len(ta["pf_row"])
            ta["n_groups"] = np.ones(na, np.uint16)
            both = {k: np.concatenate([psm[k], np.asarray(ta[k]).astype(np.asarray(psm[k]).dtype)]) for k in psm}
            fb2, finfo2 = fb.run_v2(l1, both, mix, mf, nce=finfo["nce"])
            n0 = len(psm["pf_row"])
            if not all(np.array_equal(fb2[k][:n0].view(np.uint8), np.asarray(feats[k]).view(np.uint8)) for k in feats):
                raise SystemExit(f"step4b: A-syn rows changed the real rows' features in {name}")
            fa_ = {k: v[n0:] for k, v in fb2.items()}
            mpa = model_psm({k: v[n0:] for k, v in both.items()})
            pq.write_table(pa.table({"spec_idx": ta["spec_idx"], "group": ta["group"], "chim_rank": ta["chim_rank"],
                                     "rank": ta["rank"], "pf_row": ta["pf_row"], "pf_id": ta["pf_id"],
                                     "charge": ta["charge"], "hyperscore": ta["hyperscore"].astype(np.float64),
                                     "v2": v2.apply(mpa, fa_).astype(np.float64),
                                     "v2_noRT": v2.apply(mpa, neutral_rt(fa_)).astype(np.float64),
                                     "rt_pred": np.asarray(rt_of(l1, ta["pf_row"]), np.float32)}),
                           os.path.join(d, name + ".asyn.scored.parquet"))


def phase_v12(cfg):  # pragma: no cover
    import pyarrow as pa
    import pyarrow.parquet as pq
    import m3_common as C
    from vulcan.encoders.score import load_model
    from vulcan.index.build import L1Index
    v2 = models(cfg)
    ce = load_model(os.path.join(ROOT, cfg["scorers"]["ce_s"]["model"]))
    l1 = L1Index(C.index_dir("human_contam", "prod"))
    d = os.path.join(OUT, "psms")
    for f in sorted(glob.glob(os.path.join(MIX, "*.spectra.pkl"))):
        name = os.path.basename(f)[:-len(".spectra.pkl")]
        out = os.path.join(d, name + ".v12.parquet")
        if os.path.exists(out):
            continue
        mix = pickle.load(open(f, "rb"))
        prov = np.load(os.path.join(MIX, name + ".prov.npy"))
        single = np.load(os.path.join(MIX, name + ".single.npy"))
        a = np.flatnonzero(prov == P_A)
        r = exhaustive_v12_scores(mix, l1, single, a, True, os.path.join(A3W, "predcache_dn2"), v2_model=v2,
                                  ce_model=ce)
        pq.write_table(pa.table({"spec_idx": r["spec_idx"], "pf_row": r["pf_row"], "hyperscore": r["hyperscore"],
                                 "v2": r["v2"], "ce_s": r["ce_s"], "in_top5": r["in_top5"]}), out)


def phase_record(cfg, out_path):  # pragma: no cover
    import pyarrow.parquet as pq
    d = os.path.join(OUT, "psms")
    rec = {"mixtures": {}, "ms1": json.load(open(os.path.join(OUT, "ms1_record.json"))),
           "asyn": json.load(open(os.path.join(OUT, "asyn_record.json")))}
    for f in sorted(glob.glob(os.path.join(MIX, "*.spectra.pkl"))):
        name = os.path.basename(f)[:-len(".spectra.pkl")]
        prov = np.load(os.path.join(MIX, name + ".prov.npy"))
        r = {}
        for suf in (".chim.psms.parquet", ".chim.scored.parquet", ".asyn.psms.parquet", ".asyn.scored.parquet",
                    ".v12.parquet"):
            p_ = os.path.join(d, name + suf)
            r[suf.strip(".")] = sha256_file(p_) if os.path.exists(p_) else None
        t = pq.read_table(os.path.join(d, name + ".chim.psms.parquet"), columns=["spec_idx", "group", "chim_rank"])
        r["unit_counts"] = unit_counts(prov, t.column("chim_rank").to_numpy(), t.column("spec_idx").to_numpy(),
                                       t.column("group").to_numpy())
        fa = os.path.join(d, name + ".asyn.psms.parquet")
        if os.path.exists(fa):
            ta = pq.read_table(fa, columns=["spec_idx", "group", "chim_rank"])
            key = ta.column("spec_idx").to_numpy().astype(np.int64) * 64 + ta.column("group").to_numpy()
            u, first = np.unique(key, return_index=True)
            rk = ta.column("chim_rank").to_numpy()[first]
            r["asyn_units"] = {str(k): int((rk == k).sum()) for k in (2, 3)}
        rec["mixtures"][name] = r
    sw = os.path.join(OUT, "swap_test.json")
    rec["swap_test"] = json.load(open(sw)) if os.path.exists(sw) else None
    json.dump(rec, open(out_path, "w"), indent=1, default=float)
    return rec


def phase_swap_test(cfg):  # pragma: no cover
    import pyarrow.parquet as pq
    import tempfile
    import m3_common as C
    from vulcan.chimeric import search as CS
    from vulcan.gates import a2_config
    from vulcan.gates.a2_identity import GateDB
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index, build
    gate = C.index_dir("human_contam", "prod")
    meta = json.load(open(os.path.join(gate, "meta.json")))
    dd = tempfile.mkdtemp(prefix="a3_swap_", dir=A3W)
    fasta = meta["fasta"] if isinstance(meta["fasta"], list) else [meta["fasta"]]
    build(fasta, dd, params=meta["params"], swap_labels=True)
    same = {f: sha256_file(os.path.join(gate, f)) == sha256_file(os.path.join(dd, f))
            for f in ("peptidoforms.arrow", "residues.npy")}
    info, src = admitted_sources()
    name = sorted(os.path.basename(f)[:-len(".spectra.pkl")] for f in glob.glob(os.path.join(MIX, "*.spectra.pkl")))[0]
    mix = pickle.load(open(os.path.join(MIX, name + ".spectra.pkl"), "rb"))
    prov = np.load(os.path.join(MIX, name + ".prov.npy"))
    srcrun = np.load(os.path.join(MIX, name + ".srcrun.npy"))
    hosts = {"store": [C.E2.store_path(r) for r in C.DEV["D0"]["raw"]]}
    hosts.update({h: os.path.basename(C.E2.store_path(r))[:-len(".vspec.arrow")] for h, r in enumerate(C.DEV["D0"]["raw"])})
    p = E.SearchParams(report_psms=5)
    l1 = L1Index(dd)
    gi = E.GpuIndex(l1, p, STATIC, VAR)
    from vulcan.gates.a2_fs import a2_hypotheses
    from vulcan.gpu.host import pre_hypotheses
    single = np.load(os.path.join(MIX, name + ".single.npy"))
    valid, _ = pre_hypotheses(mix, gi, p)
    prim, _, _ = E.search(mix, gi, l1, p, pre_hyp=(valid, a2_hypotheses(mix, valid, gi.mass32_np, p, single)))
    groups, _ = mixture_groups(mix, prov, srcrun, src, hosts)
    tbl, _ = CS.search_chimeric(mix, gi, l1, p, groups, primary=prim)
    ref = pq.read_table(os.path.join(OUT, "psms", name + ".chim.psms.parquet"))
    cols = {c: bool(np.array_equal(np.asarray(ref.column(c).to_numpy(zero_copy_only=False)).view(np.uint8)
                                   if ref.column(c).type != "string" else ref.column(c).to_numpy(zero_copy_only=False),
                                   np.asarray(tbl.column(c).to_numpy(zero_copy_only=False)).view(np.uint8)
                                   if tbl.column(c).type != "string" else tbl.column(c).to_numpy(zero_copy_only=False)))
            for c in ref.column_names}
    # A-syn construction of the first mixture's A spectra on the swapped index
    cfg2 = a2_config.load()
    offs = AS.load_offsets(os.path.join(ROOT, cfg["asyn"]["offsets_table"]), cfg["asyn"]["offsets_sha256"])
    db = GateDB.from_l1(dd)
    gi0 = E.GpuIndex(l1, E.SearchParams(), STATIC, VAR)
    a_runs = sorted({src[int(k)]["run"] for k in srcrun[prov == P_A]})
    accs = sorted({s["accession"] for s in src.values() if s["run"] in a_runs})
    asyn_same = {}
    for acc in accs:
        cert = certified_A(acc, cfg2)
        for run, c in cert.items():
            if run not in a_runs:
                continue
            c["run"] = run
            tab = asyn_for_run(c, db, gi0, l1, offs, cfg["asyn"]["seed"])
            ref_h = json.load(open(os.path.join(OUT, "asyn", run + ".asyn.npz.json")))["table_sha256"]
            asyn_same[run] = AS.table_hash(tab) == ref_h
    res = {"mixture": name, "label_free_files_identical": same, "chimeric_columns_identical": cols,
           "asyn_tables_identical": asyn_same,
           "pass": bool(all(same.values()) and all(cols.values()) and all(asyn_same.values()))}
    json.dump(res, open(os.path.join(OUT, "swap_test.json"), "w"), indent=1)
    return res


def main():  # pragma: no cover
    from vulcan.gates import a3_config
    ap = argparse.ArgumentParser()
    for k in ("ms1", "asyn", "search", "score", "v12", "swap-test"):
        ap.add_argument("--" + k, action="store_true")
    ap.add_argument("--record", metavar="OUT_JSON")
    a = ap.parse_args()
    cfg = a3_config.load()
    require_freeze(cfg)
    os.makedirs(OUT, exist_ok=True)
    if a.ms1:
        phase_ms1()
    if a.asyn:
        phase_asyn(cfg)
    if a.search:
        phase_search(cfg)
    if a.score:
        phase_score(cfg)
    if a.v12:
        phase_v12(cfg)
    if a.swap_test:
        print("P1 swap test", "PASS" if phase_swap_test(cfg)["pass"] else "FAIL")
    if a.record:
        phase_record(cfg, a.record)


if __name__ == "__main__":
    main()
