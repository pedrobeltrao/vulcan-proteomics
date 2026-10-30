#!/usr/bin/env python3
"""A6.4 spectrum admission for HN1 (per admitted run / project).  Step 2-4 code; LABEL-FREE on every null candidate.

Label guard (A6.4.6.3): the ONLY functions here that receive full-index labels are `load_full_labels`,
`certify_run` and `swap_recompute_I` (the certification module).  They index the label arrays only at rows of the HOST
full-index certification search (all reported PSMs of the run's own spectra, i.e. rank-1 PSMs and their competitors)
and never at a row that came from a reduced-index search or a null copy (asserted by `LabelView`, which records every
row it was indexed at; tests/a6/test_admit.py).  Everything after certification (`project_removal`, `null_stage`,
`freeze`) receives no label array.  No output of this module carries a per-label key.

Per project (A6.3.4(e), A6.4):
  1. certify_run (A6.4.1 hyperscore part): production closed search of the run on the FULL prod index; paired double
     competition per run (hash tie-break seed 11); certified = rank-1 human target, PSM q <= 0.001, peptidoform
     q <= 0.01, not A1.3-near the contaminant FASTA (vulcan.gates.a2_admission.certify_run, as scripts/c9/pilot_yield).
     Run threshold thr = lowest accepted-target score at peptidoform q <= 0.01 (vulcan.gates.stats.tail_threshold).
     Symmetric I (A6.4.3): unordered pairs (label-free key min(row, l1.partner[row])) whose double-competition pair
     winner has score >= thr, WHICHEVER label won.  Count assertion: |I| equals the number of distinct label-free pair
     keys whose max rank-1 score is >= thr, computed without any label (psm_competition + max per pair key).
  2. sage_agree (A6.4.1 second engine): Sage (tools/bin/sage, T1 plain comparator config, DB = prod T + its decoys
     prod_1001.fasta) rank-1 target with spectrum_q <= 0.001, same I/L-collapsed stripped sequence and same charge.
  3. project_removal (A6.4.4): R_P = pairs with either member (i) A1.3-near I (both member sequences of every I pair,
     plus isobaric variants), (ii) in the single-hypothesis precursor window of any certified p of the project with
     ladder overlap f >= 0.05 (z=1, +-20 ppm, cross-type, all modforms; pilot_yield.ladders/f_scores), or (iii) a
     contaminant / Panel B neighbourhood.  Evaluated over every index row that lies in the precursor window of ANY
     certified hypothesis of the project (window widened by 0.5 ppm: stricter) and their pair partners: rows outside
     those windows can never be a candidate of any null search (single hypothesis, iso 0), so every search output is
     identical to removing R_P from the whole index (reduced_index exactness test).  Pair-closed.
  4. freeze (A6.4.6.1): certified list, I, thresholds, R_P written + sha256 manifest BEFORE the reduced index is built.
  5. null_stage (A6.4.5): reduced index (reduced_index.ReducedL1), production search of every A spectrum with the
     single certified hypothesis; A_cand (>= 1 candidate with >= 4 matched peaks); tail proxy = label-collapsed max
     hyperscore >= the run's own thr.
  6. swap_test (A6.4.6.2): rebuild the index with swap_labels=True (G0 knob; peptidoforms/residues asserted
     bit-identical), keep frozen only the certified list and thr scalars, recompute double competitions, I, R_P, the
     reduced index, A, A_cand and tail-proxy memberships from the RAW full-index PSM table and the swapped index; all
     must be bit-identical.  A mismatch is an implementation stop.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import sys
import time
from dataclasses import replace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
for _p in (ROOT, os.path.join(ROOT, "scripts", "m3"), os.path.join(ROOT, "scripts", "c9"), HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

IDX = "/tmp/vulcan_e2/index/human_contam_prod"
IDX_SWAPPED = "/tmp/vulcan_a6/index_prod_swapped"
CONTAM_FASTA = os.path.join(ROOT, "data", "fasta", "contaminants_frankenfield2022.fasta")
PANELS_FASTA = os.path.join(ROOT, "data", "fasta", "a2_panels.fasta")
TB_SEED = 11
Q_PSM, Q_PEP, Q_ID = 0.001, 0.01, 0.01
F_MAX = 0.05
SAGE_Q = 1e-3
WINDOW_MARGIN_PPM = 0.5
MIN_MATCHED = 4
PROTON = 1.0072764


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------------------------------------ label view (guard)
class LabelView:
    """Full-index labels, indexable only through `at(rows, origin)`; records every row read (label guard test).
    origin must be "host_full_search"; any other origin (reduced-index search, null copy) raises."""

    def __init__(self, is_decoy, tdc_pair):
        self._d, self._p = is_decoy, tdc_pair
        self.rows_read = []

    def at(self, rows, origin):
        if origin != "host_full_search":
            raise PermissionError(f"label read refused for rows of origin {origin!r}")
        rows = np.asarray(rows, np.int64)
        self.rows_read.append(rows)
        return self._d[rows], self._p[rows]


def load_full_labels(idx=IDX):
    """Certification module only: (LabelView, peptide -> sources)."""
    from pilot_yield import full_labels
    lab, src = full_labels(idx)
    return LabelView(lab["is_decoy"], lab["tdc_pair"]), src


def pair_key(l1, rows):
    rows = np.asarray(rows, np.int64)
    par = l1.partner[rows].astype(np.int64)
    par = np.where(par < 0, rows, par)
    return np.minimum(rows, par)


# ------------------------------------------------------------------------------------------------ 1. certification
def identified_pairs(dc, pk, thr):
    """Symmetric I (A6.4.3): label-free pair keys of the double-competition pair winners with score >= thr,
    WHICHEVER label won.  dc: DoubleCompetitionResult over rank-1 rows; pk: label-free pair key per row."""
    return np.unique(np.asarray(pk)[dc.row[dc.score >= thr]])


def check_I_count(I_keys, spec, hs, pk, tb, thr):
    """A6.4.6.2 count assertion: |I| (and the set) equals the distinct label-free pair keys whose max rank-1 score is
    >= thr, computed without any label."""
    from vulcan.stats.tdc import psm_competition
    r1 = psm_competition(spec, hs, tiebreak=tb)
    uk, inv = np.unique(np.asarray(pk)[r1], return_inverse=True)
    mx = np.full(uk.size, -np.inf)
    np.maximum.at(mx, inv, np.asarray(hs, float)[r1])
    lf = np.unique(uk[mx >= thr])
    if not np.array_equal(lf, np.asarray(I_keys)):
        raise AssertionError(f"|I| count assertion failed: {len(I_keys)} (double competition) vs {lf.size} "
                             "(label-free pair maxima)")
    return int(lf.size)


def certify_run(sp, gi, l1, p, labels: LabelView, src, near_contam, raw=None):
    """Full-index certification search + certification + thr + symmetric I.  raw: optional precomputed raw PSM
    table dict (swap test re-uses the frozen raw table).  Returns (cert, I_keys, thr, raw, info)."""
    from vulcan.gates import a2_admission as AD
    from vulcan.gates import stats as GS
    from vulcan.gpu import engine as E
    from vulcan.gpu.host import NEUTRON
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import double_competition, psm_competition
    st_wall = None
    if raw is None:
        tbl, st, _ = E.search(sp, gi, l1, p)
        st_wall = st.wall_s["total"]
        raw = {c: tbl.column(c).to_numpy() for c in ("spec_idx", "pf_row", "pf_id", "hyperscore", "charge",
                                                      "isotope_error", "rank", "matched_peaks")}
    # rank-1 PSM per spectrum first (label-free): labels are then indexed ONLY at rank-1 rows of the host search.
    # double_competition / certify_run start with the same psm_competition, so restricting to its winners is exact.
    r1 = psm_competition(raw["spec_idx"].astype(np.int64), raw["hyperscore"].astype(float),
                         tiebreak=hash_tiebreak(raw["pf_id"], TB_SEED))
    spec = raw["spec_idx"].astype(np.int64)[r1]
    rows = raw["pf_row"].astype(np.int64)[r1]
    pf = raw["pf_id"][r1]
    hs = raw["hyperscore"].astype(float)[r1]
    dec, pair = labels.at(rows, "host_full_search")
    tb = hash_tiebreak(pf, TB_SEED)
    useq = {int(r): l1.sequence(int(r)) for r in np.unique(rows)}
    il = np.array([useq[int(r)].replace("I", "L") for r in rows], dtype=object)
    res = AD.certify_run(spec, pf, pair, hs, dec, il, src, near_contam, tb, Q_PSM, Q_PEP)
    w = res["rows"]
    dc = double_competition(spec, pf, pair, hs, dec, tiebreak=tb)
    thr = GS.tail_threshold(dc.score, dc.is_decoy, dc.q, Q_ID)
    pk = pair_key(l1, rows)
    I_keys = identified_pairs(dc, pk, thr)
    check_I_count(I_keys, spec, hs, pk, tb, thr)
    z = raw["charge"][r1][w].astype(np.int64)
    iso = raw["isotope_error"][r1][w].astype(np.int64)
    s_idx = spec[w]
    cert = {"spec": s_idx, "row": rows[w], "z": z,
            "mz_mono": sp.precursor_mz[s_idx] - iso * float(NEUTRON) / z, "hs_p": hs[w]}
    info = {"n_ms2": len(sp), "n_rank1": res["n_rank1"], "n_certified_hs": res["n_certified"],
            "n_I_pairs_run": int(I_keys.size), "thr": thr, "search_full_gpu_wall_s": st_wall}
    return cert, I_keys, thr, raw, info


# ------------------------------------------------------------------------------------------------ 2. Sage
def sage_rank1(path):
    """(run stem, scan) -> (I/L stripped peptide, charge) for rank-1 TARGET PSMs with spectrum_q <= SAGE_Q.
    Sage's own target/decoy column is the second engine's certification (host spectra), not a gate label."""
    import csv
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["rank"] != "1" or r["label"] != "1" or float(r["spectrum_q"]) > SAGE_Q:
                continue
            m = re.search(r"scan=(\d+)", r["scannr"])
            pep = re.sub(r"[^A-Z]", "", r["peptide"]).replace("I", "L")
            out[(r["filename"][:-len(".mzML")], int(m.group(1)))] = (pep, int(r["charge"]))
    return out


def sage_agree(cert, l1, store_path, run, s1):
    from vulcan.io import store
    t = store.read_table(store_path)
    lvl = t.column("ms_level").to_numpy()
    scans = t.column("scan").to_numpy()[lvl == 2].astype(np.int64)
    ok = np.zeros(len(cert["spec"]), bool)
    for i, (s, r, z) in enumerate(zip(cert["spec"], cert["row"], cert["z"])):
        h = s1.get((run, int(scans[int(s)])))
        ok[i] = h is not None and h[0] == l1.sequence(int(r)).replace("I", "L") and h[1] == int(z)
    return ok


def run_sage_project(mzmls, out_dir, fasta, threads):
    import bench.t1.comparators as CMP
    return CMP.run_sage("a6", fasta, mzmls, out_dir, {"prec_ppm": 10.0, "frag_ppm": 20.0}, threads, tag="sage")


# ------------------------------------------------------------------------------------------------ 3. removal set
class ContamFilter:
    """A6.4.4(iii): occurs (I/L) in any protein of the contaminant FASTA or Panel B, or A1.3-near a NON-human protein
    of either.  Exact occurrence in a human contaminant protein: length-L windows (I/L collapsed)."""

    def __init__(self):
        from vulcan.gates.homolog import HomologFilter
        prots = []
        for path, sel in ((CONTAM_FASTA, lambda h: True), (PANELS_FASTA, lambda h: h.startswith(">PNL|B|"))):
            hdr, cur = None, []
            for line in open(path):
                if line.startswith(">"):
                    if hdr is not None and sel(hdr):
                        prots.append((hdr, "".join(cur)))
                    hdr, cur = line.strip(), []
                else:
                    cur.append(line.strip())
            if hdr is not None and sel(hdr):
                prots.append((hdr, "".join(cur)))
        human = [s for h, s in prots if "OS=Homo sapiens" in h]
        nonhuman = [s for h, s in prots if "OS=Homo sapiens" not in h]
        self.n_human, self.n_nonhuman = len(human), len(nonhuman)
        self.hf = HomologFilter(nonhuman)
        self.hum = "\0".join(s.replace("I", "L") for s in human)

    def near(self, seqs):
        a = np.asarray(self.hf.near(seqs), bool)
        b = np.array([s.replace("I", "L") in self.hum for s in seqs], bool)
        return a | b


CONTAM_CACHE = "/tmp/vulcan_a6/contam_rows.npy"


def contam_rows(l1, contam, cache=None, threads=8):
    """A6.4.4(iii) row flag for EVERY full-index row (before pair closure).  Project-independent and label-free, so it
    is computed once (fixed FASTAs x index) and cached with its provenance (index content hash, FASTA sha256s)."""
    prov = {"index_content_hash": l1.meta.get("content_hash"), "contam_fasta": sha256_file(CONTAM_FASTA),
            "panels_fasta": sha256_file(PANELS_FASTA), "n_rows": len(l1)}
    cache = cache or CONTAM_CACHE.replace(".npy", f"_{prov['index_content_hash']}.npy")
    mp = cache + ".json"
    if os.path.exists(cache) and os.path.exists(mp):
        m = json.load(open(mp))
        if m.get("prov") == prov and m.get("sha256") == sha256_file(cache):
            return np.load(cache)
    seqs = [l1.sequence(i) for i in range(len(l1))]
    useq = sorted(set(seqs))
    flags = near_parallel("__contam__", useq, "contam")
    nf = dict(zip(useq, flags))
    out = np.fromiter((nf[x] for x in seqs), bool, len(seqs))
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    np.save(cache, out)
    json.dump({"prov": prov, "sha256": sha256_file(cache), "n_flagged": int(out.sum())}, open(mp, "w"), indent=1)
    return out


_POOL = None
_W = {}


def _w_near(args):
    """Worker: HomologFilter over `ref` (cached per worker by key) .near(chunk)."""
    key, ref, chunk = args
    from vulcan.gates.homolog import HomologFilter
    if key not in _W:
        _W.clear()
        _W[key] = ContamFilter() if ref == "__contam__" else HomologFilter(ref)
    return np.asarray(_W[key].near(chunk), bool)


def start_pool(n=8):
    """Process pool for the near-homology tests; started BEFORE any CUDA initialisation (fork)."""
    global _POOL
    if _POOL is None and n > 1:
        import multiprocessing as mp
        _POOL = mp.get_context("fork").Pool(n)
    return _POOL


def near_parallel(ref, seqs, key, chunk=20000, n_parts=None):
    """bool per query: near(ref) (ref = list of sequences for HomologFilter, or "__contam__" for ContamFilter),
    evaluated in the process pool on chunks.  Element-wise identical to the serial call (queries are independent)."""
    from vulcan.gates.homolog import HomologFilter
    seqs = list(seqs)
    if not seqs:
        return np.zeros(0, bool)
    if _POOL is None or len(seqs) <= chunk:
        f = ContamFilter() if ref == "__contam__" else HomologFilter(ref)
        return np.asarray(f.near(seqs), bool)
    # one part per worker: the reference list is pickled once per part, so few large parts beat many small ones
    n = n_parts or _POOL._processes
    step = -(-len(seqs) // n)
    parts = [(key, ref, seqs[i:i + step]) for i in range(0, len(seqs), step)]
    return np.concatenate(_POOL.map(_w_near, parts, chunksize=1))


def window_rows(l1, p, mz_mono, z):
    """Row ranges [lo, hi) of every single certified hypothesis (iso 0), engine semantics, window widened by
    WINDOW_MARGIN_PPM (stricter)."""
    from vulcan.gpu.host import make_hypotheses, Spectra
    n = len(z)
    sp = Spectra(run="w", scan=np.arange(n), index=np.arange(n), rt=np.zeros(n),
                 precursor_mz=np.asarray(mz_mono, np.float64), charge=np.asarray(z).astype(np.int8),
                 mz=np.zeros(0, np.float32), intensity=np.zeros(0, np.float32), offsets=np.zeros(n + 1, np.int64))
    lo_ppm, hi_ppm = p.precursor_ppm
    pw = replace(p, isotope_errors=(0, 0), precursor_ppm=(lo_ppm - WINDOW_MARGIN_PPM, hi_ppm + WINDOW_MARGIN_PPM))
    h = make_hypotheses(sp, np.ones(n, bool), l1.mass32, pw)
    o = np.argsort(h["spec"], kind="stable")
    return h["lo"][o], h["hi"][o]


def project_removal(l1, p, cert_all, I_keys, contam_flag, f_max=F_MAX, threads=8):
    """R_P (bool over full-index rows, pair-closed) from the frozen certified list (rows, mz, z), I (pair keys) and
    the fixed FASTAs.  LABEL-FREE (receives no label array)."""
    from pilot_yield import f_scores, ladders
    from reduced_index import pair_close
    from vulcan.gates.a2_identity import isobaric_variants
    from vulcan.gates.homolog import HomologFilter
    t0 = time.perf_counter()
    lo, hi = window_rows(l1, p, cert_all["mz_mono"], cert_all["z"])
    cnt = hi - lo
    W = np.unique(np.concatenate([np.arange(a, b) for a, b in zip(lo, hi)])) if cnt.sum() else np.zeros(0, np.int64)
    par = l1.partner[W].astype(np.int64)
    W2 = np.unique(np.concatenate([W, par[par >= 0]]))
    # (i) near I (both members of every I pair) + isobaric variants
    I_rows = np.unique(np.concatenate([I_keys, l1.partner[I_keys][l1.partner[I_keys] >= 0].astype(np.int64)]))
    I_seqs = sorted({l1.sequence(int(r)).replace("I", "L") for r in I_rows})
    extra = []
    for q in I_seqs:
        v = isobaric_variants(q, 64)
        if v:
            extra += [x for x in v if x != q]
    I_ref = I_seqs + extra
    I_key = hashlib.sha256("\n".join(I_ref).encode()).hexdigest()
    seqs = {int(r): l1.sequence(int(r)) for r in W2}
    useq = sorted(set(seqs.values()))
    t1 = time.perf_counter()
    nI = dict(zip(useq, near_parallel(I_ref, useq, I_key)))
    t2 = time.perf_counter()
    t3 = time.perf_counter()
    flag = np.zeros(len(l1), bool)
    fi = np.array([nI[seqs[int(r)]] for r in W2], bool)
    fc = np.asarray(contam_flag, bool)[W2]
    flag[W2[fi | fc]] = True
    flag_ic = pair_close(l1, flag)
    # (ii) ladder, project-wide: every row in the window of certified p (rows not yet removed), f(c) >= f_max
    fl = np.zeros(len(l1), bool)
    p_ions_all, p_k, p_off = ladders(l1, cert_all["row"])
    Lp = l1.length[cert_all["row"]].astype(np.int64)
    n_lad_rows = 0
    for i in range(len(lo)):
        cr = np.arange(lo[i], hi[i])
        if cr.size == 0:
            continue
        cpar = l1.partner[cr].astype(np.int64)
        both = np.unique(np.concatenate([cr, cpar[cpar >= 0]]))
        both = both[~flag_ic[both] & ~fl[both]]
        if both.size == 0:
            continue
        n_lad_rows += both.size
        pi = p_ions_all[p_off[i]:p_off[i + 1]]
        pk = p_k[p_off[i]:p_off[i + 1]]
        pi = pi[(pk >= 2) & (pk <= Lp[i] - 1)]
        c_ions, _, c_off = ladders(l1, both)
        f = f_scores(pi, c_ions, c_off)
        fl[both[f >= f_max]] = True
    t4 = time.perf_counter()
    R = pair_close(l1, flag_ic | fl)
    info = {"n_window_rows": int(W.size), "n_window_rows_with_partners": int(W2.size),
            "n_unique_window_seqs": len(useq), "n_I_pairs": int(I_keys.size), "n_I_seqs": len(I_seqs),
            "R_rows": int(R.sum()), "R_rows_near_I_or_contam": int(flag_ic.sum()),
            "R_rows_ladder_only": int((R & ~flag_ic).sum()), "ladder_rows_tested": int(n_lad_rows),
            "wall_s": {"windows+I_filter_build": t1 - t0, "near_I": t2 - t1, "contam": t3 - t2, "ladder": t4 - t3}}
    return R, info


# ------------------------------------------------------------------------------------------------ 4. freeze
def freeze(out_dir, project, per_run, I_keys_proj, R):
    """Write the frozen admission inputs and their sha256 manifest (before the reduced index is built)."""
    os.makedirs(out_dir, exist_ok=True)
    files = {}
    for run, d in per_run.items():
        c = d["cert"]
        f = os.path.join(out_dir, f"{run}.cert.npz")
        np.savez(f, spec=c["spec"], row=c["row"], z=c["z"], mz_mono=c["mz_mono"], sage_agree=d["agree"])
        files[os.path.basename(f)] = sha256_file(f)
    f = os.path.join(out_dir, "thresholds.json")
    json.dump({r: d["thr"] for r, d in per_run.items()}, open(f, "w"), indent=1, sort_keys=True)
    files["thresholds.json"] = sha256_file(f)
    f = os.path.join(out_dir, "I_pairs.npy")
    np.save(f, I_keys_proj)
    files["I_pairs.npy"] = sha256_file(f)
    f = os.path.join(out_dir, "R_rows.npy")
    np.save(f, np.flatnonzero(R).astype(np.int64))
    files["R_rows.npy"] = sha256_file(f)
    json.dump({"project": project, "files": files}, open(os.path.join(out_dir, "freeze_manifest.json"), "w"),
              indent=1, sort_keys=True)
    return files


# ------------------------------------------------------------------------------------------------ 5. null stage
def null_stage(sp, cert, agree, thr, gir, red, p):
    """A spectra (certified by both engines) searched on the reduced index with the single certified hypothesis
    (production search).  LABEL-FREE.  Returns membership arrays over the certified list."""
    from reduced_index import map_back
    from vulcan.gates.a2_fs import certified_precursor
    from vulcan.gpu import engine as E
    A = np.asarray(agree, bool)
    idx = np.flatnonzero(A)
    mx = np.full(len(cert["spec"]), -np.inf)
    if idx.size:
        sub = certified_precursor(sp, cert["spec"][idx], cert["mz_mono"][idx], cert["z"][idx])
        p0 = replace(p, isotope_errors=(0, 0))
        try:
            tbl, _, _ = E.search(sub, gir, red, p0)
            t = map_back(tbl, red)
            ok = t["matched_peaks"] >= MIN_MATCHED
            m = np.full(idx.size, -np.inf)
            np.maximum.at(m, t["spec_idx"][ok].astype(np.int64), t["hyperscore"][ok].astype(float))
            mx[idx] = m
        except TypeError:
            pass
    A_cand = A & np.isfinite(mx)
    tail = A_cand & (mx >= thr)
    return {"A": A, "A_cand": A_cand, "tail_proxy": tail, "max_hs": mx}


# ------------------------------------------------------------------------------------------------ project driver
class Ctx:
    def __init__(self, idx=IDX, pool=8):
        start_pool(pool)                       # fork the near-test workers before CUDA is initialised
        import m3_common as C
        from pilot_yield import CONTAM_FASTA as CF
        from vulcan.gates.homolog import HomologFilter, read_fasta_seqs
        from vulcan.gpu import engine as E
        from vulcan.index.build import L1Index
        self.idx = idx
        self.l1 = L1Index(idx)
        self.p = E.SearchParams(report_psms=C.REPORT)
        self.C = C
        self.gi = E.GpuIndex(self.l1, self.p, C.STATIC, C.VAR)
        self.labels, self.src = load_full_labels(idx)
        self.near_contam = HomologFilter(read_fasta_seqs([CF])).near
        self.contam = ContamFilter()
        self.contam_flag = contam_rows(self.l1, self.contam)


def admit_project(ctx, project, runs, sage_s1, out_dir, raws=None, clk=None):
    """runs: [(run, store_path)].  sage_s1: dict from sage_rank1 (project search).  raws: optional frozen raw PSM
    tables per run (swap test).  Returns (per-run membership + counts, timing)."""
    from reduced_index import ReducedL1
    from vulcan.gpu import engine as E
    T = {}

    def tick(k, t0):
        T[k] = T.get(k, 0.0) + time.perf_counter() - t0
    per = {}
    I_proj = []
    for run, path in runs:
        t0 = time.perf_counter()
        sp = E.load_store(path)
        tick("load_store", t0)
        t0 = time.perf_counter()
        cert, I_keys, thr, raw, info = certify_run(sp, ctx.gi, ctx.l1, ctx.p, ctx.labels, ctx.src, ctx.near_contam,
                                                   raw=None if raws is None else raws[run])
        tick("certify", t0)
        t0 = time.perf_counter()
        agree = sage_agree(cert, ctx.l1, path, run, sage_s1)
        tick("sage_agree", t0)
        per[run] = {"sp": sp, "cert": cert, "thr": thr, "raw": raw, "info": info, "agree": agree}
        I_proj.append(I_keys)
    I_keys_proj = np.unique(np.concatenate(I_proj)) if I_proj else np.zeros(0, np.int64)
    # removal uses the certified list of the project (certified = both engines agree, A6.4.1)
    cert_all = {k: np.concatenate([d["cert"][k][d["agree"]] for d in per.values()]) for k in ("row", "mz_mono", "z")}
    t0 = time.perf_counter()
    R, rinfo = project_removal(ctx.l1, ctx.p, cert_all, I_keys_proj, ctx.contam_flag)
    tick("removal", t0)
    t0 = time.perf_counter()
    man = freeze(out_dir, project, per, I_keys_proj, R)
    tick("freeze", t0)
    t0 = time.perf_counter()
    red = ReducedL1(ctx.l1, R)
    gir = E.GpuIndex(red, ctx.p, ctx.C.STATIC, ctx.C.VAR)
    tick("reduced_index", t0)
    out = {}
    t0 = time.perf_counter()
    for run, d in per.items():
        m = null_stage(d["sp"], d["cert"], d["agree"], d["thr"], gir, red, ctx.p)
        out[run] = {"m": m, "counts": {"n_ms2": d["info"]["n_ms2"], "n_certified_hs": d["info"]["n_certified_hs"],
                                       "A": int(m["A"].sum()), "A_cand": int(m["A_cand"].sum()),
                                       "tail_proxy": int(m["tail_proxy"].sum()),
                                       "A_cand_hs_only": None}}
    tick("null_search", t0)
    del gir
    return {"per_run": out, "I_keys": I_keys_proj, "R": R, "removal": rinfo, "freeze": man, "timing_s": T,
            "raw": {r: d["raw"] for r, d in per.items()}, "thr": {r: d["thr"] for r, d in per.items()}}


def swap_test(ctx_swapped, project, runs, sage_s1, out_dir, ref):
    """A6.4.6.2: recompute on the swapped index from the raw PSM tables; compare I, R_P, memberships bit for bit.
    Frozen and reused: the certified list and the thr scalars (A6.4.6.2) -- certification is recomputed from the raw
    table, but the certified list and thr are taken from the frozen files and compared, not replaced."""
    from reduced_index import ReducedL1
    from vulcan.gpu import engine as E
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import double_competition, psm_competition
    t0 = time.perf_counter()
    same_files = {f: sha256_file(os.path.join(IDX, f)) == sha256_file(os.path.join(ctx_swapped.idx, f))
                  for f in ("peptidoforms.arrow", "residues.npy")}
    l1 = ctx_swapped.l1
    I_sw = []
    frozen = {}
    for run, path in runs:
        z = np.load(os.path.join(out_dir, f"{run}.cert.npz"))
        frozen[run] = {k: z[k] for k in z.files}
        raw = ref["raw"][run]
        r1 = psm_competition(raw["spec_idx"].astype(np.int64), raw["hyperscore"].astype(float),
                             tiebreak=hash_tiebreak(raw["pf_id"], TB_SEED))
        spec = raw["spec_idx"].astype(np.int64)[r1]
        rows = raw["pf_row"].astype(np.int64)[r1]
        hs = raw["hyperscore"].astype(float)[r1]
        pf = raw["pf_id"][r1]
        dec, pair = ctx_swapped.labels.at(rows, "host_full_search")
        dc = double_competition(spec, pf, pair, hs, dec, tiebreak=hash_tiebreak(pf, TB_SEED))
        I_sw.append(identified_pairs(dc, pair_key(l1, rows), ref["thr"][run]))
    I_sw = np.unique(np.concatenate(I_sw)) if I_sw else np.zeros(0, np.int64)
    cert_all = {k: np.concatenate([frozen[r][k][frozen[r]["sage_agree"]] for r, _ in runs])
                for k in ("row", "mz_mono", "z")}
    R_sw, _ = project_removal(l1, ctx_swapped.p, cert_all, I_sw, ctx_swapped.contam_flag)
    red = ReducedL1(l1, R_sw)
    gir = E.GpuIndex(red, ctx_swapped.p, ctx_swapped.C.STATIC, ctx_swapped.C.VAR)
    mem_ok = {}
    for run, path in runs:
        sp = E.load_store(path)
        fz = frozen[run]
        cert = {"spec": fz["spec"], "row": fz["row"], "z": fz["z"], "mz_mono": fz["mz_mono"]}
        m = null_stage(sp, cert, fz["sage_agree"], ref["thr"][run], gir, red, ctx_swapped.p)
        r0 = ref["per_run"][run]["m"]
        mem_ok[run] = {k: bool(np.array_equal(m[k], r0[k])) for k in ("A", "A_cand", "tail_proxy")}
    del gir
    res = {"project": project, "label_free_index_files_identical": same_files,
           "I_identical": bool(np.array_equal(I_sw, ref["I_keys"])),
           "R_identical": bool(np.array_equal(R_sw, ref["R"])), "memberships_identical": mem_ok,
           "wall_s": time.perf_counter() - t0}
    res["pass"] = bool(all(same_files.values()) and res["I_identical"] and res["R_identical"]
                       and all(all(v.values()) for v in mem_ok.values()))
    return res


def build_swapped_index(out=IDX_SWAPPED, src=IDX):
    """A6.4.6.2: rebuild the gate index from the same FASTA and params with swap_labels=True (G0 knob).  The
    label-free files must be bit-identical to the prod index (checked in swap_test)."""
    from vulcan.index.build import build
    meta = json.load(open(os.path.join(src, "meta.json")))
    fasta = meta["fasta"] if isinstance(meta["fasta"], list) else [meta["fasta"]]
    m = build(fasta, out, params=meta["params"], swap_labels=True)
    same = {f: sha256_file(os.path.join(src, f)) == sha256_file(os.path.join(out, f))
            for f in ("peptidoforms.arrow", "residues.npy")}
    if not all(same.values()):
        raise AssertionError(f"swapped index label-free files differ: {same}")
    return {"out": out, "timing_s": m.get("timing_s"), "label_free_identical": same}


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-swapped", action="store_true")
    a = ap.parse_args()
    if a.build_swapped:
        print(json.dumps(build_swapped_index(), default=float))
