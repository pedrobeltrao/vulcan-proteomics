#!/usr/bin/env python3
"""C9 pilot (D28.2): label-free yield of human 'removed-peptide' nulls on human DEV runs (D0, D8h).

Per run (one process, resident GPU index; DEV data only; never a sealed set):
  1. Certification search (NOT a null search): production closed search of the run against the FULL human+contaminant
     prod index (E2 SearchParams, report_psms 5), hyperscore (fixed, G0-verified, not under test), paired double
     competition per run (hash tie-break seed 11).  Certified = rank-1 PSM is a HUMAN target with PSM-level q <= 0.001
     and peptidoform q <= 0.01, peptide not A1.3-near the contaminant FASTA (A2.2a analogue; vulcan.gates.a2_admission.
     certify_run).  The labels of the full index are read here ONLY for this certification, the run's own hyperscore
     1% threshold (host threshold, as scripts/a2/host_threshold.py) and the identified set I.
  2. Identified set I (per PROJECT): I/L-collapsed stripped target peptides accepted at peptidoform q <= 0.01 by the
     same double competition (rank-1 PSMs only; see NOTES: secondary hits on the null spectrum itself are NOT used,
     because selecting them on the null spectrum would be a label-asymmetric (X-type, A2.1) removal).
  3. Null search: every certified spectrum, single precursor hypothesis (certified mono m/z, z, isotope 0; A2.2(e)),
     EXHAUSTIVE K2 over every T∪D peptidoform in the +-10 ppm window (vulcan.gpu.engine.search_exhaustive).
  4. Removal (label-free, pair-closed): a candidate row c is removed iff c OR its pair partner
       (i)  is I/L-identical to, Hamming <= 2 from any window of, or shares a 6-mer with any peptide in I or any
            isobaric variant of it (N<->GG, Q<->AG) (vulcan.gates.homolog.HomologFilter over I ∪ variants), or
       (ii) shares >= 20% of the certified peptide p's b2..b(L-1)/y2..y(L-1) ions (z=1, +-20 ppm, cross-type) with
            its b/y ladder (A2.2(b) f(c), exact numpy), all peptidoforms (mods) included.
     The remaining candidates are the reduced database T∪D \\ R for that spectrum.  No target/decoy label of any null
     candidate is read: the index labels are never joined to the exhaustive table.
  5. Label-free counts: certified; certified with a single MS1 feature in the isolation window (E11/E7 MS1 finder
     output, if present); A_cand = certified spectra with >= 1 surviving candidate (>= 4 matched peaks);
     tail proxy = A_cand spectra whose label-collapsed max surviving hyperscore >= the host 1% hyperscore threshold
     (own run, and the A2 convention max over the two D0 hosts); per-stage wall and process-CPU seconds.

    conda_envs/pe/bin/python scripts/c9/pilot_yield.py [--runs D0 D8h] [--out results/2026-10-28_C9_pilot/per_run]
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
T1_WORK = os.environ.get("VULCAN_T1_WORK", f"/tmp/vulcan_t1_{os.getuid()}")  # bench/t1_driver.py default work root
import resource
import sys
import time
from dataclasses import replace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))

IDX = "/tmp/vulcan_e2/index/human_contam_prod"
CONTAM_FASTA = os.path.join(ROOT, "data", "fasta", "contaminants_frankenfield2022.fasta")
D0_STORES = sorted(glob.glob("/tmp/vulcan_e2/store/b19*_293T_proteinID_*.vspec.arrow"))
D0_MS1 = T1_WORK + "/D0/ms1"
D8H_STORES = sorted(glob.glob(T1_WORK + "/D8h/stage/E9_*.vspec.arrow"))
D8H_MS1 = T1_WORK + "/D8h/ms1"
OUT = os.path.join(ROOT, "results", "2026-10-28_C9_pilot", "per_run")
TB_SEED = 11
Q_PSM, Q_PEP, Q_ID = 0.001, 0.01, 0.01
F_THR, FRAG_PPM = 0.20, 20.0
F_DIAG = (0.15, 0.10, 0.05)          # label-free sensitivity of the ladder screen (diagnostic only)
PROTON, H2O = 1.0072764, 18.010565


def cpu_s():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


class Clock:
    def __init__(self):
        self.t = {}

    def __call__(self, name):
        clk = self

        class _C:
            def __enter__(s):
                s.w, s.c = time.perf_counter(), cpu_s()

            def __exit__(s, *a):
                d = clk.t.setdefault(name, [0.0, 0.0])
                d[0] += time.perf_counter() - s.w
                d[1] += cpu_s() - s.c
        return _C()


# ---------------------------------------------------------------------------------------------- full-index labels
def full_labels(idx):
    """Full prod index labels: CERTIFICATION + host threshold + identified set only (never joined to null rows)."""
    import pyarrow as pa
    with pa.OSFile(os.path.join(idx, "labels.arrow"), "rb") as src:
        lab = pa.ipc.open_file(src).read_all()
    with pa.OSFile(os.path.join(idx, "proteins.arrow"), "rb") as src:
        prot = pa.ipc.open_file(src).read_all()
    acc = open(os.path.join(idx, "accessions.txt")).read().split("\n")
    pidx = np.load(os.path.join(idx, "protein_idx.npy"))
    src = {}
    for pep, o, n in zip(prot.column("peptide").to_numpy(zero_copy_only=False), prot.column("prot_off").to_numpy(),
                         prot.column("prot_len").to_numpy()):
        src[pep] = {("CON|" if "|Cont_" in acc[int(i)] else "OWN|") + acc[int(i)]
                    for i in pidx[int(o):int(o) + int(n)]}
    return {"is_decoy": lab.column("is_decoy").to_numpy(zero_copy_only=False),
            "tdc_pair": lab.column("tdc_pair").to_numpy(zero_copy_only=False)}, src


# ---------------------------------------------------------------------------------------------- ladders (vectorised)
_MONO = None


def _mono():
    global _MONO
    if _MONO is None:
        from vulcan.gates.a2_identity import SAGE_MONO
        st = SAGE_MONO.copy()
        st[ord("C") - 65] += 57.021464
        vm = np.zeros(26)
        vm[ord("M") - 65] = 15.9949
        _MONO = (st, vm)
    return _MONO


def ladders(l1, rows):
    """For each row: z=1 b_k (k=1..L-1) and y_k (k=1..L-1) m/z with the row's static/variable mods.
    Returns (ions, ion_index k, offsets) flat over rows (2(L-1) ions per row: b first then y)."""
    st, vm = _mono()
    rows = np.asarray(rows, np.int64)
    L = l1.length[rows].astype(np.int64)
    n = rows.size
    starts = np.cumsum(L) - L
    pos = np.arange(L.sum()) - np.repeat(starts, L)
    seg = np.repeat(np.arange(n), L)
    res = l1.residues[np.repeat(l1.seq_off[rows].astype(np.int64), L) + pos].astype(np.int64) - 65
    bit = (l1.mask[rows].astype(np.uint64)[seg] >> pos.astype(np.uint64)) & np.uint64(1)
    m = st[res] + vm[res] * bit.astype(np.float64)
    cs = np.cumsum(m)
    pre = cs - np.repeat(cs[starts] - m[starts], L)            # inclusive prefix sum within the row
    tot = np.repeat(pre[starts + L - 1], L)
    keep = pos <= np.repeat(L - 2, L)                         # k = pos + 1 in 1..L-1
    b = pre[keep] + PROTON
    y = tot[keep] - pre[keep] + H2O + PROTON                  # y index = L-1-pos
    kb = pos[keep] + 1
    ky = (np.repeat(L, L) - 1 - pos)[keep]
    segk = seg[keep]
    nk = L - 1
    off = np.zeros(n + 1, np.int64)
    np.cumsum(2 * nk, out=off[1:])
    ions = np.empty(off[-1])
    kidx = np.empty(off[-1], np.int64)
    # place b then y per row
    within = np.arange(segk.size) - np.repeat(np.cumsum(nk) - nk, nk)
    ions[off[segk] + within] = b
    kidx[off[segk] + within] = kb
    ions[off[segk] + nk[segk] + within] = y
    kidx[off[segk] + nk[segk] + within] = ky
    return ions, kidx, off


def f_scores(p_ions, cand_ions, cand_off, ppm=FRAG_PPM):
    """f(c) for every candidate: fraction of p ions with any c ion within +-ppm (of the p ion)."""
    nc = cand_off.size - 1
    if p_ions.size == 0 or nc == 0:
        return np.zeros(nc)
    cid = np.repeat(np.arange(nc), np.diff(cand_off))
    o = np.argsort(cand_ions, kind="stable")
    s, sc = cand_ions[o], cid[o]
    tol = p_ions * ppm * 1e-6
    lo = np.searchsorted(s, p_ions - tol, "left")
    hi = np.searchsorted(s, p_ions + tol, "right")
    cnt = hi - lo
    if cnt.sum() == 0:
        return np.zeros(nc)
    j = np.repeat(np.arange(p_ions.size), cnt)
    k = np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt) + np.repeat(lo, cnt)
    pairs = np.unique(j.astype(np.int64) * nc + sc[k])
    return np.bincount(pairs % nc, minlength=nc) / p_ions.size


# ---------------------------------------------------------------------------------------------- MS1 co-isolation
def single_feature_flags(store_path, ms1_dir, run):
    """bool per MS2 spectrum (load_store order): n_features_in_window <= 1 (label-free MS1).  None if no MS1 table."""
    import pyarrow as pa
    from vulcan.io import store
    hp = os.path.join(ms1_dir, run + ".hypotheses.arrow")
    if not os.path.exists(hp):
        return None
    t = store.read_table(store_path)
    lvl = t.column("ms_level").to_numpy()
    sid = t.column("spectrum_id").to_numpy()[lvl == 2]
    h = pa.ipc.open_file(pa.memory_map(hp)).read_all()
    hs = h.column("spectrum_id").to_numpy()
    nf = h.column("n_features_in_window").to_numpy().astype(np.int64)
    u, inv = np.unique(hs, return_inverse=True)
    mx = np.zeros(u.size, np.int64)
    np.maximum.at(mx, inv, nf)
    pos = np.searchsorted(u, sid)
    pos = np.minimum(pos, u.size - 1)
    found = u[pos] == sid
    return np.where(found, mx[pos] <= 1, False), float(found.mean())


# ---------------------------------------------------------------------------------------------- per run
def certify(run, sp, gi, l1, p, lab, src, near_contam, clk):
    """Stage 1+2 (full DB; labels read here only)."""
    from vulcan.gates import a2_admission as AD
    from vulcan.gates import stats as GS
    from vulcan.gpu import engine as E
    from vulcan.gpu.host import NEUTRON
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import double_competition
    with clk("search_full"):
        tbl, st, _ = E.search(sp, gi, l1, p)
    with clk("certify"):
        spec = tbl.column("spec_idx").to_numpy().astype(np.int64)
        rows = tbl.column("pf_row").to_numpy().astype(np.int64)
        pf = tbl.column("pf_id").to_numpy()
        hs = tbl.column("hyperscore").to_numpy().astype(float)
        dec = lab["is_decoy"][rows]
        pair = lab["tdc_pair"][rows]
        tb = hash_tiebreak(pf, TB_SEED)
        useq = {int(r): l1.sequence(int(r)) for r in np.unique(rows)}
        il = np.array([useq[int(r)].replace("I", "L") for r in rows], dtype=object)
        res = AD.certify_run(spec, pf, pair, hs, dec, il, src, near_contam, tb, Q_PSM, Q_PEP)
        w = res["rows"]
        dc = double_competition(spec, pf, pair, hs, dec, tiebreak=tb)
        thr = GS.tail_threshold(dc.score, dc.is_decoy, dc.q, 0.01)
        ok = (dc.q <= Q_ID) & ~dc.is_decoy
        ident = {useq[int(rows[i])].replace("I", "L") for i in dc.row[ok]}
        z = tbl.column("charge").to_numpy()[w].astype(np.int64)
        iso = tbl.column("isotope_error").to_numpy()[w].astype(np.int64)
        s_idx = spec[w]
        cert = {"spec": s_idx, "row": rows[w], "z": z,
                "mz_mono": sp.precursor_mz[s_idx] - iso * float(NEUTRON) / z, "hs_p": hs[w]}
    return cert, ident, thr, {"n_ms2": len(sp), "n_searched": int(st.counts["n_searched"]),
                              "n_rank1": res["n_rank1"], "n_certified": res["n_certified"],
                              "n_identified_peptides_run": len(ident), "search_full_gpu_wall_s": st.wall_s["total"]}


def null_stage(sp, gi, l1, p, cert, hf_I, clk):
    """Stages 3-4.  LABEL-FREE: receives no label array."""
    from vulcan.gates.a2_fs import certified_precursor
    from vulcan.gpu import engine as E
    with clk("search_null_exhaustive"):
        sub = certified_precursor(sp, cert["spec"], cert["mz_mono"], cert["z"])
        p0 = replace(p, isotope_errors=(0, 0))
        try:
            tbl, _, _ = E.search_exhaustive(sub, gi, l1, p0, report=1 << 30)
            s = tbl.column("spec_idx").to_numpy().astype(np.int64)    # index into sub == index into cert arrays
            r = tbl.column("pf_row").to_numpy().astype(np.int64)
            h = tbl.column("hyperscore").to_numpy().astype(float)
        except TypeError:                                             # engine: no (hypothesis, row) pair at all
            s = r = np.zeros(0, np.int64)
            h = np.zeros(0)
    with clk("removal_homolog"):
        par = l1.partner[r].astype(np.int64)
        par = np.where(par < 0, r, par)
        allrows = np.unique(np.concatenate([r, par]))
        seqs = {int(x): l1.sequence(int(x)) for x in allrows}
        useq = sorted(set(seqs.values()))
        nr = dict(zip(useq, hf_I.near(useq)))
        flag_row = {x: nr[seqs[x]] for x in seqs}
        rem_h = np.array([flag_row[int(a)] or flag_row[int(b)] for a, b in zip(r, par)], bool)
    with clk("removal_ladder"):
        n = len(cert["spec"])
        o = np.argsort(s, kind="stable")
        bnd = np.searchsorted(s[o], np.arange(n + 1))
        pr = cert["row"]
        p_ions_all, p_k, p_off = ladders(l1, pr)
        Lp = l1.length[pr].astype(np.int64)
        fpair = np.zeros(r.size)
        need = np.flatnonzero(~rem_h)
        ur = np.unique(np.concatenate([r[need], par[need]])) if need.size else np.zeros(0, np.int64)
        c_ions, _, c_off = ladders(l1, ur) if ur.size else (np.zeros(0), None, np.zeros(1, np.int64))
        for i in range(n):
            sel = o[bnd[i]:bnd[i + 1]]
            sel = sel[~rem_h[sel]]
            if sel.size == 0:
                continue
            pi = p_ions_all[p_off[i]:p_off[i + 1]]
            pk = p_k[p_off[i]:p_off[i + 1]]
            pi = pi[(pk >= 2) & (pk <= Lp[i] - 1)]                   # b2..b(L-1), y2..y(L-1)
            cr = np.unique(np.concatenate([r[sel], par[sel]]))
            ci = np.searchsorted(ur, cr)
            cnt = c_off[ci + 1] - c_off[ci]
            loc_off = np.zeros(cr.size + 1, np.int64)
            np.cumsum(cnt, out=loc_off[1:])
            take = np.repeat(c_off[ci], cnt) + (np.arange(cnt.sum()) - np.repeat(loc_off[:-1], cnt))
            f = f_scores(pi, c_ions[take], loc_off)
            fm = dict(zip(cr.tolist(), f.tolist()))
            fpair[sel] = [max(fm[int(a)], fm[int(b)]) for a, b in zip(r[sel], par[sel])]
    rem_f = fpair >= F_THR
    keep = ~(rem_h | rem_f)
    mx = np.full(len(cert["spec"]), -np.inf)
    np.maximum.at(mx, s[keep], h[keep])
    # label-free diagnostics: stricter ladder thresholds, and the ladder overlap f of each spectrum's winning row
    mx_f = {}
    for t in F_DIAG:
        k2 = ~rem_h & (fpair < t)
        m = np.full(len(cert["spec"]), -np.inf)
        np.maximum.at(m, s[k2], h[k2])
        mx_f[t] = m
    win_f = np.full(len(cert["spec"]), np.nan)
    kk = np.flatnonzero(keep)
    if kk.size:
        o2 = kk[np.lexsort((-h[kk], s[kk]))]
        first = np.r_[True, s[o2][1:] != s[o2][:-1]]
        win_f[s[o2[first]]] = fpair[o2[first]]
    return mx, mx_f, win_f, {"n_candidate_rows": int(r.size), "removed_rows_homolog": int(rem_h.sum()),
                             "removed_rows_ladder": int((rem_f & ~rem_h).sum()), "kept_rows": int(keep.sum()),
                             "spectra_with_any_candidate_before_removal": int(np.unique(s).size)}


def run_project(name, runs, gi, l1, p, lab, src, near_contam, d0_thr, ms1_dir, out_dir):
    """runs: list of (run_name, store_path).  Certification of all runs first (project-level I), then nulls."""
    from vulcan.gates.homolog import HomologFilter
    from vulcan.gates.a2_identity import isobaric_variants
    from vulcan.gpu import engine as E
    per = {}
    ident_proj = set()
    for run, path in runs:
        clk = Clock()
        with clk("load_store"):
            sp = E.load_store(path)
        cert, ident, thr, info = certify(run, sp, gi, l1, p, lab, src, near_contam, clk)
        ident_proj |= ident
        per[run] = (sp, cert, thr, info, clk, path)
    if d0_thr is None and name == "PXD001468":      # A2 convention: max over the two D0 hosts
        d0_thr = max(v[2] for v in per.values())
    clkp = Clock()
    with clkp("build_I_filter"):
        extra = []
        for q in ident_proj:
            v = isobaric_variants(q, 64)
            if v:
                extra += [x for x in v if x != q]
        hf_I = HomologFilter(sorted(ident_proj) + extra)
    res = {}
    for run, (sp, cert, thr, info, clk, path) in per.items():
        mx, mx_f, win_f, rinfo = null_stage(sp, gi, l1, p, cert, hf_I, clk)
        sf = single_feature_flags(path, ms1_dir, run) if ms1_dir else None
        cand = np.isfinite(mx)
        single = sf[0][cert["spec"]] if sf is not None else None
        rec = {"project": name, "run": run, **info, **rinfo,
               "n_identified_peptides_project": len(ident_proj),
               "host_thr_own_run": thr, "host_thr_D0max": d0_thr,
               "A_cand": int(cand.sum()),
               "frac_certified_with_candidate": float(cand.mean()) if cand.size else 0.0,
               "tail_proxy_own_thr": int((mx >= thr).sum()),
               "tail_proxy_D0max_thr": int((mx >= d0_thr).sum()) if d0_thr is not None else None,
               "certified_single_ms1_feature": int(single.sum()) if single is not None else None,
               "A_cand_single_ms1_feature": int((cand & single).sum()) if single is not None else None,
               "tail_proxy_own_thr_single_ms1_feature": int(((mx >= thr) & single).sum()) if single is not None
               else None,
               "ms1_mapped_fraction": sf[1] if sf is not None else None,
               "timing_wall_cpu_s": {k: [round(v[0], 2), round(v[1], 2)] for k, v in clk.t.items()},
               "project_I_filter_wall_cpu_s": {k: [round(v[0], 2), round(v[1], 2)] for k, v in clkp.t.items()}}
        tail = mx >= thr
        rec["diag_ladder_threshold"] = {str(t): {"A_cand": int(np.isfinite(m).sum()), "tail_proxy_own_thr":
                                                 int((m >= thr).sum())} for t, m in mx_f.items()}
        bins = [0, 0.05, 0.10, 0.15, 0.20]
        rec["diag_winner_f_hist_bins"] = bins
        rec["diag_winner_f_hist_A_cand"] = np.histogram(win_f[cand], bins + [1.01])[0].tolist()[:4]
        rec["diag_winner_f_hist_tail"] = np.histogram(win_f[tail], bins + [1.01])[0].tolist()[:4]
        if single is not None:
            rec["diag_coisolation"] = {"single_A_cand": int((cand & single).sum()),
                                       "single_tail": int((tail & single).sum()),
                                       "multi_A_cand": int((cand & ~single).sum()),
                                       "multi_tail": int((tail & ~single).sum())}
        # max-hyperscore quantiles of A_cand (label-collapsed) for the proxy-vs-threshold sensitivity
        if cand.any():
            rec["max_hs_A_cand_quantiles"] = {str(q): float(np.quantile(mx[cand], q))
                                              for q in (0.5, 0.9, 0.99, 0.999)}
        json.dump(rec, open(os.path.join(out_dir, run + ".json"), "w"), indent=1, default=float)
        print(json.dumps({k: rec[k] for k in ("run", "n_ms2", "n_certified", "certified_single_ms1_feature",
                                              "A_cand", "tail_proxy_own_thr", "tail_proxy_D0max_thr",
                                              "host_thr_own_run")}), flush=True)
        res[run] = rec
    return res, d0_thr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sets", nargs="+", default=["D0", "D8h"])
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--resume", action="store_true", help="skip projects whose run JSONs exist; D0 threshold from them")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    import m3_common as C
    from vulcan import hw
    from vulcan.gates.homolog import HomologFilter, read_fasta_seqs
    from vulcan.gpu import engine as E
    from vulcan.index.build import L1Index
    t0 = time.perf_counter()
    meta = {"hw": hw.fingerprint(), "loadavg_start": os.getloadavg(), "index": IDX}
    clk = Clock()
    with clk("setup_index"):
        l1 = L1Index(IDX)
        p = E.SearchParams(report_psms=C.REPORT)
        gi = E.GpuIndex(l1, p, C.STATIC, C.VAR)
        lab, src = full_labels(IDX)
        near_contam = HomologFilter(read_fasta_seqs([CONTAM_FASTA])).near
    meta["setup_wall_cpu_s"] = {k: [round(v[0], 2), round(v[1], 2)] for k, v in clk.t.items()}
    e9 = {r["id"]: r["accession"] for r in csv.DictReader(open(os.path.join(ROOT, "data/manifests/e9_d8.tsv")),
                                                           delimiter="\t")}
    projects = []
    if "D0" in a.sets:
        projects.append(("PXD001468", [(os.path.basename(x)[:-len(".vspec.arrow")], x) for x in D0_STORES], D0_MS1))
    if "D8h" in a.sets:
        for x in D8H_STORES:
            run = os.path.basename(x)[:-len(".vspec.arrow")]
            projects.append((e9[run.split("__")[0]], [(run, x)], D8H_MS1))
    if a.limit:
        projects = projects[:a.limit]
    # D0 host thresholds first (A2 convention: max over the two D0 hosts) -- from D0's own certification search
    d0_thr = None
    allres = {}
    if a.resume:
        done = {os.path.basename(f)[:-5] for f in glob.glob(os.path.join(a.out, "*.json"))}
        d0 = [json.load(open(os.path.join(a.out, r + ".json"))) for r, _ in projects[0][1] if r in done] \
            if projects and projects[0][0] == "PXD001468" else []
        d0_thr = max(x["host_thr_own_run"] for x in d0) if d0 else None
        projects = [pr for pr in projects if not all(r in done for r, _ in pr[1])]
        meta["resumed"] = True
    for name, runs, ms1 in projects:
        r, d0_thr = run_project(name, runs, gi, l1, p, lab, src, near_contam, d0_thr, ms1, a.out)
        allres.update(r)
    meta["loadavg_end"] = os.getloadavg()
    meta["total_wall_s"] = time.perf_counter() - t0
    meta["d0_host_thr_max"] = d0_thr
    mp = os.path.join(a.out, "_meta.json")
    if a.resume and os.path.exists(mp):
        meta = {"first": json.load(open(mp)).get("first", json.load(open(mp))), "resume": meta}
    else:
        meta = {"first": meta}
    json.dump(meta, open(mp, "w"), indent=1, default=float)
    print("done", round(time.perf_counter() - t0), "s")


if __name__ == "__main__":
    main()
