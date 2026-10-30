"""L0.5 pilot pass (DESIGN §4): seeded subsample -> wide-tolerance closed search -> calibration and run diagnostics.

    res = run_pilot(store_path, fasta, workdir, searcher=SageSearcher(annotate=True))
    res.calibration.to_json(...)      # offsets + empirical tolerances (vulcan.calib.calibration)
    res.report                        # JSON-able diagnostics

Outputs (per run):
  * precursor and fragment ppm-offset curves vs RT and m/z (Calibration; offsets, never baked into an index);
  * empirical tolerances before and after recalibration;
  * analyzer/activation from metadata, cross-checked against the data (fragment error scale, ion-series intensity,
    low-mass cutoff);
  * reporter-ion detection (TMT/TMTpro/iTRAQ) on MS2 and MS3, against off-reporter control positions;
  * a label-blind mass-offset histogram from an open search of the same subsample (rank-1 PSMs, no labels used);
  * charge and isotope-error priors (overall, and per hypothesis source when a hypothesis table is given).

Label use: the pilot's IDs are never reported. Anchors are all rank-1 PSMs whose q-value (TDC on the searcher's own
decoys) is <= q_anchor, decoys included: the selection is a score threshold applied identically to both labels.
"""
from __future__ import annotations

import os
import re
import time
from dataclasses import asdict, dataclass, field

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pcsv

from ..io import schema as S
from ..io import store
from .calibration import Calibration, OffsetModel, core_tolerance, fragment_tolerance, precursor_tolerance
from .searcher import SageSearcher, SearchConfig, Searcher, SpectraBatch, load_batch

PROTON = 1.007276466812
NEUTRON = 1.0033548378
H2O = 18.0105646837
NH3 = 17.0265491
AA = {"G": 57.02146372, "A": 71.03711379, "S": 87.03203535, "P": 97.05276385, "V": 99.06841392, "T": 101.0476785,
      "C": 103.0091845, "L": 113.0840640, "I": 113.0840640, "N": 114.0429275, "D": 115.0269431, "Q": 128.0585775,
      "K": 128.0949630, "E": 129.0425931, "M": 131.0404846, "H": 137.0589119, "F": 147.0684139, "R": 156.1011110,
      "Y": 163.0633285, "W": 186.0793130, "U": 150.9536, "O": 237.1477}
REPORTERS = {
    "TMT10/11": [126.127726, 127.124761, 127.131081, 128.128116, 128.134436, 129.131471, 129.137790, 130.134825,
                 130.141145, 131.138180, 131.144499],
    "TMTpro_extra": [132.141535, 132.147855, 133.144890, 133.151210, 134.148245, 134.154565, 135.151600],
    "iTRAQ4": [114.1112, 115.1083, 116.1116, 117.1150],
    "iTRAQ8_extra": [113.1078, 118.1112, 119.1153, 121.1220],
}
# offsets (Da) annotated in the mass-offset histogram (Unimod monoisotopic deltas)
KNOWN_OFFSETS = {
    "isotope +1": 1.00335, "deamidation": 0.98402, "oxidation": 15.99491, "dioxidation": 31.98983,
    "carbamidomethyl": 57.02146, "acetyl": 42.01057, "carbamyl": 43.00581, "formyl": 27.99491,
    "methyl": 14.01565, "dimethyl": 28.03130, "trimethyl": 42.04695, "phospho": 79.96633, "GlyGly": 114.04293,
    "TMT6/10/11": 229.16293, "TMTpro": 304.20715, "iTRAQ4": 144.10206, "iTRAQ8": 304.20536,
    "SILAC K+8": 8.01420, "SILAC R+10": 10.00827, "SILAC K+4": 4.02511, "SILAC R+6": 6.02013,
    "Na adduct": 21.98194, "K adduct": 37.95588, "ammonia loss": -17.02655, "water loss": -18.01056,
    "Met loss": -131.04049, "Met loss+acetyl": -89.02992, "amidation": -0.98402, "pyro-carbamidomethyl": 39.99492,
    "Cys->Dha": -33.98772, "sulfo": 79.95682, "HexNAc": 203.07937, "Hex": 162.05282, "propionamide": 71.03711,
    "+Gly (missed/extra residue)": 57.02146, "+Ala": 71.03711, "+Ser": 87.03203,
}


@dataclass
class PilotParams:
    fraction: float = 0.10                 # seeded subsample of MS2 (DESIGN: 5-10%)
    seed: int = 11
    q_anchor: float = 0.01
    wide_precursor_ppm: tuple = (-50.0, 50.0)
    wide_fragment_ppm: tuple = (-50.0, 50.0)   # high-resolution MS2
    wide_fragment_da: tuple = (-0.6, 0.6)      # ion-trap MS2
    isotope_errors: tuple = (-1, 3)
    precursor_charge: tuple = (2, 4)
    open_search: bool = True
    open_da: tuple = (-150.0, 500.0)       # open-search error window: delta = exp - calc in [-150, +500] Da
    open_bin_da: float = 0.01
    reporter_ppm: float = 10.0
    max_annotate: int = 3000               # anchors used for the ion-series check
    tol_k: float = 4.0                     # tolerance = mu +- tol_k * sd of the sigma-clipped error core
    floor_precursor_ppm: float = 1.0
    floor_fragment: float = 1.0            # ppm (Da runs: 0.01)


@dataclass
class PilotResult:
    calibration: Calibration
    report: dict
    psms: pa.Table
    anchors: pa.Table
    timings: dict = field(default_factory=dict)


# ------------------------------------------------------------------------------------------------ peptides
_MOD = re.compile(r"\[([+-][\d.]+)\]")


def parse_peptide(s: str):
    """Sage notation -> (residues, per-residue mass incl. mods, n-term mod mass)."""
    nterm = 0.0
    m = re.match(r"^\[([+-][\d.]+)\]-?", s)
    if m:
        nterm = float(m.group(1))
        s = s[m.end():]
    res, mass = [], []
    i = 0
    while i < len(s):
        c = s[i]
        if c == "[":
            j = s.index("]", i)
            mass[-1] += float(s[i + 1:j])
            i = j + 1
            continue
        if c == "-":
            i += 1
            continue
        res.append(c)
        mass.append(AA.get(c, 0.0))
        i += 1
    return "".join(res), np.asarray(mass), nterm


def fragment_ladders(peptides, kinds=("b", "y", "c", "z")):
    """Flat singly-charged fragment m/z per peptide: dict kind -> (values, offsets)."""
    out = {k: ([], [0]) for k in kinds}
    for pep in peptides:
        _, m, nt = parse_peptide(pep)
        pre = np.cumsum(m)[:-1] + nt
        tot = m.sum() + nt
        ions = {"b": pre + PROTON, "y": (tot - pre) + H2O + PROTON, "c": pre + NH3 + PROTON,
                "z": (tot - pre) + H2O - NH3 + PROTON + 1.007825}   # z-dot (z+1), the dominant ETD z ion
        for k in kinds:
            v, o = out[k]
            v.append(ions[k][::-1] if k in ("y", "z") else ions[k])
            o.append(o[-1] + len(ions[k]))
    return {k: (np.concatenate(v) if v else np.empty(0), np.asarray(o, np.int64)) for k, (v, o) in out.items()}


def match_ladders(batch: SpectraBatch, rows: np.ndarray, ions: np.ndarray, ion_off: np.ndarray, tol_ppm: float,
                  tol_da: float | None = None):
    """Nearest peak for every theoretical ion of PSM k (spectrum batch row rows[k]). Returns (ion idx, peak idx, err)."""
    n_k = len(rows)
    owner = np.repeat(np.arange(n_k), np.diff(ion_off))
    spec = rows[owner]
    W = 1e5
    pk_spec = np.repeat(np.arange(len(batch)), np.diff(batch.offsets))
    key = pk_spec * W + np.asarray(batch.mz, np.float64)
    q = spec * W + ions
    pos = np.searchsorted(key, q)
    best = np.full(len(q), -1, np.int64)
    bd = np.full(len(q), np.inf)
    for c in (pos - 1, pos):
        cc = np.clip(c, 0, len(key) - 1)
        ok = (c >= 0) & (c < len(key)) & (pk_spec[cc] == spec)
        d = np.where(ok, np.abs(key[cc] - q), np.inf)
        b = d < bd
        best, bd = np.where(b, cc, best), np.where(b, d, bd)
    lim = tol_da if tol_da is not None else ions * tol_ppm * 1e-6
    hit = bd <= lim
    err = (np.asarray(batch.mz, np.float64)[np.clip(best, 0, None)] - ions)
    return np.flatnonzero(hit), best[hit], err[hit]


# ------------------------------------------------------------------------------------------------ pieces
def subsample(n: int, fraction: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    k = max(1, int(round(n * fraction)))
    return np.sort(rng.choice(n, size=k, replace=False))


def anchor_select(psms: pa.Table, q: float) -> pa.Table:
    """Label-symmetric anchors: rank-1 PSMs with q <= q (targets and decoys alike)."""
    return psms.filter(pc.less_equal(psms.column("q_value"), q))


def precursor_errors(anchors: pa.Table) -> np.ndarray:
    exp = anchors.column("exp_mass").to_numpy()
    calc = anchors.column("calc_mass").to_numpy()
    iso = anchors.column("isotope_error").to_numpy().astype(float)
    return (exp - iso * NEUTRON - calc) / calc * 1e6


def read_fragments(path: str) -> pa.Table:
    return pcsv.read_csv(path, parse_options=pcsv.ParseOptions(delimiter="\t"))


def analyzer_activation_meta(path: str) -> dict:
    t = store.read_table(path, columns=["ms_level", "analyzer", "activation", "supplemental_activation", "nce",
                                        "filter_string"])
    out = {}
    for lv in (2, 3):
        u = t.filter(pc.equal(t.column("ms_level"), lv))
        if u.num_rows == 0:
            continue
        an = u.column("analyzer").to_numpy()
        ac = np.asarray(store.column_numpy(u, "activation"), float)
        sa = np.asarray(store.column_numpy(u, "supplemental_activation"), float)
        nce = np.asarray(store.column_numpy(u, "nce"), float)
        cnt = lambda a, names: {names.get(int(k), str(k)): int(v) for k, v in zip(*np.unique(a[~np.isnan(a)].astype(int)
                                                                                          if a.dtype.kind == "f" else a,
                                                                                          return_counts=True))}
        out[f"ms{lv}"] = {"n": int(u.num_rows), "analyzer": cnt(an, S.ANALYZER), "activation": cnt(ac, S.ACTIVATION),
                          "supplemental_activation": cnt(sa, S.ACTIVATION),
                          "nce_median": float(np.nanmedian(nce)) if np.isfinite(nce).any() else None}
    return out


def reporter_check(batches: dict, ppm: float, frag_calib=None) -> dict:
    """Fraction of spectra with a peak within ppm of each reporter vs at control positions (reporter + d)."""
    ctrl_d = np.array([-0.045, -0.035, 0.035, 0.045, 0.055])
    out = {}
    for lv, b in batches.items():
        if b is None or len(b) == 0:
            continue
        mz = np.asarray(b.mz, np.float64)
        spec = np.repeat(np.arange(len(b)), np.diff(b.offsets))
        sel = mz < 140
        mz, spec = mz[sel], spec[sel]
        key = spec * 1e3 + mz

        def frac(t):
            lo = np.searchsorted(key, spec_ids * 1e3 + t * (1 - ppm * 1e-6))
            hi = np.searchsorted(key, spec_ids * 1e3 + t * (1 + ppm * 1e-6), side="right")
            return float(np.mean(hi > lo))
        spec_ids = np.arange(len(b))
        res = {}
        for plex, ions in REPORTERS.items():
            f = [frac(t) for t in ions]
            c = [frac(t + d) for t in ions for d in ctrl_d]
            res[plex] = {"median_frac": float(np.median(f)), "control_frac": float(np.mean(c)),
                         "per_ion": [round(x, 3) for x in f]}
        tmt = res["TMT10/11"]
        itr = res["iTRAQ4"]
        res["detected"] = {
            "TMT": bool(tmt["median_frac"] >= 0.3 and tmt["median_frac"] >= 5 * max(tmt["control_frac"], 1e-3)),
            "TMTpro": bool(res["TMTpro_extra"]["median_frac"] >= 0.3 and res["TMTpro_extra"]["median_frac"]
                           >= 5 * max(res["TMTpro_extra"]["control_frac"], 1e-3)),
            "iTRAQ": bool(itr["median_frac"] >= 0.3 and itr["median_frac"] >= 5 * max(itr["control_frac"], 1e-3)),
        }
        out[f"ms{lv}"] = res
    return out


def ion_series_check(batch: SpectraBatch, anchors: pa.Table, row_of: dict, tol_ppm: float, tol_da=None,
                     max_n: int = 3000, prec_mz=None) -> dict:
    """Matched-intensity share of b/y vs c/z ions and the low-mass-cutoff signature on anchor spectra."""
    a = anchors.filter(pc.invert(anchors.column("is_decoy"))) if anchors.num_rows else anchors
    a = a.slice(0, max_n)
    if a.num_rows < 20:
        return {"n": int(a.num_rows)}
    sid = a.column("spectrum_id").to_numpy()
    rows = np.array([row_of[int(s)] for s in sid.tolist()], np.int64)
    lad = fragment_ladders(a.column("peptide").to_pylist())
    tot_int = np.add.reduceat(np.asarray(batch.intensity, np.float64), batch.offsets[:-1])[rows] \
        if len(batch.intensity) else np.zeros(len(rows))
    share = {}
    for k, (ions, off) in lad.items():
        idx, pk, _ = match_ladders(batch, rows, ions, off, tol_ppm, tol_da)
        it = np.asarray(batch.intensity, np.float64)[np.unique(pk)].sum() if len(pk) else 0.0
        share[k] = float(it / max(tot_int.sum(), 1e-30))
    by, cz = share["b"] + share["y"], share["c"] + share["z"]
    # low-mass cutoff: fraction of peaks below 0.28 * precursor m/z (ion-trap CID 1/3 rule leaves ~none)
    lens = np.diff(batch.offsets)[rows]
    pk_rows = np.repeat(rows, lens)
    pm = np.repeat(prec_mz[rows], lens)
    src = np.concatenate([np.arange(batch.offsets[r], batch.offsets[r + 1]) for r in rows])
    low = float(np.mean(np.asarray(batch.mz)[src] < 0.28 * pm)) if len(src) else float("nan")
    return {"n": int(a.num_rows), "intensity_share": share, "by_share": by, "cz_share": cz,
            "frac_peaks_below_0.28_precursor_mz": low, "inferred_activation":
                ("ETD-like" if cz > 1.5 * by else "collisional (b/y)"),
            "low_mass_cutoff_present": bool(low < 0.01)}


def offset_histogram(psms: pa.Table, bin_da: float, lo=-150.0, hi=500.0, top: int = 25) -> dict:
    """Label-blind histogram of delta = exp_mass - calc_mass over all rank-1 open-search PSMs."""
    d = (psms.column("exp_mass").to_numpy() - psms.column("calc_mass").to_numpy())
    d = d[(d >= lo) & (d <= hi)]
    edges = np.arange(lo, hi + bin_da, bin_da)
    h, _ = np.histogram(d, edges)
    cen = 0.5 * (edges[:-1] + edges[1:])
    # peaks: local maxima of a 3-bin sum over a robust background (running median over 1 Da)
    s3 = np.convolve(h, np.ones(3), "same")
    w = int(round(1.0 / bin_da))
    from scipy.ndimage import median_filter
    bg = median_filter(s3.astype(float), size=2 * w + 1, mode="nearest") + 1.0
    z = (s3 - bg) / np.sqrt(bg)
    is_max = (s3 >= np.roll(s3, 1)) & (s3 > np.roll(s3, -1))
    cand = np.flatnonzero(is_max & (z >= 5) & (s3 >= 5))
    cand = cand[np.argsort(-s3[cand])][:top]
    peaks = []
    names = np.array(list(KNOWN_OFFSETS))
    vals = np.array(list(KNOWN_OFFSETS.values()))
    for c in cand:
        m = d[np.abs(d - cen[c]) <= 1.5 * bin_da]
        mc = float(np.median(m)) if len(m) else float(cen[c])
        j = int(np.argmin(np.abs(vals - mc)))
        ann = str(names[j]) if abs(vals[j] - mc) <= 0.01 else ("unmodified" if abs(mc) <= 0.01 else None)
        if ann is None and abs(mc - round(mc / NEUTRON) * NEUTRON) <= 0.01 and abs(round(mc / NEUTRON)) <= 3:
            ann = f"isotope {int(round(mc / NEUTRON)):+d}"
        peaks.append({"delta_da": round(mc, 4), "count_3bin": int(s3[c]), "z": round(float(z[c]), 1),
                      "annotation": ann})
    return {"n_psms": int(len(d)), "bin_da": bin_da, "peaks": peaks,
            "frac_unshifted_|d|<0.02": float(np.mean(np.abs(d) < 0.02)) if len(d) else None}


def labels_from_offsets(hist: dict) -> dict:
    ann = {p["annotation"]: p["count_3bin"] for p in hist.get("peaks", []) if p["annotation"]}
    return {"TMT_label_offset": "TMT6/10/11" in ann or "TMTpro" in ann, "iTRAQ_label_offset": any("iTRAQ" in a for a in ann),
            "SILAC_offset": any(a.startswith("SILAC") for a in ann), "dimethyl_offset": "dimethyl" in ann}


# ------------------------------------------------------------------------------------------------ driver
def run_pilot(path: str, fasta: str, workdir: str, searcher: Searcher | None = None, p: PilotParams | None = None,
              hypotheses: pa.Table | None = None) -> PilotResult:
    p = p or PilotParams()
    searcher = searcher or SageSearcher(annotate=True)
    os.makedirs(workdir, exist_ok=True)
    tm = {}
    t0 = time.perf_counter()
    meta = store.run_meta(path)
    full = load_batch(path, 2)
    meta_aa = analyzer_activation_meta(path)
    an = meta_aa.get("ms2", {}).get("analyzer", {})
    lowres = an.get("IT", 0) > an.get("OT", 0) + an.get("TOF", 0) + an.get("Astral", 0)
    sel = subsample(len(full), p.fraction, p.seed)
    sub = full.take(sel)
    tm["load_subsample"] = time.perf_counter() - t0
    cfg = SearchConfig(fasta=fasta, precursor_ppm=p.wide_precursor_ppm,
                       fragment_ppm=None if lowres else p.wide_fragment_ppm,
                       fragment_da=p.wide_fragment_da if lowres else None, isotope_errors=p.isotope_errors,
                       precursor_charge=p.precursor_charge)
    psms, acct = searcher.search(sub, cfg, os.path.join(workdir, "wide"))
    tm["wide_search"] = acct
    psms = psms.filter(pc.is_finite(psms.column("calc_mass")))
    anchors = anchor_select(psms, p.q_anchor)
    row_of = {int(s): i for i, s in enumerate(sub.spectrum_id.tolist())}
    a_rows = np.array([row_of[int(s)] for s in anchors.column("spectrum_id").to_numpy().tolist()], np.int64)
    rt_a = sub.rt[a_rows]
    mz_a = sub.precursor_mz[a_rows]
    e_pre = precursor_errors(anchors)
    pre_model = OffsetModel.fit(rt_a, mz_a, e_pre)
    e_pre_cal = e_pre - pre_model(rt_a, mz_a)

    # fragments (Sage --annotate-matches; matched within the wide window)
    frag_unit = "Da" if lowres else "ppm"
    fr_model = OffsetModel()
    ef = ef_cal = np.empty(0)
    w = abs(p.wide_fragment_da[1]) if lowres else abs(p.wide_fragment_ppm[1])
    if "fragments" in acct:
        ft = read_fragments(acct["fragments"])
        a_pid = anchors.column("psm_id").to_numpy()
        ft = ft.filter(pa.array(np.isin(ft.column("psm_id").to_numpy(), a_pid)))
        calc = ft.column("fragment_mz_calculated").to_numpy().astype(float)
        expm = ft.column("fragment_mz_experimental").to_numpy().astype(float)
        rt_of = dict(zip(a_pid.tolist(), rt_a.tolist()))
        f_rt = np.array([rt_of[x] for x in ft.column("psm_id").to_numpy().tolist()])
        if lowres:
            ef = ef_cal = expm - calc
        else:
            ef = (expm - calc) / calc * 1e6
            fr_model = OffsetModel.fit(f_rt, calc, ef, n_min_per_bin=200)
            ef_cal = ef - fr_model(f_rt, calc)
    cal = Calibration(run_id=int(meta["run_id"]), run_name=meta["run_name"], precursor=pre_model, fragment=fr_model,
                      fragment_unit=frag_unit)
    tol_info = {}
    tols = {}
    for tag, e_b, e_a, win, unit, floor in (("precursor", e_pre, e_pre_cal, abs(p.wide_precursor_ppm[1]), "ppm",
                                             p.floor_precursor_ppm),
                                            ("fragment", ef, ef_cal, w, frag_unit,
                                             p.floor_fragment if not lowres else 0.01)):
        for when, e in (("before", e_b), ("after", e_a)):
            core, ci = core_tolerance(e, k=p.tol_k, floor=floor)
            cov, vi = fragment_tolerance(e, win, coverage=0.99) if len(e) >= 100 else (None, {"n": int(len(e))})
            tols[f"{tag}_{unit}_{when}"] = core
            tols[f"{tag}_{unit}_{when}_cov99"] = cov
            tol_info[f"{tag}_{when}"] = {"core": ci, "cov99": vi}
    tols["search_wide"] = {"precursor_ppm": list(p.wide_precursor_ppm),
                           "fragment": list(p.wide_fragment_da if lowres else p.wide_fragment_ppm), "unit": frag_unit}
    tols["rule"] = f"core: mu +- {p.tol_k} sd (sigma-clipped); cov99: background-subtracted 99% coverage"
    cal.tolerances = tols
    tol_frag_after = tols[f"fragment_{frag_unit}_after"]
    t1 = time.perf_counter()
    # priors
    zc = anchors.column("charge").to_numpy()
    ie = anchors.column("isotope_error").to_numpy()
    priors = {"charge": {str(int(k)): int(v) for k, v in zip(*np.unique(zc, return_counts=True))},
              "isotope_error_instrument_mz": {str(int(k)): int(v) for k, v in zip(*np.unique(ie, return_counts=True))}}
    if hypotheses is not None:
        priors["by_source"] = hypothesis_priors(hypotheses, anchors, cal)
    # analyzer / activation cross-check
    ion = ion_series_check(sub, anchors, row_of, tol_ppm=20.0 if not lowres else 0.0,
                           tol_da=0.5 if lowres else None, max_n=p.max_annotate, prec_mz=sub.precursor_mz)
    inferred_an = None
    if tol_frag_after:
        half = (tol_frag_after[1] - tol_frag_after[0]) / 2
        inferred_an = "high-res (OT/TOF/Astral)" if (frag_unit == "ppm" and half <= 30) else "low-res (IT)"
    meta_an = max(an, key=an.get) if an else None
    meta_ac = meta_aa.get("ms2", {}).get("activation", {})
    meta_ac = max(meta_ac, key=meta_ac.get) if meta_ac else None
    consistent = None
    if inferred_an and meta_an:
        consistent = (inferred_an.startswith("high") == (meta_an in ("OT", "TOF", "Astral")))
        if "inferred_activation" in ion and meta_ac:
            consistent = consistent and ((ion["inferred_activation"] == "ETD-like") == (meta_ac in ("ETD", "ECD")))
    aa_check = {"metadata": meta_aa, "data": {"fragment_tolerance_after": tol_frag_after, "fragment_unit": frag_unit,
                                              "inferred_analyzer": inferred_an, **ion},
                "metadata_major": {"analyzer": meta_an, "activation": meta_ac}, "consistent": consistent}
    # reporters (MS2 + MS3 of the same subsample fraction)
    batches = {2: sub}
    try:
        b3 = load_batch(path, 3)
        if len(b3):
            batches[3] = b3.take(subsample(len(b3), p.fraction, p.seed))
    except Exception:
        pass
    rep = reporter_check(batches, max(p.reporter_ppm, 0.0) if not lowres else 400.0)
    tm["diagnostics_s"] = time.perf_counter() - t1
    # open search: label-blind offset histogram
    hist = None
    if p.open_search:
        ocfg = SearchConfig(fasta=fasta, precursor_da=p.open_da, fragment_ppm=None if lowres else (-20.0, 20.0),
                            fragment_da=(-0.5, 0.5) if lowres else None, isotope_errors=(0, 0),
                            precursor_charge=p.precursor_charge, variable_mods={}, max_variable_mods=0,
                            predict_rt=False)
        opsms, oacct = SageSearcher(annotate=False).search(sub, ocfg, os.path.join(workdir, "open"))
        tm["open_search"] = oacct
        hist = offset_histogram(opsms, p.open_bin_da, lo=p.open_da[0], hi=p.open_da[1])
    report = {
        "run_name": meta["run_name"], "run_id": int(meta["run_id"]), "instrument_model": meta.get("instrument_model"),
        "params": asdict(p), "n_ms2": int(len(full)), "n_subsample": int(len(sub)),
        "n_psms": int(psms.num_rows), "n_anchors": int(anchors.num_rows),
        "n_anchor_decoys": int(pc.sum(anchors.column("is_decoy")).as_py() or 0),
        "precursor_model": pre_model.to_dict(), "fragment_model": fr_model.to_dict(),
        "tolerance_info": tol_info,
        "tolerances": cal.tolerances, "priors": priors, "analyzer_activation": aa_check, "reporters": rep,
        "offset_histogram": hist, "labels": labels_from_offsets(hist) if hist else None,
        "precursor_error_quantiles_before": np.quantile(e_pre, [0.01, 0.25, 0.5, 0.75, 0.99]).round(3).tolist()
        if len(e_pre) else None,
        "precursor_error_quantiles_after": np.quantile(e_pre_cal, [0.01, 0.25, 0.5, 0.75, 0.99]).round(3).tolist()
        if len(e_pre) else None,
    }
    cal.extra = {"priors": priors, "reporters": rep.get("ms3", rep.get("ms2", {})).get("detected"),
                 "analyzer_activation": aa_check["metadata_major"], "labels": report["labels"]}
    tm["total_s"] = time.perf_counter() - t0
    return PilotResult(cal, report, psms, anchors, tm)


def hypothesis_priors(hyp: pa.Table, anchors: pa.Table, cal: Calibration | None = None, tol_ppm: float = 20.0) -> dict:
    """For each anchor, which hypothesis rows (same spectrum) explain its peptide mass, by source and isotope error.

    Label-symmetric (anchors as selected). Reports, per source, the isotope-error histogram of the best-explaining row
    and the fraction of anchors explained by any row of the table (MS1-feature validation by identifications)."""
    sid_h = hyp.column("spectrum_id").to_numpy()
    o = np.argsort(sid_h, kind="stable")
    sid_h = sid_h[o]
    mass = hyp.column("neutral_mass").to_numpy()[o]
    src = hyp.column("source").to_numpy()[o]
    ie = hyp.column("isotope_error").to_numpy()[o]
    z = hyp.column("charge").to_numpy()[o]
    sid_a = anchors.column("spectrum_id").to_numpy()
    calc = anchors.column("calc_mass").to_numpy()
    za = anchors.column("charge").to_numpy()
    lo = np.searchsorted(sid_h, sid_a)
    hi = np.searchsorted(sid_h, sid_a, side="right")
    n = len(sid_a)
    best_src = np.full(n, -1)
    best_ie = np.full(n, -99)
    n_rows = hi - lo
    owner = np.repeat(np.arange(n), n_rows)
    idx = np.concatenate([np.arange(a, b) for a, b in zip(lo, hi)]) if n else np.empty(0, np.int64)
    ok = (np.abs(mass[idx] - calc[owner]) / calc[owner] * 1e6 <= tol_ppm) & (z[idx] == za[owner])
    owner, idx = owner[ok], idx[ok]
    # best row per anchor: source priority, then |isotope error|
    if len(owner):
        oo = np.lexsort((np.abs(ie[idx]), src[idx], owner))
        owner, idx = owner[oo], idx[oo]
        f = np.r_[True, owner[1:] != owner[:-1]]
        best_src[owner[f]] = src[idx[f]]
        best_ie[owner[f]] = ie[idx[f]]
    from ..ms1.hypotheses import SOURCE_NAME
    out = {"n_anchors": int(n), "frac_explained": float(np.mean(best_src >= 0)) if n else None,
           "frac_explained_same_charge_by_source": {SOURCE_NAME[c]: float(np.mean(best_src == c)) for c in SOURCE_NAME}}
    for c, name in SOURCE_NAME.items():
        m = best_src == c
        out[name] = {str(int(k)): int(v) for k, v in zip(*np.unique(best_ie[m], return_counts=True))}
    return out
