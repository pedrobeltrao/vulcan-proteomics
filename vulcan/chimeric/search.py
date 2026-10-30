"""E7 chimeric search (DESIGN §4 L2 "chimeric: candidates restricted to MS1-feature precursors in the isolation window",
§10 E7).  LABEL-BLIND: inputs are the spectra, the vulcan.ms1 hypothesis/feature tables and the index; labels.arrow is
never read here.

Query units (precursor groups) of one MS2 spectrum:
  group 0 (chimeric rank 1)  the closed search exactly as vulcan.gpu.engine.search (instrument precursor m/z + charge,
                             isotope errors per SearchParams) -- the primary precursor.
  group k >= 1 (rank k + 1)  every co-isolated MS1 feature group of the hypothesis table that is NOT the instrument's
                             precursor (source ms1_feature, instrument_agrees false; charge 2-6, MS1 chimeric_rank <= 5,
                             isolation fraction >= 1%: the table defaults), ordered by MS1 chimeric_rank (label-free),
                             searched at the feature mono m/z and charge with isotope errors `sec_isotope_errors`.
Each group is its own K1 query unit with its own top-K, pair closure and K2 ranking.  This is implemented by giving
every secondary group a *virtual spectrum* that carries a copy of the real spectrum's peaks and the group's precursor,
and running the unchanged engine (K0/K1/K2 kernels and psm_table) on the virtual spectra.  The primary search is the
unchanged closed call, so with chimeric mode off nothing changes (bit-identical by construction and by test).

Output: one PSM table per run = closed PSM table (group 0) + secondary PSM tables, with PSM_SCHEMA §1 columns
(`spec_idx` = the REAL spectrum's index; `rank` is the rank within the group) plus, appended:
  group (u16)       0 primary, 1.. secondary group order
  chim_rank (u16)   1 + group (label-free chimeric rank; pool >= 3 in L4)
  n_groups (u16)    precursor groups of the spectrum (1 + #secondary)
and the MS1 context of the group's precursor (null-free, 0/NaN if none):
  ms1_feature_id (i64, -1 none), ms1_isolation_fraction, ms1_feature_intensity, ms1_iso_cos, ms1_feature_rt_apex (f32),
  n_features_in_window (u16), inst_charge (i8: the instrument charge of the spectrum, 0 unknown).
For group 0 the MS1 context is the instrument-agreeing feature group (if any).
"""
from __future__ import annotations

import dataclasses
import time

import numpy as np
import pyarrow as pa
import pyarrow.ipc as ipc

from vulcan.gpu.host import SearchParams, Spectra

SEC_ISOTOPE_ERRORS = (0, 1)
MS1_COLS = ["ms1_feature_id", "ms1_isolation_fraction", "ms1_feature_intensity", "ms1_iso_cos", "ms1_feature_rt_apex",
            "n_features_in_window"]


def read_hypotheses(path: str) -> pa.Table:
    with pa.OSFile(path, "rb") as f:
        return ipc.open_file(f).read_all()


def precursor_groups(hyp: pa.Table, sp: Spectra) -> dict:
    """Per-spectrum precursor groups from a hypothesis table, mapped onto `sp` (engine MS2 order).

    Returns dict with
      prim_*  [n_spectra] MS1 context of the primary (instrument-agreeing feature group; -1/NaN if none)
      sec     dict of arrays, one row per secondary group, sorted by (spec, MS1 chimeric_rank):
              spec, group (1..), mz, charge, feature_id, isolation_fraction, feature_intensity, iso_cos, rt_apex, nfw
    """
    cols = ["spectrum_id", "precursor_group", "source", "charge", "isotope_error", "precursor_mz", "feature_id",
            "feature_intensity", "isolation_fraction", "chimeric_rank", "feature_iso_cos", "feature_rt_apex",
            "instrument_agrees", "n_features_in_window"]
    t = hyp.select(cols)
    g = {c: t.column(c).to_numpy(zero_copy_only=False) for c in cols}
    # one row per (spectrum, group): the isotope-error-0 row (every group has e = 0 in the table)
    m0 = (g["isotope_error"] == 0) & (g["source"] == 0)
    g = {k: v[m0] for k, v in g.items()}
    store_idx = (g["spectrum_id"].astype(np.uint64) & np.uint64(0xFFFFFFFF)).astype(np.int64)
    sp_index = np.asarray(sp.index, np.int64)
    order = np.argsort(sp_index, kind="stable")
    pos = np.searchsorted(sp_index[order], store_idx)
    posc = np.clip(pos, 0, len(order) - 1)
    ok = (pos < len(order)) & (sp_index[order][posc] == store_idx)
    spec = np.where(ok, order[posc], -1)
    keep = spec >= 0
    g = {k: v[keep] for k, v in g.items()}
    spec = spec[keep]
    n = len(sp)
    fid = np.nan_to_num(np.asarray(g["feature_id"], np.float64), nan=-1).astype(np.int64)
    f32 = lambda c: np.asarray(g[c], np.float32)  # noqa: E731
    agree = np.asarray(g["instrument_agrees"], bool)
    prim = {"feature_id": np.full(n, -1, np.int64), "isolation_fraction": np.zeros(n, np.float32),
            "feature_intensity": np.zeros(n, np.float32), "iso_cos": np.zeros(n, np.float32),
            "rt_apex": np.full(n, np.nan, np.float32), "nfw": np.zeros(n, np.uint16)}
    nfw = np.asarray(g["n_features_in_window"], np.uint16)
    prim["nfw"][spec] = nfw                                       # repeated on every row of a spectrum
    a = np.flatnonzero(agree)
    prim["feature_id"][spec[a]] = fid[a]
    prim["isolation_fraction"][spec[a]] = f32("isolation_fraction")[a]
    prim["feature_intensity"][spec[a]] = f32("feature_intensity")[a]
    prim["iso_cos"][spec[a]] = f32("feature_iso_cos")[a]
    prim["rt_apex"][spec[a]] = f32("feature_rt_apex")[a]
    s = np.flatnonzero(~agree)
    crank = np.asarray(g["chimeric_rank"], np.int64)[s]
    o = np.lexsort((fid[s], crank, spec[s]))
    s = s[o]
    ss = spec[s]
    first = np.searchsorted(ss, ss, side="left")
    grp = (np.arange(len(ss)) - first + 1).astype(np.uint16)
    sec = {"spec": ss.astype(np.int64), "group": grp, "mz": np.asarray(g["precursor_mz"], np.float64)[s],
           "charge": np.asarray(g["charge"], np.int8)[s], "feature_id": fid[s],
           "isolation_fraction": f32("isolation_fraction")[s], "feature_intensity": f32("feature_intensity")[s],
           "iso_cos": f32("feature_iso_cos")[s], "rt_apex": f32("feature_rt_apex")[s], "nfw": nfw[s]}
    n_groups = np.ones(n, np.uint16)
    np.add.at(n_groups, ss, 1)
    return {"prim": prim, "sec": sec, "n_groups": n_groups}


def virtual_spectra(sp: Spectra, sec: dict, mass_offset: np.ndarray | None = None) -> Spectra:
    """One virtual spectrum per secondary group: a copy of the real spectrum's peaks with the group's precursor.
    mass_offset (optional, Da per REAL spectrum) shifts every group's precursor MASS (G2 offset nulls)."""
    idx = sec["spec"]
    lens = np.diff(sp.offsets)[idx]
    starts = sp.offsets[:-1][idx]
    gi = np.repeat(starts, lens) + (np.arange(lens.sum()) - np.repeat(np.cumsum(lens) - lens, lens))
    off = np.r_[0, np.cumsum(lens)].astype(np.int64)
    z = sec["charge"].astype(np.int8)
    mz = sec["mz"].astype(np.float64)
    if mass_offset is not None:
        mz = mz + mass_offset[idx] / z.astype(np.float64)
    return Spectra(sp.run, sp.scan[idx], sp.index[idx], sp.rt[idx], mz, z, sp.mz[gi], sp.intensity[gi], off)


def _append(tbl: pa.Table, cols: dict) -> pa.Table:
    for k, (v, typ) in cols.items():
        tbl = tbl.append_column(k, pa.array(v, typ))
    return tbl


def _ms1_cols(n, fid, iso_frac, fint, iso_cos, rt_apex, nfw, inst_z):
    return {"ms1_feature_id": (fid, pa.int64()), "ms1_isolation_fraction": (iso_frac, pa.float32()),
            "ms1_feature_intensity": (fint, pa.float32()), "ms1_iso_cos": (iso_cos, pa.float32()),
            "ms1_feature_rt_apex": (rt_apex, pa.float32()), "n_features_in_window": (nfw, pa.uint16()),
            "inst_charge": (inst_z, pa.int8())}


def search_chimeric(sp: Spectra, gi, l1, p: SearchParams, groups: dict, report: int | None = None,
                    sec_isotope_errors=SEC_ISOTOPE_ERRORS, primary=None, mass_offset=None):
    """Chimeric search of one run.  primary: optional precomputed closed PSM table (engine.search(sp, ...)[0]) to
    reuse.  Returns (table, stats dict)."""
    from vulcan.gpu import engine as E
    t0 = time.perf_counter()
    st = {}
    if primary is None:
        primary, s0, _ = E.search(sp, gi, l1, p, report)
        st["primary"] = {"gpu_s": sum(s0.gpu_ms.values()) / 1e3, "wall_s": s0.wall_s.get("total"), **s0.counts}
    sec = groups["sec"]
    prim = groups["prim"]
    ng = groups["n_groups"]
    s = primary.column("spec_idx").to_numpy()
    n0 = primary.num_rows
    tbl0 = _append(primary, {"group": (np.zeros(n0, np.uint16), pa.uint16()), "chim_rank": (np.ones(n0, np.uint16), pa.uint16()),
                             "n_groups": (ng[s], pa.uint16())})
    tbl0 = _append(tbl0, _ms1_cols(n0, prim["feature_id"][s], prim["isolation_fraction"][s],
                                   prim["feature_intensity"][s], prim["iso_cos"][s], prim["rt_apex"][s], prim["nfw"][s],
                                   sp.charge[s].astype(np.int8)))
    parts = [tbl0]
    if len(sec["spec"]):
        vsp = virtual_spectra(sp, sec, mass_offset)
        p2 = dataclasses.replace(p, isotope_errors=tuple(sec_isotope_errors))
        t2, s2, _ = E.search(vsp, gi, l1, p2, report)
        st["secondary"] = {"gpu_s": sum(s2.gpu_ms.values()) / 1e3, "wall_s": s2.wall_s.get("total"), **s2.counts}
        v = t2.column("spec_idx").to_numpy()
        real = sec["spec"][v]
        t2 = t2.set_column(t2.column_names.index("spec_idx"), "spec_idx", pa.array(real.astype(np.int64)))
        n2 = t2.num_rows
        t2 = _append(t2, {"group": (sec["group"][v], pa.uint16()),
                          "chim_rank": ((sec["group"][v] + 1).astype(np.uint16), pa.uint16()),
                          "n_groups": (ng[real], pa.uint16())})
        t2 = _append(t2, _ms1_cols(n2, sec["feature_id"][v], sec["isolation_fraction"][v], sec["feature_intensity"][v],
                                   sec["iso_cos"][v], sec["rt_apex"][v], sec["nfw"][v], sp.charge[real].astype(np.int8)))
        parts.append(t2)
    tbl = pa.concat_tables(parts)
    st["n_secondary_groups"] = int(len(sec["spec"]))
    st["n_spectra"] = len(sp)
    st["wall_s"] = time.perf_counter() - t0
    return tbl, st
