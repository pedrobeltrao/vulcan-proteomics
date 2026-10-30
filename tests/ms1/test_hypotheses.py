"""Precursor-hypothesis table: pinned schema and semantics on a synthetic run with known co-isolation."""
import os
import sys

import numpy as np
import pyarrow as pa
import pytest

sys.path.insert(0, os.path.dirname(__file__))

from vulcan.calib.calibration import Calibration, Curve, OffsetModel  # noqa: E402
from vulcan.ms1.averagine import NEUTRON, PROTON, envelope  # noqa: E402
from vulcan.ms1.features import MS1Data, find_features  # noqa: E402
from vulcan.ms1.hypotheses import (HYPOTHESIS_SCHEMA, ISOTOPE_ERRORS, SOURCE, build_hypotheses,  # noqa: E402
                                   hypothesis_stats, read_table, write_table)

PINNED = [
    ("spectrum_id", "uint64", False), ("hypothesis_id", "uint16", False), ("precursor_group", "uint16", False),
    ("source", "uint8", False), ("charge", "int8", False), ("isotope_error", "int8", False),
    ("precursor_mz", "double", False), ("precursor_mz_cal", "double", False), ("neutral_mass", "double", False),
    ("isotope_error_prior", "float", False), ("feature_id", "uint32", True), ("feature_intensity", "float", True),
    ("feature_window_intensity", "float", True), ("isolation_fraction", "float", True),
    ("chimeric_rank", "uint16", True), ("feature_iso_cos", "float", True), ("feature_rt_apex", "float", True),
    ("instrument_agrees", "bool", False), ("instrument_isotope_offset", "int8", True),
    ("n_features_in_window", "uint16", False),
]


def test_schema_pinned():
    got = [(f.name, str(f.type), f.nullable) for f in HYPOTHESIS_SCHEMA]
    assert got == PINNED, "HYPOTHESIS_SCHEMA changed: update vulcan/ms1/SCHEMA.md and bump HYPOTHESIS_VERSION"
    assert SOURCE == {"ms1_feature": 0, "instrument": 1, "fallback": 2}
    assert ISOTOPE_ERRORS == (-1, 0, 1, 2, 3)


# species: (neutral mass, charge, height); A and B co-elute and overlap in m/z, C elutes alone
SPECIES = [(1600.80, 2, 1e8), (2402.90, 3, 3e7), (1100.55, 2, 5e7)]
N_SCANS = 40


def _ms1():
    off, mzs, its = [0], [], []
    for s in range(N_SCANS):
        row = []
        for i, (m, z, h) in enumerate(SPECIES):
            apex = 20.0 if i < 2 else 30.0
            prof = h * np.exp(-0.5 * ((s - apex) / 4.0) ** 2)
            env = envelope(m)
            env = env / env.max()
            for k in range(8):
                if env[k] > 0.03 and prof * env[k] > 1e3:
                    row.append((m / z + PROTON + k * NEUTRON / z, prof * env[k]))
        row.sort()
        mzs += [r[0] for r in row]
        its += [r[1] for r in row]
        off.append(len(mzs))
    # MS1 at store rows 0, 10, 20, ...; scans 1, 11, 21, ...
    return MS1Data(5, np.arange(N_SCANS) * 10, np.arange(N_SCANS) * 10 + 1, np.arange(N_SCANS) * 0.05,
                   np.asarray(mzs, np.float32), np.asarray(its, np.float32), np.asarray(off, np.int64))


def _ms2_table(rows):
    """rows: (store index, master MS1 ordinal, target m/z, instrument mono or nan, charge or 0, half width)."""
    n = len(rows)
    r = np.array(rows, dtype=float)
    idx = r[:, 0].astype(np.uint32)
    master = (r[:, 1] * 10 + 1).astype(np.int32)
    mono = r[:, 3]
    z = r[:, 4].astype(int)
    return pa.table({
        "spectrum_id": pa.array((np.uint64(5) << np.uint64(32)) | idx.astype(np.uint64)),
        "index": pa.array(idx), "ms_level": pa.array(np.full(n, 2, np.uint8)),
        "rt_min": pa.array(r[:, 1] * 0.05 + 0.01), "precursor_mz": pa.array(np.where(np.isnan(mono), r[:, 2], mono)),
        "isolation_target_mz": pa.array(r[:, 2]), "mono_mz": pa.array(mono, mask=np.isnan(mono)),
        "isolation_lower_offset": pa.array(r[:, 5].astype(np.float32)),
        "isolation_upper_offset": pa.array(r[:, 5].astype(np.float32)),
        "precursor_charge": pa.array(z.astype(np.int8), mask=z == 0),
        "master_scan": pa.array(master),
    })


mzA = SPECIES[0][0] / 2 + PROTON
mzB = SPECIES[1][0] / 3 + PROTON
mzC = SPECIES[2][0] / 2 + PROTON


@pytest.fixture(scope="module")
def built():
    ms1 = _ms1()
    fr = find_features(ms1)
    rows = [
        (101, 20, mzA + 0.2, mzA, 2, 1.5),          # A dominant, B co-isolated (chimeric), instrument agrees
        (102, 20, mzA + 0.5, mzA + NEUTRON / 2, 2, 1.0),   # instrument picked A's M+1 -> offset +1, feature mono wins
        (103, 30, mzC, np.nan, 0, 0.7),           # no instrument charge, feature C present
        (104, 5, 700.0, np.nan, 0, 0.7),          # nothing there: fallback z = 2, 3, 4
        (105, 5, 900.0, 900.0, 3, 0.7),           # instrument only
    ]
    tbl, st = build_hypotheses(_ms2_table(rows), fr)
    return fr, tbl, st


def _rows(tbl, idx):
    sid = int((np.uint64(5) << np.uint64(32)) | np.uint64(idx))
    m = tbl.column("spectrum_id").to_numpy() == sid
    return tbl.filter(pa.array(m)).to_pandas()


def test_features_found(built):
    fr, _, _ = built
    F = fr.features.to_pandas()
    for m, z, _ in SPECIES:
        assert ((F.charge == z) & (np.abs(F.neutral_mass - m) < 0.01)).sum() == 1


def test_chimeric_window(built):
    _, tbl, _ = built
    d = _rows(tbl, 101)
    g = d[d.isotope_error == 0].sort_values("precursor_group")
    assert list(g.source) == [0, 0]
    assert list(g.chimeric_rank) == [1, 2]
    assert list(g.charge) == [2, 3]
    assert abs(g.neutral_mass.iloc[0] - SPECIES[0][0]) < 0.01 and abs(g.neutral_mass.iloc[1] - SPECIES[1][0]) < 0.01
    assert g.instrument_agrees.iloc[0] and not g.instrument_agrees.iloc[1]
    assert g.instrument_isotope_offset.iloc[0] == 0
    assert np.isclose(g.isolation_fraction.sum(), 1.0, atol=1e-5)
    assert (g.n_features_in_window == 2).all()


def test_isotope_expansion_and_mass(built):
    _, tbl, _ = built
    d = _rows(tbl, 101)
    assert len(d) == 2 * len(ISOTOPE_ERRORS)
    assert list(d.hypothesis_id) == list(range(len(d)))
    for _, r in d.iterrows():
        exp = (r.precursor_mz_cal - PROTON) * r.charge - r.isotope_error * NEUTRON
        assert abs(r.neutral_mass - exp) < 1e-9
    for g, x in d.groupby("precursor_group"):
        assert list(x.isotope_error) == list(ISOTOPE_ERRORS)
        assert abs(x.isotope_error_prior.sum() - 1.0) < 1e-5


def test_instrument_isotope_offset(built):
    _, tbl, _ = built
    d = _rows(tbl, 102)
    a = d[(d.isotope_error == 0) & d.instrument_agrees]
    assert len(a) == 1
    assert a.instrument_isotope_offset.iloc[0] == 1
    assert abs(a.precursor_mz.iloc[0] - mzA) < 1e-3       # the feature mono, not the instrument's M+1
    assert not (d.source == SOURCE["instrument"]).any()


def test_no_charge_with_feature_and_fallback(built):
    _, tbl, _ = built
    d = _rows(tbl, 103)
    assert set(d.source) == {0} and set(d.charge) == {2}
    d = _rows(tbl, 104)
    assert set(d.source) == {SOURCE["fallback"]}
    assert sorted(set(d.charge)) == [2, 3, 4]
    assert len(d) == 3 * len(ISOTOPE_ERRORS)
    assert d.feature_id.isna().all()
    d = _rows(tbl, 105)
    assert set(d.source) == {SOURCE["instrument"]} and set(d.charge) == {3}


def test_calibration_applied(built):
    fr, _, _ = built
    cal = Calibration(precursor=OffsetModel(Curve([0.0], [5.0]), Curve()))   # +5 ppm everywhere
    tbl, _ = build_hypotheses(_ms2_table([(101, 20, mzA + 0.2, mzA, 2, 1.5)]), fr, calib=cal)
    d = tbl.to_pandas()
    assert np.allclose(d.precursor_mz_cal, d.precursor_mz / (1 + 5e-6))
    meta = tbl.schema.metadata[b"vulcan_hypotheses"]
    assert b'"calibration"' in meta


def test_roundtrip_and_stats(built, tmp_path):
    _, tbl, st = built
    p = str(tmp_path / "h.arrow")
    write_table(tbl, p)
    t2 = read_table(p)
    assert t2.equals(tbl)
    assert t2.schema.remove_metadata() == HYPOTHESIS_SCHEMA
    s = hypothesis_stats(tbl, st)
    assert s["n_ms2"] == 5 and s["n_rows"] == tbl.num_rows
    assert s["groups_by_source"]["fallback"] == 3
