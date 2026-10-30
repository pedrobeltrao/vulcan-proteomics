"""vulcan.calib: offset curves, tolerance rules, pilot helpers (no search engine needed)."""
import numpy as np
import pyarrow as pa

from vulcan.calib.calibration import Calibration, Curve, OffsetModel, core_tolerance, fragment_tolerance
from vulcan.calib.pilot import (H2O, PROTON, REPORTERS, offset_histogram, parse_peptide, reporter_check, subsample)
from vulcan.calib.searcher import PSM_SCHEMA, SpectraBatch, read_sage, write_mgf


def test_offset_model_recovers_rt_drift():
    rng = np.random.default_rng(0)
    rt = rng.uniform(0, 120, 3000)
    mz = rng.uniform(400, 1200, 3000)
    truth = 3.0 * np.sin(rt / 30.0) + 1.0
    e = truth + rng.normal(0, 1.0, 3000)
    m = OffsetModel.fit(rt, mz, e)
    assert m.used["rt"]
    assert np.median(np.abs(m(rt, mz) - truth)) < 0.4
    assert m.mad_after < 0.7 * m.mad_before
    m2 = OffsetModel.from_dict(m.to_dict())
    assert np.allclose(m2(rt, mz), m(rt, mz))


def test_offset_model_constant_when_no_structure():
    rng = np.random.default_rng(1)
    rt, mz = rng.uniform(0, 100, 500), rng.uniform(300, 1500, 500)
    e = 2.0 + rng.normal(0, 1, 500)
    m = OffsetModel.fit(rt, mz, e)
    assert not m.used["mz"]
    assert abs(np.median(m(rt, mz)) - 2.0) < 0.2


def test_core_tolerance_ignores_flat_tail():
    rng = np.random.default_rng(2)
    e = np.r_[rng.normal(0.5, 1.0, 2000), rng.uniform(-50, 50, 200)]
    tol, info = core_tolerance(e, k=4)
    assert abs(info["mu"] - 0.5) < 0.1 and abs(info["sd"] - 1.0) < 0.15
    assert 3.4 < tol[1] - 0.5 < 4.6
    cov, vi = fragment_tolerance(e, 50.0, coverage=0.99)
    assert cov[1] - cov[0] < 10      # background subtraction removes the flat tail


def test_calibration_correction_roundtrip(tmp_path):
    cal = Calibration(run_id=3, run_name="x", precursor=OffsetModel(Curve([0, 100], [2.0, 6.0]), Curve()),
                      fragment=OffsetModel(Curve([50.0], [-3.0]), Curve()))
    mz = np.array([500.0, 1000.0])
    rt = np.array([0.0, 100.0])
    true = mz / (1 + np.array([2e-6, 6e-6]))
    assert np.allclose(cal.correct_precursor(mz, rt), true)
    assert np.allclose(cal.correct_fragment(mz, rt), mz / (1 - 3e-6))
    p = str(tmp_path / "c.json")
    cal.to_json(p)
    c2 = Calibration.from_json(p)
    assert np.allclose(c2.correct_precursor(mz, rt), true)


def test_parse_peptide_sage_notation():
    res, m, nt = parse_peptide("[+42.0106]-AC[+57.0215]M[+15.9949]K")
    assert res == "ACMK" and abs(nt - 42.0106) < 1e-9
    assert abs(m.sum() - (71.03711379 + 103.0091845 + 57.0215 + 131.0404846 + 15.9949 + 128.0949630)) < 1e-6
    _, m2, _ = parse_peptide("PEPTIDEK")
    assert abs(m2.sum() + H2O - 927.4549) < 0.01


def _batch(spectra, prec=None):
    off = np.cumsum([0] + [len(s) for s in spectra]).astype(np.int64)
    mz = np.concatenate([np.sort(np.asarray(s, float)) for s in spectra])
    n = len(spectra)
    return SpectraBatch(np.arange(n, dtype=np.uint64), np.zeros(n), np.full(n, 500.0) if prec is None else prec,
                        np.full(n, 2), mz, np.ones_like(mz), off)


def test_reporter_check_positive_and_negative():
    rng = np.random.default_rng(3)
    tmt = [list(REPORTERS["TMT10/11"]) + list(rng.uniform(150, 1500, 50)) for _ in range(200)]
    plain = [list(rng.uniform(110, 1500, 80)) for _ in range(200)]
    r = reporter_check({2: _batch(tmt), 3: _batch(plain)}, 10.0)
    assert r["ms2"]["detected"]["TMT"] and not r["ms2"]["detected"]["iTRAQ"]
    assert not r["ms3"]["detected"]["TMT"]


def test_offset_histogram_finds_injected_offsets():
    rng = np.random.default_rng(4)
    calc = rng.uniform(800, 3000, 4000)
    d = np.r_[np.zeros(2000), np.full(400, 15.9949), np.full(300, 79.9663), rng.uniform(-150, 500, 1300)]
    d = d + rng.normal(0, 0.002, len(d))
    t = pa.table({"exp_mass": calc + d, "calc_mass": calc})
    h = offset_histogram(t, 0.01)
    ann = [p["annotation"] for p in h["peaks"][:3]]
    assert ann[0] == "unmodified" and "oxidation" in ann and "phospho" in ann


def test_subsample_seeded():
    a = subsample(1000, 0.1, 11)
    assert len(a) == 100 and np.array_equal(a, subsample(1000, 0.1, 11)) and np.all(np.diff(a) > 0)


def test_mgf_and_sage_reader(tmp_path):
    b = _batch([[200.0, 300.0], [250.0]], prec=np.array([600.5, np.nan]))
    p = str(tmp_path / "x.mgf")
    assert write_mgf(b, p) == 1
    txt = open(p).read()
    assert "TITLE=0" in txt and "PEPMASS=600.5000000" in txt and "CHARGE=2+" in txt
    tsv = tmp_path / "results.sage.tsv"
    tsv.write_text("psm_id\tpeptide\tproteins\tscannr\tlabel\texpmass\tcalcmass\tcharge\tisotope_error\thyperscore\t"
                   "spectrum_q\tpeptide_q\trt\tmatched_peaks\n"
                   "1\tPEPTIDEK\tsp|X\t0\t1\t928.46\t927.4549\t2\t1.00335\t30.0\t0.001\t0.002\t1.5\t8\n")
    t = read_sage(str(tsv))
    assert t.select(PSM_SCHEMA.names).schema == PSM_SCHEMA
    assert t.column("isotope_error").to_pylist() == [1] and t.column("is_decoy").to_pylist() == [False]
    assert abs(PROTON - 1.00727646) < 1e-6


def test_sage_config_window_convention():
    from vulcan.calib.searcher import SearchConfig, sage_config
    c = sage_config(SearchConfig(fasta="x", precursor_ppm=(-2.0, 6.0), fragment_ppm=(-5.0, 9.0)), "a.mgf", "o")
    assert c["precursor_tol"] == {"ppm": [-6.0, 2.0]}      # error = exp - calc in [-2, 6]
    assert c["fragment_tol"] == {"ppm": [-9.0, 9.0]}
    c = sage_config(SearchConfig(fasta="x", precursor_da=(-150.0, 500.0)), "a.mgf", "o")
    assert c["precursor_tol"] == {"da": [-500.0, 150.0]}
