"""MS1 feature finder on synthetic isotope envelopes with known truth."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from synth import make_run, match_truth  # noqa: E402

from vulcan.ms1 import hills as H  # noqa: E402
from vulcan.ms1.averagine import NEUTRON, PROTON, envelope  # noqa: E402
from vulcan.ms1.features import FEATURE_SCHEMA, FeatureParams, MS1Data, find_features  # noqa: E402


@pytest.fixture(scope="module")
def easy():
    ms1, truth = make_run(n_feat=150, n_scans=250, seed=1)
    return ms1, truth, find_features(ms1)


def test_schema(easy):
    _, _, fr = easy
    assert fr.features.schema.remove_metadata() == FEATURE_SCHEMA
    ids = fr.features.column("feature_id").to_numpy()
    assert np.array_equal(ids, np.arange(len(ids)))


def test_recall_precision_charge_mono(easy):
    _, truth, fr = easy
    hit, used = match_truth(fr.features, truth, ppm=5.0)
    assert (hit >= 0).mean() >= 0.97
    assert used.mean() >= 0.97          # noise peaks (uniform m/z) must not form features
    F = fr.features
    mz = F.column("mono_mz").to_numpy()[hit[hit >= 0]]
    t = truth["mono_mz"][hit >= 0]
    assert np.median(np.abs(mz - t) / t * 1e6) < 1.5   # 1.5 ppm simulated scatter


def test_apex_and_bounds(easy):
    ms1, truth, fr = easy
    hit, _ = match_truth(fr.features, truth)
    F = fr.features
    ap = F.column("ms1_apex").to_numpy()[hit[hit >= 0]].astype(float)
    assert np.median(np.abs(ap - truth["apex_scan"][hit >= 0])) <= 1.0
    assert np.all(F.column("rt_start").to_numpy() <= F.column("rt_apex").to_numpy() + 1e-6)


def test_peak_feature_map(easy):
    ms1, _, fr = easy
    pf = fr.peak_feature
    assert len(pf) == len(ms1.mz)
    F = fr.features
    # every feature owns its mono peaks: peaks mapped to iso 0 lie at the feature mono m/z
    m = (pf >= 0) & (fr.peak_iso == 0)
    mono = F.column("mono_mz").to_numpy()[pf[m]]
    assert np.median(np.abs(ms1.mz[m] - mono) / mono * 1e6) < 3


def test_linker_matches_numpy_reference():
    ms1, _ = make_run(n_feat=60, n_scans=80, seed=3, noise_per_scan=400)
    lens = np.diff(ms1.off)
    scan = np.repeat(np.arange(len(lens), dtype=np.int32), lens)
    mz = ms1.mz.astype(np.float64)
    ref = H.link_numpy(scan, mz, 8.0)
    soff = np.searchsorted(scan, np.arange(len(lens) + 1)).astype(np.int64)
    got = np.full(len(mz), -1, np.int64)
    H._link_adjacent(soff, mz, 8.0, got)
    assert np.array_equal(ref, got)
    assert np.array_equal(H.roots_numpy(ref), H._roots_seq(ref))


def _single(mass, z, n_scans=30, heights=(1e7,), apexes=(15.0,), sig=3.0, drop_iso=()):
    """One or more co-located elution peaks of the same species (for valley splitting / pattern tests)."""
    env = envelope(mass)
    env = env / env.max()
    mono = mass / z + PROTON
    off = [0]
    mzs, its = [], []
    for s in range(n_scans):
        prof = sum(h * np.exp(-0.5 * ((s - a) / sig) ** 2) for h, a in zip(heights, apexes))
        row = [(mono + k * NEUTRON / z, prof * env[k]) for k in range(len(env)) if env[k] > 0.02 and k not in drop_iso
               and prof * env[k] > 1e3]
        row.sort()
        mzs += [r[0] for r in row]
        its += [r[1] for r in row]
        off.append(len(mzs))
    return MS1Data(1, np.arange(n_scans) * 2, np.arange(n_scans) * 2 + 1, np.arange(n_scans) / 60.0,
                   np.asarray(mzs, np.float32), np.asarray(its, np.float32), np.asarray(off, np.int64))


def test_charge_harmonic_not_reported_as_lower_charge():
    fr = find_features(_single(2000.0, 2))
    F = fr.features
    assert F.num_rows == 1
    assert F.column("charge").to_numpy()[0] == 2
    assert abs(F.column("neutral_mass").to_numpy()[0] - 2000.0) < 0.01


def test_valley_split_two_elution_peaks():
    fr = find_features(_single(1500.0, 2, n_scans=60, heights=(1e7, 8e6), apexes=(15.0, 42.0), sig=3.0))
    F = fr.features
    assert F.num_rows == 2
    ap = np.sort(F.column("ms1_apex").to_numpy())
    assert abs(ap[0] - 15) <= 1 and abs(ap[1] - 42) <= 1


def test_missing_mono_shifts_by_one_isotope():
    # an envelope without its mono peak is reported (at most) one isotope higher, never at a wrong charge
    fr = find_features(_single(1200.0, 2, drop_iso=(0,)))
    F = fr.features
    if F.num_rows:
        assert F.column("charge").to_numpy()[0] == 2
        d = F.column("neutral_mass").to_numpy()[0] - 1200.0
        assert abs(d - NEUTRON) < 0.02


def test_explained_fraction_rejects_heavy_two_peak_pattern():
    # a 5 kDa species whose only visible peaks are M0 and M1 (< 50% of its averagine envelope)
    fr = find_features(_single(5000.0, 4, drop_iso=tuple(range(2, 12))))
    assert fr.features.num_rows == 0
    fr2 = find_features(_single(5000.0, 4, drop_iso=tuple(range(2, 12))), FeatureParams(min_explained=0.0))
    assert fr2.features.num_rows == 1


def test_deterministic(easy):
    ms1, _, fr = easy
    fr2 = find_features(ms1)
    assert fr.features.equals(fr2.features)
