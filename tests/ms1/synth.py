"""Synthetic MS1 runs with known isotope-envelope features (truth for tests/ms1)."""
from __future__ import annotations

import numpy as np

from vulcan.ms1.averagine import NEUTRON, PROTON, envelope
from vulcan.ms1.features import MS1Data


def make_run(n_feat=200, n_scans=300, seed=0, noise_per_scan=200, ppm_sd=1.5, charges=(1, 2, 3, 4),
             mass_range=(800.0, 4500.0), min_sep_ppm=30.0, run_id=7, return_peaks=False):
    """MS1 scans (1 s apart) with Gaussian-elution averagine features + uniform noise.

    Truth columns: mono_mz, charge, mass, apex_scan, sigma (scans), height. Features are placed so that no two
    isotope peaks of different features fall within `min_sep_ppm` while co-eluting (unambiguous truth).
    """
    rng = np.random.default_rng(seed)
    rt = np.arange(n_scans) / 60.0
    truth = []
    occupied = []   # (mz, scan_lo, scan_hi)
    while len(truth) < n_feat:
        z = int(rng.choice(charges))
        mass = rng.uniform(*mass_range)
        mono = mass / z + PROTON
        if not 350 < mono < 1600:
            continue
        apex = rng.uniform(15, n_scans - 15)
        sig = rng.uniform(2.0, 5.0)
        lo, hi = apex - 4 * sig, apex + 4 * sig
        env = envelope(mass)
        n_iso = int(np.sum(env / env.max() > 0.05))
        mzs = mono + np.arange(n_iso) * NEUTRON / z
        clash = False
        for (m, a, b) in occupied:
            if b >= lo and a <= hi and np.any(np.abs(mzs - m) / m * 1e6 < min_sep_ppm):
                clash = True
                break
        if clash:
            continue
        for m in mzs:
            occupied.append((m, lo, hi))
        height = 10 ** rng.uniform(5, 8)
        truth.append((mono, z, mass, apex, sig, height, n_iso))
    truth = np.array(truth, dtype=[("mono_mz", "f8"), ("charge", "i4"), ("mass", "f8"), ("apex_scan", "f8"),
                                   ("sigma", "f8"), ("height", "f8"), ("n_iso", "i4")])
    scans_mz, scans_int = [[] for _ in range(n_scans)], [[] for _ in range(n_scans)]
    for f in truth:
        env = envelope(f["mass"])
        env = env / env.max()
        s0, s1 = int(max(0, f["apex_scan"] - 4 * f["sigma"])), int(min(n_scans - 1, f["apex_scan"] + 4 * f["sigma"]))
        for s in range(s0, s1 + 1):
            prof = f["height"] * np.exp(-0.5 * ((s - f["apex_scan"]) / f["sigma"]) ** 2)
            for k in range(f["n_iso"]):
                it = prof * env[k]
                if it < 1e4:          # detection limit
                    continue
                mz = (f["mono_mz"] + k * NEUTRON / f["charge"]) * (1 + rng.normal(0, ppm_sd) * 1e-6)
                scans_mz[s].append(mz)
                scans_int[s].append(it * rng.uniform(0.85, 1.15))
    for s in range(n_scans):
        k = rng.poisson(noise_per_scan)
        scans_mz[s].extend(rng.uniform(350, 1600, k))
        scans_int[s].extend(10 ** rng.uniform(4, 5.5, k))
    off = np.zeros(n_scans + 1, np.int64)
    mz_l, it_l = [], []
    for s in range(n_scans):
        m = np.asarray(scans_mz[s]); i = np.asarray(scans_int[s])
        o = np.argsort(m, kind="stable")
        mz_l.append(m[o]); it_l.append(i[o]); off[s + 1] = off[s] + len(m)
    ms1 = MS1Data(run_id, np.arange(n_scans, dtype=np.int64) * 10, np.arange(n_scans, dtype=np.int64) * 10 + 1, rt,
                  np.concatenate(mz_l).astype(np.float32), np.concatenate(it_l).astype(np.float32), off)
    return ms1, truth


def match_truth(features, truth, ppm=10.0, scan_tol=None):
    """Greedy truth matching: same charge, mono within ppm, truth apex inside the feature's scan range."""
    f_mz = features.column("mono_mz").to_numpy()
    f_z = features.column("charge").to_numpy()
    f_s0 = features.column("ms1_start").to_numpy().astype(float)
    f_s1 = features.column("ms1_end").to_numpy().astype(float)
    hit_t = np.full(len(truth), -1)
    used = np.zeros(len(f_mz), bool)
    for i, t in enumerate(truth):
        c = np.flatnonzero((f_z == t["charge"]) & (np.abs(f_mz - t["mono_mz"]) / t["mono_mz"] * 1e6 <= ppm)
                           & (f_s0 - 1 <= t["apex_scan"]) & (f_s1 + 1 >= t["apex_scan"]) & ~used)
        if len(c):
            hit_t[i] = c[0]
            used[c[0]] = True
    return hit_t, used
