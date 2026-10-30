"""Averagine isotope envelopes (nominal-mass bins), tabulated on a mass grid.

Averagine (Senko 1995): C4.9384 H7.7583 N1.3577 O1.4773 S0.0417 per 111.1254 Da.
The table holds, for neutral masses 0..MAX_MASS in STEP-Da rows, the relative abundance of
isotopes M+0..M+(N_ISO-1), normalised to sum 1. Lookup is linear interpolation between rows.
"""
from __future__ import annotations

import numpy as np

PROTON = 1.007276466812
NEUTRON = 1.0033548378          # 13C - 12C; spacing used for isotope matching (as biosaur2)
AVERAGINE_MASS = 111.1254
AVERAGINE = {"C": 4.9384, "H": 7.7583, "N": 1.3577, "O": 1.4773, "S": 0.0417}
# natural isotope abundances by nominal-mass shift
ISOTOPES = {
    "C": [0.9893, 0.0107],
    "H": [0.999885, 0.000115],
    "N": [0.99636, 0.00364],
    "O": [0.99757, 0.00038, 0.00205],
    "S": [0.9499, 0.0075, 0.0425, 0.0, 0.0001],
}
N_ISO = 12
STEP = 10.0
MAX_MASS = 15000.0


def _pow_dist(p: np.ndarray, n: int, keep: int) -> np.ndarray:
    """Distribution of the sum of n iid draws from p (convolution power), truncated to keep bins."""
    out = np.zeros(keep)
    out[0] = 1.0
    base = np.zeros(keep)
    base[:min(len(p), keep)] = p[:keep]
    while n:
        if n & 1:
            out = np.convolve(out, base)[:keep]
        base = np.convolve(base, base)[:keep]
        n >>= 1
    return out


def envelope(mass: float, n_iso: int = N_ISO) -> np.ndarray:
    """Averagine isotope envelope (sum 1) for a neutral mass."""
    units = max(mass, 0.0) / AVERAGINE_MASS
    d = np.zeros(n_iso)
    d[0] = 1.0
    for el, per in AVERAGINE.items():
        n = int(round(per * units))
        if n:
            d = np.convolve(d, _pow_dist(np.asarray(ISOTOPES[el]), n, n_iso))[:n_iso]
    return d / d.sum()


def table(step: float = STEP, max_mass: float = MAX_MASS, n_iso: int = N_ISO) -> np.ndarray:
    masses = np.arange(0.0, max_mass + step, step)
    return np.stack([envelope(m, n_iso) for m in masses]).astype(np.float64)


_TABLE: np.ndarray | None = None


def get_table() -> np.ndarray:
    global _TABLE
    if _TABLE is None:
        _TABLE = table()
    return _TABLE
