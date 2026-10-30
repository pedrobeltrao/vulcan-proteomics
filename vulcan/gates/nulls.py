"""Null-spectrum helpers for the G2 mixture design (DESIGN §6.4, A1.3)."""
from __future__ import annotations

import numpy as np


def quantile_map(src_rt, host_rt) -> np.ndarray:
    """A1.3 RT mapping: each source retention time goes to the host-run RT at the same empirical quantile
    (mid-rank of the source MS2 RT distribution -> host MS2 RT quantile, linear interpolation).  Monotone, so the
    within-source elution order is kept; the mapped distribution equals the host's."""
    s = np.asarray(src_rt, float)
    h = np.sort(np.asarray(host_rt, float))
    if s.size == 0:
        return s.copy()
    order = np.argsort(s, kind="stable")
    u = np.empty(s.size)
    u[order] = (np.arange(s.size) + 0.5) / s.size
    # ties in the source get the mean of their quantiles
    _, inv, cnt = np.unique(s, return_inverse=True, return_counts=True)
    if (cnt > 1).any():
        u = np.bincount(inv, weights=u) [inv] / cnt[inv]
    return np.interp(u, (np.arange(h.size) + 0.5) / h.size, h)
