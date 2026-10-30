"""Small shared helpers: seeding, segment operations, key encoding."""
from __future__ import annotations

import numpy as np

# Pre-registered seeds (DESIGN.md §6.1). Do not edit without a dated RESEARCH_TREE note.
SEEDS = (11, 23, 47, 83, 131)
DECOY_SEEDS = (1001, 1002, 1003, 1004, 1005)
ENTRAP_SEEDS_DEV = (2001, 2002, 2003)
ENTRAP_SEEDS_SEALED = (2004, 2005, 2006)
CALIB_DECOY_SEED = 3001
FOLD_SEED = 11


def make_rng(seed=None, *keys) -> np.random.Generator:
    """Deterministic Generator from a base seed plus optional integer sub-keys.

    ``make_rng(11, 3)`` and ``make_rng(11, 4)`` are independent streams.
    A Generator passed as ``seed`` is returned unchanged (keys ignored).
    """
    if isinstance(seed, np.random.Generator):
        return seed
    if seed is None:
        return np.random.default_rng()
    return np.random.default_rng(np.random.SeedSequence([int(seed), *[int(k) for k in keys]]))


def as_codes(x) -> np.ndarray:
    """Map an arbitrary key array to dense int64 codes (order of first sort)."""
    x = np.asarray(x)
    if x.dtype.kind in "iu" and x.ndim == 1:
        return x.astype(np.int64, copy=False)
    _, inv = np.unique(x, return_inverse=True)
    return inv.astype(np.int64).ravel()


def segment_starts(sorted_keys: np.ndarray) -> np.ndarray:
    """Boolean mask marking the first element of each run of equal keys."""
    n = sorted_keys.shape[0]
    m = np.ones(n, dtype=bool)
    if n > 1:
        m[1:] = sorted_keys[1:] != sorted_keys[:-1]
    return m


def segment_cumsum(values: np.ndarray, starts: np.ndarray) -> np.ndarray:
    """Cumulative sum that restarts at every ``starts`` position (inclusive)."""
    c = np.cumsum(values, dtype=np.int64 if values.dtype.kind in "biu" else np.float64)
    seg = np.cumsum(starts) - 1
    first = np.flatnonzero(starts)
    base = np.concatenate([[0], c])[first]  # cumulative value *before* each segment
    return c - base[seg]


def segment_reverse_cummin(values: np.ndarray, starts: np.ndarray, cap: float = 1.0) -> np.ndarray:
    """Reverse cumulative minimum within segments; values are clipped to [0, cap].

    Implemented with an offset trick (values clipped to ``cap`` so offsets of
    2*cap per segment separate segments exactly in float64).
    """
    v = np.minimum(values, cap)
    seg = np.cumsum(starts) - 1
    off = (2.0 * cap + 1.0) * seg.astype(np.float64)
    # reversed traversal visits later (higher-index) segments first; giving them a
    # larger offset guarantees their values never win the running minimum of an
    # earlier segment (all values are in [0, cap]).
    w = (v + off)[::-1]
    w = np.minimum.accumulate(w)[::-1]
    return w - off


def _splitmix64(x: np.ndarray) -> np.ndarray:
    with np.errstate(over="ignore"):
        z = x.astype(np.uint64) + np.uint64(0x9E3779B97F4A7C15)
        z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        return z ^ (z >> np.uint64(31))


def hash_tiebreak(keys, seed: int) -> np.ndarray:
    """Storage-order-free tie-break key = hash(key, seed) (DESIGN v0.3 §4 L1).

    Integer keys use a vectorised splitmix64; other keys (e.g. peptide
    sequences) use blake2b. Pass the result as ``tiebreak=`` to
    :func:`vulcan.stats.tdc.compete`; unlike a seeded RNG stream it does not
    depend on row order, so permuting the input never changes a winner.
    """
    k = np.asarray(keys)
    if k.dtype.kind in "iu":
        return _splitmix64(k.astype(np.uint64) ^ _splitmix64(np.array([seed], dtype=np.uint64))[0])
    import hashlib

    s = str(int(seed)).encode()
    return np.fromiter((int.from_bytes(hashlib.blake2b(str(x).encode() + b"|" + s, digest_size=8).digest(), "little")
                        for x in k.ravel()), dtype=np.uint64, count=k.size)
