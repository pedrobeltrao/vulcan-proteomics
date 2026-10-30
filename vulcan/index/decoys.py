"""Peptide-level decoys and paired shuffled entrapment for the L1 index (DESIGN.md §4 L1, §6.3).

The recipes are those of scripts/entrapment/peptides.py (M1e databases), re-implemented here so that the
permutation of every decoy/entrapment is known (vulcan_core.expand_pairs carries variable-mod sites
through it). tests/index/test_decoys.py checks sequence-identity with peptides.py.

* decoy  = pseudo-reverse with the protease C-term residue fixed; if it is a palindrome or collides with any
           forward peptide (T or E) or an earlier decoy, hash(peptide, DECOY_SEED, k) C-term-fixed shuffles.
* entrapment E_s (paired, r = 1) = hash(target, s, k) C-term-fixed shuffle of each target, re-drawn until it
           is new (not in T, not an earlier entrapment, not a pseudo-reverse of T, not a palindrome, its own
           pseudo-reverse not a forward peptide). Same checks as build_dbs.arm_paired's repair loop.

Everything is processed in sorted order and seeded by hash(sequence, seed, attempt), so the result does not
depend on FASTA/database order.
"""
from __future__ import annotations

import hashlib
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np

DECOY_SEED = 1001
STD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def pseudo_reverse(pep: str) -> str:
    return pep[-2::-1] + pep[-1]


def _rng(*keys) -> np.random.Generator:
    h = hashlib.blake2b("|".join(map(str, keys)).encode(), digest_size=8).digest()
    return np.random.default_rng(int.from_bytes(h, "little"))


def hash_shuffle_perm(pep: str, seed: int, attempt: int = 0) -> np.ndarray:
    """Permutation p (dec[j] = pep[p[j]]) of the C-term-fixed shuffle seeded by hash(pep, seed, attempt).
    Identical draws to scripts/entrapment/peptides.hash_shuffle (same Generator, same array length)."""
    n = len(pep)
    idx = np.arange(n - 1)
    _rng(pep, seed, attempt).shuffle(idx)
    return np.concatenate([idx, [n - 1]]).astype(np.uint8)


def apply_perm(pep: str, perm) -> str:
    return "".join(pep[i] for i in perm)


def hash_shuffle(pep: str, seed: int, attempt: int = 0) -> str:
    return apply_perm(pep, hash_shuffle_perm(pep, seed, attempt))


def _first_shuffles(args):
    peps, seed = args
    return [hash_shuffle_perm(p, seed, 0) for p in peps]


def _parallel_first(peps: list[str], seed: int, workers: int) -> list[np.ndarray]:
    if workers <= 1 or len(peps) < 20000:
        return _first_shuffles((peps, seed))
    step = (len(peps) + workers * 4 - 1) // (workers * 4)
    parts = [(peps[i:i + step], seed) for i in range(0, len(peps), step)]
    out: list[np.ndarray] = []
    import multiprocessing as mp
    with ProcessPoolExecutor(workers, mp_context=mp.get_context("forkserver")) as ex:
        for r in ex.map(_first_shuffles, parts):
            out.extend(r)
    return out


def make_entrapment(targets: list[str], seed: int, taken: set[str], workers: int = 1, max_tries: int = 200):
    """Paired shuffled entrapment for `targets` (sorted, I/L collapsed).

    taken: forward sequences that an entrapment must not equal (T and entrapments of earlier seeds); it is
    updated in place with the new entrapments. Returns ({target: (entrapment, perm)}, n_redrawn, n_dropped).
    """
    targets = sorted(targets)
    tset = set(targets)
    rev_t = {pseudo_reverse(t) for t in targets}
    first = _parallel_first(targets, seed, workers)
    out: dict[str, tuple[str, np.ndarray]] = {}
    n_redraw = n_drop = 0

    def bad(e: str) -> bool:
        pe = pseudo_reverse(e)
        return (e in taken or e in tset or e in rev_t or pe == e or pe in taken or pe in tset)

    for t, perm in zip(targets, first):
        e = apply_perm(t, perm)
        k = 0
        while bad(e) or e == t:
            k += 1
            if k > max_tries:
                e = None
                break
            perm = hash_shuffle_perm(t, seed, k)
            e = apply_perm(t, perm)
        if e is None:
            n_drop += 1
            continue
        if k:
            n_redraw += 1
        out[t] = (e, perm)
        taken.add(e)
    return out, n_redraw, n_drop


def make_decoys(forward: list[str], seed: int = DECOY_SEED, max_tries: int = 50):
    """Pseudo-reverse decoys (C-term fixed) for all forward peptides; collisions -> hash shuffles.
    Same rule as scripts/entrapment/peptides.build_decoys. Returns ({fwd: (decoy, perm or None)}, n_shuffled).
    perm None = pseudo-reverse."""
    fwd = set(forward)
    used: set[str] = set()
    out: dict[str, tuple[str, np.ndarray | None]] = {}
    n_repl = 0
    for p in sorted(forward):
        d, perm = pseudo_reverse(p), None
        k = 0
        while d == p or d in fwd or d in used:
            if k >= max_tries:
                d = None
                break
            perm = hash_shuffle_perm(p, seed, k)
            d = apply_perm(p, perm)
            k += 1
        if k:
            n_repl += 1
        if d is not None:
            out[p] = (d, perm)
            used.add(d)
    return out, n_repl


def default_workers() -> int:
    try:
        from vulcan.hw import n_threads
        return n_threads()
    except Exception:
        return len(os.sched_getaffinity(0))


def make_shuffled_decoys(forward: list[str], seed: int, used: set[str], workers: int = 1, max_tries: int = 50):
    """One shuffled-decoy set (DESIGN §6.2 variability arm, DECOY_SEEDS): hash(peptide, seed, k) C-term-fixed shuffle
    of every forward peptide (T and E), re-drawn (k+1) while it equals the peptide, any forward peptide, or any decoy
    already drawn (this set or an earlier set; `used` is updated in place). Sorted order + hash seeding => order-free.
    Returns ({fwd: (decoy, perm)}, n_redrawn, n_dropped)."""
    fwd = set(forward)
    forward = sorted(forward)
    first = _parallel_first(forward, seed, workers)
    out: dict[str, tuple[str, np.ndarray]] = {}
    n_redraw = n_drop = 0
    for p, perm in zip(forward, first):
        d = apply_perm(p, perm)
        k = 0
        while d == p or d in fwd or d in used:
            k += 1
            if k > max_tries:
                d = None
                break
            perm = hash_shuffle_perm(p, seed, k)
            d = apply_perm(p, perm)
        if d is None:
            n_drop += 1
            continue
        n_redraw += k > 0
        out[p] = (d, perm)
        used.add(d)
    return out, n_redraw, n_drop
