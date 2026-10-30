"""Amendment A3.2.2: synthetic secondary hypotheses on A spectra (pool "A-syn"), label-free and outcome-independent.

Docstring = specification (A3.7; stricter reading where the text is open; listed in scripts/a3/IMPLEMENTATION_NOTES.md).

Offsets table.  config/a3_secondary_offsets.tsv (built by scripts/a3/build_offsets.py from the D0 and D1 vulcan.ms1
  hypothesis tables, label-free, sha256 in config/a3.yaml): one row per dev secondary precursor group, columns
  delta_mz (Th, %.6f) and charge, rows in the committed order (sorted by delta_mz, charge, ...).  The empirical CDF is
  the row order: a draw u in [0, 1) selects row floor(u * N).
Draw.  For A spectrum key K (label-free: "<source run>|<store spectrum index>", independent of the mixture), chimeric
  rank r in {2, 3} and try t = 0, 1, ..., 31 (the "hash chain"):
      u = int(sha256(f"{K}|{r}|{seed}|{t}")[:16 hex], 16) / 2**64,   seed = 47 (config/a3.yaml)
  row = floor(u * N); the hypothesis is (mz' = mz_A + delta_mz[row], z' = charge[row]) where mz_A is the A spectrum's
  precursor m/z as carried in the mixture (= the certified mono m/z, A2.2(e)).  The first try that is not excluded
  is used; if all 32 are excluded the spectrum has no synthetic hypothesis of rank r.
Exclusions (a draw is excluded if ANY holds):
  (env) Isotope envelope of the certified precursor.  M = (mz_A - PROTON) z_A.  The draw is excluded if
        (i)  any searched neutral mass of the draw, M' - i NEUTRON with M' = (mz' - PROTON) z' and i in the searched
             isotope set {0, 1} (E7 sec_isotope_errors), lies within +-20 ppm of M + k NEUTRON for some k in -1..+3;
             or
        (ii) mz' lies within +-20 ppm of an envelope peak (M + k NEUTRON) / zz + PROTON for some k in -1..+3 and some
             charge zz in 1..6 ("at any z <= 6").
        Both readings are applied (stricter).
  (b)   Identity screen at W' (A2.2(b), unchanged code vulcan.gates.a2_identity.max_f_in_window): W' = (mz', z'),
        +-10 ppm, isotope errors {0, 1, 2} (the A2 screen convention: wider than the searched {0, 1}); excluded iff
        max_c f(c) >= 0.20.  The Hamming clause of A2.2(b) does not depend on the window and was already passed by
        every A spectrum; it is not recomputed.
  (d)   Family screen at W' (A2.2(d)): excluded iff the max exhaustive hyperscore of the A spectrum's peaks over
        (T∪D) ∩ W' (same W') is >= hyperscore(p) of the certified peptide (ties excluded).
  The screens are label-free functions of the unordered T∪D set (P1).  The callback order is env -> b -> d (the
  exclusion set does not depend on the order).
A-syn units are searched as separate virtual spectra and never enter the mixture's competition (step 4b).
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass

import numpy as np

PROTON = 1.0072764
NEUTRON = 1.00335
SEED = 47
RANKS = (2, 3)
MAX_TRIES = 32
ENV_PPM = 20.0
ENV_K = (-1, 0, 1, 2, 3)
ENV_ZMAX = 6
SEARCHED_ISO = (0, 1)
W_ISO = (0, 1, 2)
W_PPM = 10.0


@dataclass
class Offsets:
    delta_mz: np.ndarray
    charge: np.ndarray
    sha256: str

    def __len__(self):
        return self.delta_mz.size


def load_offsets(path: str, expected_sha256: str | None = None) -> Offsets:
    raw = open(path, "rb").read()
    sha = hashlib.sha256(raw).hexdigest()
    if expected_sha256 is not None and sha != expected_sha256:
        raise RuntimeError(f"A3 offsets table sha256 {sha} != recorded {expected_sha256}")
    lines = raw.decode().splitlines()
    hdr = lines[0].split("\t")
    i_d, i_z = hdr.index("delta_mz"), hdr.index("charge")
    d = np.array([float(x.split("\t")[i_d]) for x in lines[1:]], np.float64)
    z = np.array([int(x.split("\t")[i_z]) for x in lines[1:]], np.int64)
    return Offsets(d, z, sha)


def draw_u(spec_key: str, rank: int, t: int, seed: int = SEED) -> float:
    h = hashlib.sha256(f"{spec_key}|{rank}|{seed}|{t}".encode()).hexdigest()
    return int(h[:16], 16) / 2.0 ** 64


def draw_row(spec_key: str, rank: int, t: int, n_rows: int, seed: int = SEED) -> int:
    return min(int(draw_u(spec_key, rank, t, seed) * n_rows), n_rows - 1)


def envelope_excluded(mz_a: float, z_a: int, mz_s: float, z_s: int, ppm: float = ENV_PPM, ks=ENV_K,
                      zmax: int = ENV_ZMAX, searched_iso=SEARCHED_ISO) -> bool:
    M = (mz_a - PROTON) * z_a
    env = np.array([M + k * NEUTRON for k in ks])
    Ms = (mz_s - PROTON) * z_s
    for i in searched_iso:
        m = Ms - i * NEUTRON
        if np.any(np.abs(m - env) <= ppm * 1e-6 * env):
            return True
    for zz in range(1, zmax + 1):
        pk = env / zz + PROTON
        if np.any(np.abs(mz_s - pk) <= ppm * 1e-6 * pk):
            return True
    return False


def construct(spec_keys, mz_a, z_a, offsets: Offsets, screen, seed: int = SEED, ranks=RANKS,
              max_tries: int = MAX_TRIES) -> dict:
    """A-syn hypotheses for every A spectrum.

    spec_keys, mz_a, z_a: per A spectrum (label-free key, mixture precursor m/z = certified mono m/z, certified z).
    screen(i, mz_s, z_s) -> (excluded_b: bool, excluded_d: bool): the label-free W' screens for A spectrum i.
    Returns dict of arrays, one row per (A spectrum, rank) with a hypothesis: a_idx, rank, mz, charge, tries, row;
    plus counts of every exclusion reason (label-free)."""
    out = {k: [] for k in ("a_idx", "rank", "mz", "charge", "tries", "row")}
    cnt = {"env": 0, "b": 0, "d": 0, "none_after_max_tries": 0, "accepted": 0}
    n = len(offsets)
    for i, key in enumerate(spec_keys):
        for r in ranks:
            got = False
            for t in range(max_tries):
                row = draw_row(str(key), r, t, n, seed)
                mz_s = float(mz_a[i]) + float(offsets.delta_mz[row])
                z_s = int(offsets.charge[row])
                if envelope_excluded(float(mz_a[i]), int(z_a[i]), mz_s, z_s):
                    cnt["env"] += 1
                    continue
                eb, ed = screen(i, mz_s, z_s)
                if eb:
                    cnt["b"] += 1
                    continue
                if ed:
                    cnt["d"] += 1
                    continue
                for k, v in (("a_idx", i), ("rank", r), ("mz", mz_s), ("charge", z_s), ("tries", t + 1), ("row", row)):
                    out[k].append(v)
                cnt["accepted"] += 1
                got = True
                break
            if not got:
                cnt["none_after_max_tries"] += 1
    res = {"a_idx": np.asarray(out["a_idx"], np.int64), "rank": np.asarray(out["rank"], np.int64),
           "mz": np.asarray(out["mz"], np.float64), "charge": np.asarray(out["charge"], np.int64),
           "tries": np.asarray(out["tries"], np.int64), "row": np.asarray(out["row"], np.int64)}
    res["counts"] = cnt
    return res


def table_hash(t: dict) -> str:
    """sha256 of the A-syn hypothesis table content (keys + arrays), order as produced."""
    h = hashlib.sha256()
    for k in ("a_idx", "rank", "mz", "charge", "tries", "row"):
        h.update(k.encode())
        h.update(np.ascontiguousarray(t[k]).tobytes())
    return h.hexdigest()


def identity_screen_fn(db, p_seq, p_mask, cache=None, f_threshold: float = 0.20, ppm_frag: float = 20.0,
                       ppm_prec: float = W_PPM, isotopes=W_ISO):
    """(b) at W' with the unchanged A2 code: returns f(i, mz_s, z_s) -> bool (excluded)."""
    from vulcan.gates.a2_identity import max_f_in_window
    cache = {} if cache is None else cache

    def f(i, mz_s, z_s):
        return max_f_in_window(db, str(p_seq[i]), int(p_mask[i]), mz_s, z_s, cache, ppm_frag, ppm_prec,
                               tuple(isotopes)) >= f_threshold
    return f
