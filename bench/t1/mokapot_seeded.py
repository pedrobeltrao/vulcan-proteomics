"""mokapot CLI with its --seed actually honoured (run with the mokapot env's python).

mokapot 0.10.0's CLI calls np.random.seed(--seed) but then brew() / the datasets use np.random.default_rng(None),
i.e. fresh OS entropy, so identical inputs give different results (D1: 37-68 peptides of ~6,030 at q <= 1% differ
between repeated runs).  Here default_rng(None) draws its seed from the (seeded) global numpy state; nothing else
changes.  Usage: python mokapot_seeded.py <mokapot CLI args>
"""
import numpy as np

_default_rng = np.random.default_rng


def _seeded_default_rng(seed=None):
    if seed is None:
        seed = int(np.random.randint(0, 2 ** 31 - 1))
    return _default_rng(seed)


np.random.default_rng = _seeded_default_rng

if __name__ == "__main__":
    from mokapot.mokapot import main
    main()
