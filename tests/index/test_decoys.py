"""L1 decoy/entrapment recipes must equal the M1e recipes (scripts/entrapment/peptides.py) and be order-free."""
import os
import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "entrapment"))

import peptides as M1E  # noqa: E402

from vulcan.index import decoys as D  # noqa: E402

PEPS = ["PEPTIDEK", "LLLLLLLK", "AGCDEFGHLKR", "MSLVHLPHALMR", "VSLMWRGRNFR", "AAAAAAAR", "QWERTYLPASDFK"]


def test_hash_shuffle_identical_to_m1e():
    for p in PEPS:
        for seed in (1001, 2001, 2002):
            for k in range(5):
                assert D.hash_shuffle(p, seed, k) == M1E.hash_shuffle(p, seed, k)


def test_perm_semantics():
    for p in PEPS:
        perm = D.hash_shuffle_perm(p, 2001, 0)
        assert perm[-1] == len(p) - 1 and sorted(perm) == list(range(len(p)))
        assert D.apply_perm(p, perm) == M1E.hash_shuffle(p, 2001, 0)
    assert D.pseudo_reverse("PEPTIDEK") == M1E.pseudo_reverse("PEPTIDEK") == "EDITPEPK"


def test_decoys_match_m1e_and_order_free():
    fwd = sorted(set(PEPS))
    ours, n1 = D.make_decoys(fwd[::-1])
    ref, n2 = M1E.build_decoys({p: "T" for p in fwd}, D.DECOY_SEED)
    assert {k: v[0] for k, v in ours.items()} == ref and n1 == n2


def test_entrapment_is_new_and_order_free():
    t = sorted(PEPS)
    a, _, _ = D.make_entrapment(t, 2001, set(t))
    b, _, _ = D.make_entrapment(t[::-1], 2001, set(t))
    assert {k: v[0] for k, v in a.items()} == {k: v[0] for k, v in b.items()}
    for tg, (e, perm) in a.items():
        assert e not in t and e[-1] == tg[-1] and sorted(e) == sorted(tg)
        assert D.apply_perm(tg, perm) == e
