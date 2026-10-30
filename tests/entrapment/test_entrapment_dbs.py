"""Tests for scripts/entrapment (DESIGN §6.1/§6.3): uniqueness vs all targets and all decoys of all decoy seeds,
C-term-fixed decoys, hash-seeded (order independent) shuffles, and the foreign-proteome homolog filter."""
import argparse
import os
import random
import shutil
import sys

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts", "entrapment"))
import build_dbs as B  # noqa: E402
import peptides as P  # noqa: E402

PFU = os.path.join(ROOT, "data", "fasta", "pfu_UP000001013.fasta")
HUMAN = os.path.join(ROOT, "data", "fasta", "human_sp_canonical.fasta")


def _head_fasta(src, dst, n):
    k = 0
    with open(src) as fi, open(dst, "w") as fo:
        for line in fi:
            if line.startswith(">"):
                k += 1
                if k > n:
                    break
            fo.write(line)
    return dst


def test_hash_shuffle_cterm_fixed_deterministic():
    for p in ["LSTPEPTIDEK", "AAAAAAAAK", "GGGGGGGR"]:
        for s in P.DECOY_SEEDS:
            a, b = P.hash_shuffle(p, s, 0), P.hash_shuffle(p, s, 0)
            assert a == b and a[-1] == p[-1] and sorted(a) == sorted(p)
    assert P.hash_shuffle("LSTPEPTIDEK", 2001, 0) != P.hash_shuffle("LSTPEPTIDEK", 2002, 0)


def test_build_decoys_order_independent_and_disjoint():
    rnd = random.Random(0)
    peps = {"".join(rnd.choice("ACDEFGHKLMNPQRSTVWY") for _ in range(8)) + "K": "T" for _ in range(500)}
    peps["ABBA"[:0] + "LSTTSLK"] = "T"          # reverse (C-term fixed) of LSTTSLK is itself -> palindrome
    d1, n1 = P.build_decoys(peps, 1001)
    items = list(peps.items())
    rnd.shuffle(items)
    d2, n2 = P.build_decoys(dict(items), 1001)
    assert d1 == d2 and n1 == n2 and n1 >= 1
    assert not (set(d1.values()) & set(peps))
    assert all(d[-1] == p[-1] for p, d in d1.items())


def test_isobaric_key_catches_ngg_qag():
    k = P.isobaric_key
    assert k("AANK") == k("AAGGK")          # N <-> GG
    assert k("AAQK") == k("AAAGK") == k("AAGAK")    # Q <-> AG / GA
    assert k("AGGK") == k("QGK")            # the greedy replacement used before missed this
    assert k("AAIK") == k("AALK")           # I/L
    assert k("AAAK") != k("AAGK")


def test_shuffle_decoy_set_matches_serial():
    seqs = ["LSTPEPTIDEK", "SAMPLEPEPR", "AVVVAAAK"]
    got = P.shuffle_decoy_set(seqs, procs=2)
    assert got == {P.hash_shuffle(s, d, 0) for s in seqs for d in P.DECOY_SEEDS}


@pytest.fixture(scope="module")
def tiny(tmp_path_factory):
    d = tmp_path_factory.mktemp("ent")
    t = _head_fasta(PFU, str(d / "t.fasta"), 60)
    f = _head_fasta(HUMAN, str(d / "f.fasta"), 300)
    return str(d), t, f


def _args(d, t, f):
    return argparse.Namespace(name="tiny", target=t, foreign=[f], out=d, seeds=[2001])


@pytest.mark.skipif(not os.path.exists(B.JAVA), reason="java env missing")
def test_paired_arm_unique_vs_all_decoys(tiny):
    d, t, f = tiny
    B.arm_paired(_args(d, t, f), 2001)
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(d, "tiny.P2001.peptides.tsv"))][1:]
    T = {r[0] for r in rows if r[1] == "T"}
    E = {r[0] for r in rows if r[1] == "E"}
    DT = {r[0] for r in rows if r[1] == "DT"}
    DE = {r[0] for r in rows if r[1] == "DE"}
    assert len(T) > 500 and len(E) > 0.95 * len(T)
    assert not (E & T) and not (DT & (T | E)) and not (DE & (T | E))
    # (a target whose own shuffle is a forward peptide is the decoy builder's business, not an entrapment constraint)
    assert not (E & P.shuffle_decoy_set(T)), "entrapment equals a shuffled decoy (seeds 1001-1005) of a target"
    assert not (P.shuffle_decoy_set(E) & (T | E)), "a shuffled decoy of an entrapment equals a forward peptide"
    for r in rows:
        if r[1] == "E":
            assert r[2] and sorted(r[0]) == sorted(r[2]) and r[0][-1] == r[2][-1] and r[0] != r[2]
        if r[1] in ("DT", "DE"):
            assert r[0][-1] == r[3][-1]


def test_foreign_arm_filter(tiny):
    d, t, f = tiny
    B.arm_foreign(_args(d, t, f))
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(d, "tiny.C.peptides.tsv"))][1:]
    T = {r[0] for r in rows if r[1] == "T"}
    E = {r[0] for r in rows if r[1] == "E"}
    assert E and not (E & T)
    tk = {P.isobaric_key(x) for x in T}
    assert all(P.isobaric_key(e) not in tk for e in E)
    heads = {(x[:i], x[i + 1:]) for x in T for i in range(len(x))}
    assert all(not any((e[:i], e[i + 1:]) in heads for i in range(len(e))) for e in E)
    assert not (E & P.shuffle_decoy_set(T))
    assert not (P.shuffle_decoy_set(E) & (T | E))


def test_shuffled_decoy_arm(tiny):
    d, t, f = tiny
    B.arm_shuffled(_args(d, t, f), 1002)
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(d, "tiny.S1002.peptides.tsv"))][1:]
    T = {r[0] for r in rows if r[1] == "T"}
    D = {r[0]: r[3] for r in rows if r[1] == "DT"}
    assert len(D) == len(T) and not (set(D) & T)
    assert all(k[-1] == v[-1] and sorted(k) == sorted(v) and k != v for k, v in D.items())
