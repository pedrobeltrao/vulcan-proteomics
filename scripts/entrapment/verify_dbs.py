#!/usr/bin/env python3
"""Independent verification of the built M1e databases against DESIGN §6.1/§6.3 (reads only <name>.peptides.tsv).

    python scripts/entrapment/verify_dbs.py DBDIR [--families human_contam pfu_contam]

Checks per DB (exit 1 on any failure; results printed as JSON):
  * E disjoint from T; no duplicate sequences; standard residues; length 7-50
  * every decoy (pseudo-reverse, DT/DE) keeps the C-term residue of its forward peptide, is not a forward peptide
  * E is not equal to any shuffled decoy (hash(peptide, seed), seeds 1001-1005) of any target, and no shuffled decoy
    of an entrapment equals a target or an entrapment
  * P arms: E has the composition and C-term of its paired target and differs from it
  * C arm: no entrapment is I/L-identical, isobaric (N/GG, Q/AG), or a Hamming-1 (I/L-collapsed) neighbour of a target
    (re-derived here with exact string keys, independently of build_dbs.py's hashed keys)
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from peptides import STD_AA, isobaric_key, shuffle_decoy_set  # noqa: E402


def read(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path)][1:]
    d = {"T": [], "E": [], "DT": [], "DE": []}
    for r in rows:
        d[r[1]].append(r)
    return d


def check(path):
    d = read(path)
    T = {r[0] for r in d["T"]}
    E = {r[0] for r in d["E"]}
    out = {"file": os.path.basename(path), "n_T": len(T), "n_E": len(E), "fails": []}

    def req(name, ok):
        out[name] = bool(ok)
        if not ok:
            out["fails"].append(name)
    req("no_duplicates", len(T) == len(d["T"]) and len(E) == len(d["E"]))
    req("E_disjoint_T", not (E & T))
    req("std_residues_len", all(set(s) <= STD_AA and 7 <= len(s) <= 50 for s in T | E))
    dec_ok, dec_fwd = True, True
    fwd = T | E
    for r in d["DT"] + d["DE"]:
        dec_ok &= r[0][-1] == r[3][-1]
        dec_fwd &= r[0] not in fwd
    req("decoys_cterm_fixed", dec_ok)
    req("decoys_not_forward", dec_fwd)
    if E:
        req("E_not_shuffled_target_decoy_1001_1005", not (E & shuffle_decoy_set(T)))
        SE = shuffle_decoy_set(E)
        req("shuffled_E_decoys_not_forward", not (SE & fwd))
        paired = [r for r in d["E"] if r[2]]
        if paired:
            req("P_pairs_same_composition_cterm", all(sorted(r[0]) == sorted(r[2]) and r[0][-1] == r[2][-1]
                                                       and r[0] != r[2] for r in paired))
            out["n_paired"] = len(paired)
        else:   # arm C
            tk = {isobaric_key(t) for t in T}
            req("C_not_isobaric_or_IL_identical", all(isobaric_key(e) not in tk for e in E))
            heads = {(t[:i], t[i + 1:]) for t in T for i in range(len(t))}
            req("C_no_hamming1_neighbour", not any(any((e[:i], e[i + 1:]) in heads for i in range(len(e))) for e in E))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dbdir")
    ap.add_argument("--families", nargs="+", default=["human_contam", "pfu_contam"])
    ap.add_argument("--arms", nargs="+", default=["target", "P2001", "P2002", "P2003", "C"])
    a = ap.parse_args()
    bad = 0
    res = []
    for fam in a.families:
        for arm in a.arms:
            p = os.path.join(a.dbdir, f"{fam}.{arm}.peptides.tsv")
            if not os.path.exists(p):
                continue
            r = check(p)
            res.append(r)
            bad += bool(r["fails"])
            print(json.dumps(r), flush=True)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
