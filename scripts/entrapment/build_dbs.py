#!/usr/bin/env python3
"""Build the M1e peptide-level search databases (DESIGN.md §4 L1, §6.1, §6.3).

All databases are *peptide-level* FASTAs (one entry per peptide, searched with no further cleavage), so
that every target, entrapment and decoy peptide is an explicit, labelled entry:

    DB = T  ∪  E  ∪  dec(T)  ∪  dec(E)

  T       tryptic peptides of the target FASTA (K/R not before P, <= 2 missed cleavages, length 7-50,
          I->L collapsed, standard residues only)
  E       entrapment peptides
  dec(.)  pseudo-reverse (C-term fixed); a palindrome/collision is replaced by hash(peptide, 1001) shuffles

Arms
  target   : T only (+ decoys). Control that the peptide-level, no-cleave set-up reproduces the
             protein-level production search.
  P<seed>  : FDRBench v1.1.1 paired shuffled entrapment (r = 1, C-term fixed, -I2L, -check, seed). Each
             FDRBench entrapment is then re-checked for uniqueness against T, all other E, dec(T) and
             dec(E); offenders are re-drawn by hash(target, seed, k) C-term-fixed shuffles (order
             independent, §6.1). Targets FDRBench could not pair are kept as unpaired targets.
  C        : natural foreign-proteome entrapment. Foreign peptides are removed if they are I/L-identical
             to a target, isobaric to one (N<->GG, Q<->AG/GA), a single-substitution neighbour of one
             (I/L-collapsed Hamming 1), or equal to a target decoy. r = |E| / |T| on the effective space
             (unique, I/L-collapsed, standard residues, length 7-50, CAM mass 500-5000 Da), also per
             target-mass decile.

Outputs per (dataset-db, arm) in --out:
  <name>.fasta          search DB (T_*, E_*, F_*, rev_* accessions)
  <name>.peptides.tsv   seq, kind (T/E/DT/DE), pair_key (target seq for paired arms), decoy_of
  <name>.json           counts, r, repairs, checksums
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from peptides import (DECOY_SEEDS, shuffle_decoy_set, shuffle_hits, STD_AA, accession, build_decoys, digest, hash_shuffle, il, in_space, isobaric_key,  # noqa: E402
                      mass, pseudo_reverse, read_fasta, sha256, write_peptide_db)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
JAVA = os.path.join(ROOT, "conda_envs", "java", "bin", "java")
FDRBENCH = os.path.join(ROOT, "tools", "fdrbench", "fdrbench-1.1.1", "fdrbench-1.1.1.jar")
DECOY_SEED = 1001


def target_peptides(fasta):
    """{I/L-collapsed peptide: set(accessions)} for standard-residue peptides; plus #dropped non-standard."""
    peps, nonstd = {}, set()
    for h, s in read_fasta(fasta):
        acc = accession(h)
        for p in digest(s):
            q = il(p)
            if not set(q) <= STD_AA:
                nonstd.add(q)
                continue
            peps.setdefault(q, set()).add(acc)
    return peps, nonstd


def write_outputs(out, name, T, E, pair_key, meta, kind_e="E"):
    """T: list of target seqs, E: list of entrapment seqs (same order as pair_key for paired arms)."""
    fwd = {t: "T" for t in T}
    for e in E:
        assert e not in fwd, e
        fwd[e] = kind_e
    dec, n_repl = build_decoys(fwd, DECOY_SEED)
    # final disjointness check (all I/L collapsed)
    D = set(dec.values())
    assert not (D & set(fwd)), "decoy equals a forward peptide"
    assert len(D) == len(dec), "duplicate decoys"
    missing = [p for p in fwd if p not in dec]
    fa = os.path.join(out, f"{name}.fasta")
    tsv = os.path.join(out, f"{name}.peptides.tsv")

    def entries():
        for i, t in enumerate(sorted(T)):
            yield f"T_{i:07d}", t
        for i, e in enumerate(sorted(E)):
            yield f"{kind_e}_{i:07d}", e
        for i, p in enumerate(sorted(dec)):
            yield f"rev_{fwd[p]}_{i:07d}", dec[p]
    write_peptide_db(fa, entries())
    with open(tsv, "w") as fh:
        fh.write("seq\tkind\tpair_key\tdecoy_of\n")
        for t in sorted(T):
            fh.write(f"{t}\tT\t{t if pair_key is not None else ''}\t\n")
        for e in sorted(E):
            fh.write(f"{e}\t{kind_e}\t{pair_key.get(e, '') if pair_key is not None else ''}\t\n")
        for p in sorted(dec):
            fh.write(f"{dec[p]}\tD{fwd[p]}\t\t{p}\n")
    meta.update({
        "name": name, "n_target": len(T), "n_entrapment": len(E), "n_decoy": len(dec),
        "n_decoy_replaced_by_shuffle": n_repl, "n_forward_without_decoy": len(missing),
        "decoy_recipe": f"pseudo-reverse C-term fixed; collisions/palindromes -> hash(peptide,{DECOY_SEED}) shuffle",
        "fasta": os.path.relpath(fa, ROOT), "fasta_sha256": sha256(fa),
        "peptides_tsv": os.path.relpath(tsv, ROOT), "peptides_tsv_sha256": sha256(tsv),
    })
    with open(os.path.join(out, f"{name}.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    print(json.dumps({k: v for k, v in meta.items() if not isinstance(v, (dict, list))}), flush=True)
    return meta


def effective_space(seqs):
    return [s for s in seqs if in_space(s)]


def r_by_decile(T, E):
    tm = np.array([mass(s) for s in effective_space(T)])
    em = np.array([mass(s) for s in effective_space(E)])
    edges = np.quantile(tm, np.linspace(0, 1, 11))
    edges[0], edges[-1] = -np.inf, np.inf
    ct = np.histogram(tm, edges)[0]
    ce = np.histogram(em, edges)[0]
    return {"r": float(len(em) / len(tm)), "n_target_eff": int(len(tm)), "n_entrapment_eff": int(len(em)),
            "decile_mass_edges_da": [float(x) for x in edges[1:-1]],
            "r_per_mass_decile": [float(e / t) for e, t in zip(ce, ct)]}


# ------------------------------------------------------------------------------------------- arms

def arm_target(args):
    peps, nonstd = target_peptides(args.target)
    T = sorted(peps)
    meta = {"arm": "target", "target_fasta": os.path.relpath(args.target, ROOT), "target_fasta_sha256": sha256(args.target),
            "n_target_nonstandard_dropped": len(nonstd), "n_target_eff": len(effective_space(T))}
    write_outputs(args.out, f"{args.name}.target", T, [], None, meta)


def arm_paired(args, seed):
    peps, nonstd = target_peptides(args.target)
    name = f"{args.name}.P{seed}"
    fd_out = os.path.join(args.out, "fdrbench", f"{name}.fdrbench.txt")
    os.makedirs(os.path.dirname(fd_out), exist_ok=True)
    cmd = [JAVA, "-Xmx48G", "-jar", FDRBENCH, "-level", "peptide", "-db", args.target, "-o", fd_out, "-I2L",
           "-fix_nc", "c", "-enzyme", "1", "-miss_c", "2", "-minLength", "7", "-maxLength", "50",
           "-seed", str(seed), "-check", "-fold", "1"]
    t0 = time.time()
    log = os.path.join(args.out, "fdrbench", f"{name}.fdrbench.log")
    with open(log, "w") as fh:
        subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.STDOUT, cwd=os.path.dirname(fd_out))
    t_fdrbench = time.time() - t0
    pairs = {}
    with open(fd_out) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pairs.setdefault(r["peptide_pair_index"], {})[r["peptide_type"]] = r["sequence"]
    fd_T = {p["target"] for p in pairs.values()}
    # digest equivalence: our digest (I/L) vs FDRBench's target list
    ours = set(peps) | nonstd
    only_ours, only_fd = ours - fd_T, fd_T - ours
    T = set(peps)                                   # all standard-residue targets are kept
    E_of = {}
    for p in pairs.values():
        t, e = p["target"], p["p_target"]
        if t in T:
            E_of[t] = e
    # uniqueness checks and order-independent repair
    Eset = Counter(E_of.values())
    revT = {pseudo_reverse(t) for t in T}
    # DESIGN 6.1: every entrapment peptide unique against all targets and all decoys of all decoy seeds
    ST = shuffle_decoy_set(T)
    E0 = set(E_of.values())
    SE0 = shuffle_decoy_set(E0)                    # shuffled decoys of the initial entrapment set
    SEadd = set()                                   # shuffled decoys of accepted replacement candidates
    coll0 = shuffle_hits(E0, T | E0)               # initial E whose own shuffled decoy hits a target / entrapment
    reasons = Counter()

    def bad(t, e, Ecur, revE, orig=False):
        if e in T:
            return "E_in_T"
        if e in ST:
            return "E_in_shuffle_decT"
        if orig:
            if e in coll0:
                return "E_shuffle_decoy_collides"
        elif e in SE0 or e in SEadd or any(hash_shuffle(e, s, 0) in T or hash_shuffle(e, s, 0) in Ecur for s in DECOY_SEEDS):
            return "E_shuffle_decoy_collides"
        if Ecur[e] > 1:
            return "E_duplicate"
        if e in revT:
            return "E_in_decT"
        if pseudo_reverse(e) == e:
            return "E_palindrome"
        if pseudo_reverse(e) in Ecur and pseudo_reverse(e) != e:
            return "E_in_decE"
        if not set(e) <= STD_AA:
            return "E_nonstandard"
        return None
    revE = None
    repaired = {}
    for t in sorted(E_of):
        e = E_of[t]
        why = bad(t, e, Eset, revE, orig=True)
        if why is None:
            continue
        reasons[why] += 1
        Eset[e] -= 1
        if Eset[e] == 0:
            del Eset[e]
        k, new = 0, None
        while k < 200:
            cand = hash_shuffle(t, seed, k)
            k += 1
            Eset[cand] += 1
            if bad(t, cand, Eset, revE) is None and cand != t:
                new = cand
                break
            Eset[cand] -= 1
            if Eset[cand] == 0:
                del Eset[cand]
        if new is None:
            del E_of[t]
            reasons["unrepairable_dropped_pair"] += 1
        else:
            E_of[t] = new
            SEadd.update(hash_shuffle(new, s, 0) for s in DECOY_SEEDS)
            repaired[t] = (e, new)
    # a second full verification pass (after all repairs)
    Efinal = set(E_of.values())
    assert len(Efinal) == len(E_of)
    assert not (Efinal & T)
    assert not (Efinal & revT)
    assert not ({pseudo_reverse(e) for e in Efinal} & (T | Efinal) - {e for e in Efinal if pseudo_reverse(e) == e})
    assert all(e[-1] == t[-1] and sorted(e) == sorted(t) for t, e in E_of.items())
    assert not (Efinal & ST), "entrapment equals a shuffled decoy of a target"
    SE = shuffle_decoy_set(Efinal)
    _bad = [(t, e) for t, e in E_of.items()
            if any(hash_shuffle(e, s, 0) in T or hash_shuffle(e, s, 0) in Efinal for s in DECOY_SEEDS)]
    assert not _bad, f"shuffled decoy of an entrapment equals a target/entrapment: {len(_bad)} e.g. {_bad[:3]}"
    unpaired = sorted(T - set(E_of))
    pair_key = {e: t for t, e in E_of.items()}
    meta = {"arm": "P", "estimator": "paired", "r": 1.0, "entrap_seed": seed,
            "generator": "FDRBench v1.1.1 (github.com/Noble-Lab/FDRBench release zip sha256 "
                         "d1740d90decd139c0e2294d766797be800c2b7bd3fd7ae3624e9822ffd7dceed)",
            "fdrbench_cmd": " ".join(os.path.relpath(c, ROOT) if c.startswith(ROOT) else c for c in cmd),
            "fdrbench_wall_s": round(t_fdrbench, 1),
            "fdrbench_pair_file": os.path.relpath(fd_out, ROOT), "fdrbench_pair_file_sha256": sha256(fd_out),
            "target_fasta": os.path.relpath(args.target, ROOT), "target_fasta_sha256": sha256(args.target),
            "n_fdrbench_pairs": len(pairs), "digest_only_ours": len(only_ours), "digest_only_fdrbench": len(only_fd),
            "digest_only_ours_examples": sorted(only_ours)[:10], "digest_only_fdrbench_examples": sorted(only_fd)[:10],
            "n_target_nonstandard_dropped": len(nonstd),
            "n_repaired": len(repaired), "repair_reasons": dict(reasons),
            "repaired_examples": [[t, a, b] for t, (a, b) in list(repaired.items())[:10]],
            "n_unpaired_targets": len(unpaired), "unpaired_targets_examples": unpaired[:20],
            "n_target_eff": len(effective_space(T)), "n_entrapment_eff": len(effective_space(Efinal))}
    write_outputs(args.out, name, sorted(T), sorted(Efinal), pair_key, meta)


def arm_shuffled(args, seed):
    """Decoy-variability arm (DESIGN 6.2): T + C-term-fixed hash(peptide, seed, k) shuffled decoys, no entrapment.
    A shuffle that equals its source, a forward peptide or an earlier decoy is re-drawn (k = 0, 1, ...). Searched with
    external decoys (Sage generate_decoys=false, Comet decoy_search=0), prefix rev_."""
    peps, nonstd = target_peptides(args.target)
    T = sorted(peps)
    fwd = set(T)
    used, dec, n_redraw = set(), {}, 0
    for p in T:
        for k in range(200):
            d = hash_shuffle(p, seed, k)
            if d != p and d not in fwd and d not in used:
                break
            n_redraw += 1
        else:
            continue
        dec[p] = d
        used.add(d)
    assert not (used & fwd) and len(used) == len(dec)
    name = f"{args.name}.S{seed}"
    fa, tsv = os.path.join(args.out, f"{name}.fasta"), os.path.join(args.out, f"{name}.peptides.tsv")
    write_peptide_db(fa, [(f"T_{i:07d}", t) for i, t in enumerate(T)] +
                     [(f"rev_T_{i:07d}", dec[p]) for i, p in enumerate(sorted(dec))])
    with open(tsv, "w") as fh:
        fh.write("seq\tkind\tpair_key\tdecoy_of\n")
        for t in T:
            fh.write(f"{t}\tT\t\t\n")
        for p in sorted(dec):
            fh.write(f"{dec[p]}\tDT\t\t{p}\n")
    meta = {"arm": "S", "decoy_seed": seed, "name": name, "n_target": len(T), "n_decoy": len(dec),
            "n_forward_without_decoy": len(T) - len(dec), "n_redraws": n_redraw,
            "target_fasta": os.path.relpath(args.target, ROOT), "target_fasta_sha256": sha256(args.target),
            "decoy_recipe": f"C-term-fixed hash(peptide,{seed},k) shuffle, k re-drawn on collision",
            "fasta": os.path.relpath(fa, ROOT), "fasta_sha256": sha256(fa),
            "peptides_tsv": os.path.relpath(tsv, ROOT), "peptides_tsv_sha256": sha256(tsv)}
    with open(os.path.join(args.out, f"{name}.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    print(json.dumps({k: v for k, v in meta.items() if not isinstance(v, (dict, list))}), flush=True)


def arm_foreign(args):
    peps, nonstd = target_peptides(args.target)
    T = set(peps)
    F_all = set()
    for f in args.foreign:
        fp, _ = target_peptides(f)
        F_all |= set(fp)
    n0 = len(F_all)
    removed = Counter()
    F = F_all - T
    removed["IL_identical_to_target"] = n0 - len(F)
    tiso = {isobaric_key(t) for t in T}
    F2 = {f for f in F if isobaric_key(f) not in tiso}
    removed["isobaric_NGG_QAG"] = len(F) - len(F2)
    F = F2
    # single-substitution (Hamming-1) neighbours, same length, I/L collapsed: masked-position keys
    tkeys = set()
    for t in T:
        for i in range(len(t)):
            tkeys.add(hash((t[:i], t[i + 1:], len(t))))
    F2 = set()
    for f in F:
        if any(hash((f[:i], f[i + 1:], len(f))) in tkeys for i in range(len(f))):
            continue
        F2.add(f)
    removed["single_substitution_neighbour"] = len(F) - len(F2)
    del tkeys
    F = F2
    revT = {pseudo_reverse(t) for t in T}
    F2 = {f for f in F if f not in revT}
    removed["equal_to_target_decoy"] = len(F) - len(F2)
    F = F2
    ST = shuffle_decoy_set(T)
    F2 = {f for f in F if f not in ST}
    removed["equal_to_shuffled_target_decoy_1001_1005"] = len(F) - len(F2)
    F = F2
    bad_f = shuffle_hits(F, T | F)        # a shuffled decoy of f equals a target or a foreign forward peptide
    F = F - bad_f
    removed["shuffled_decoy_collides_with_target_or_foreign"] = len(bad_f)
    assert not (shuffle_decoy_set(F) & (T | F))
    rinfo = r_by_decile(T, F)
    meta = {"arm": "C", "estimator": "combined", "target_fasta": os.path.relpath(args.target, ROOT),
            "target_fasta_sha256": sha256(args.target),
            "foreign_fastas": [os.path.relpath(f, ROOT) for f in args.foreign],
            "foreign_fasta_sha256": [sha256(f) for f in args.foreign],
            "n_foreign_peptides_raw": n0, "removed": dict(removed), **rinfo,
            "n_target_nonstandard_dropped": len(nonstd)}
    write_outputs(args.out, f"{args.name}.C", sorted(T), sorted(F), None, meta, kind_e="E")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True, help="db family name, e.g. human_contam")
    ap.add_argument("--target", required=True)
    ap.add_argument("--foreign", nargs="*", default=[])
    ap.add_argument("--seeds", nargs="*", type=int, default=[])
    ap.add_argument("--arms", nargs="+", default=["target", "P", "C"])
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    a.target = os.path.abspath(a.target)
    a.out = os.path.abspath(a.out)
    a.foreign = [os.path.abspath(f) for f in a.foreign]
    os.makedirs(a.out, exist_ok=True)
    if "target" in a.arms:
        arm_target(a)
    if "P" in a.arms:
        for s in a.seeds:
            arm_paired(a, s)
    if "S" in a.arms:
        for s in DECOY_SEEDS:
            arm_shuffled(a, s)
    if "C" in a.arms and a.foreign:
        arm_foreign(a)


if __name__ == "__main__":
    main()
