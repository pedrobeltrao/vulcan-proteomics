"""Semi-specific tryptic digestion with specificity classes (DESIGN §4 L1 'semi-specific', L2 mode 'semi-specific').

Enumeration follows Sage's semi-enzymatic digest (sage-core enzyme.rs, v0.14.7): specific sites (trypsin K/R, not
before P; protein termini are sites) with up to `missed` missed cleavages give the specific peptides; every prefix
and suffix of a specific peptide is a semi-specific candidate; the final length filter [min_len, max_len] is applied
after truncation (so half-tryptic peptides of long tryptic parents are kept, as in Sage).  Sequences are I/L collapsed.

Specificity class of one occurrence (a label-free property of the sequence and its protein context):
    0  specific (both ends at sites)
    1  N-semi: N-terminus not at a site, C-terminus at a site (PK cut on the N-terminal side)
    2  C-semi: C-terminus not at a site, N-terminus at a site
A sequence occurring several times gets the most specific class (min).  The decoy of a pair inherits the class of its
forward member, so decoy class counts equal target class counts exactly (DESIGN L1: "semi-specific decoys come from
the decoy proteins, so class counts match"; see vulcan/lip/index.py for why pairing replaces protein-level decoy
digestion here).  Biological exclusions (protein N-terminus, Met removal, signal/transit/propeptides) are NOT
applied at search time: they are a post-identification call (vulcan.lip.halftryptic), so the search space and the
FDR strata do not depend on annotations.
"""
from __future__ import annotations

import multiprocessing as mp
import os

import numpy as np
import pandas as pd

STD = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta(paths) -> tuple[list[str], list[str]]:
    acc, seqs, cur = [], [], []
    for p in paths:
        with open(p) as fh:
            for line in fh:
                if line.startswith(">"):
                    if cur:
                        seqs.append("".join(cur))
                        cur = []
                    acc.append(line[1:].split()[0])
                else:
                    cur.append(line.strip())
            if cur:
                seqs.append("".join(cur))
                cur = []
    assert len(acc) == len(seqs)
    return acc, seqs


def sites(seq: str, cleave="KR", restrict="P") -> list[int]:
    """Cleavage positions (peptide boundaries) incl. 0 and len(seq)."""
    n = len(seq)
    s = [0]
    for i, c in enumerate(seq):
        if c in cleave and i + 1 < n and seq[i + 1] not in restrict:
            s.append(i + 1)
    if s[-1] != n:
        s.append(n)
    return s


def _digest_chunk(args):
    seqs, base, missed, min_len, max_len, semi = args
    out_seq, out_cls, out_mc, out_prot = [], [], [], []
    for k, prot in enumerate(seqs):
        prot = prot.replace("I", "L")
        st = sites(prot)
        sset = set(st)
        n = len(prot)
        seen = {}
        for i in range(len(st) - 1):
            for j in range(i + 1, min(i + 2 + missed, len(st))):
                a, b = st[i], st[j]
                cand = [(a, b)]
                if semi:
                    cand.extend((a, e) for e in range(a + min_len, b))
                    cand.extend((s, b) for s in range(a + 1, b - min_len + 1))
                for x, y in cand:
                    L = y - x
                    if L < min_len or L > max_len:
                        continue
                    pep = prot[x:y]
                    cls = 0 if (x in sset and y in sset) else (1 if y in sset else 2)
                    old = seen.get(pep)
                    if old is None or cls < old[0]:
                        mc = sum(1 for q in range(x + 1, y) if q in sset)
                        seen[pep] = (cls, mc)
        pid = base + k
        for pep, (cls, mc) in seen.items():
            if not set(pep) <= STD:
                continue
            out_seq.append(pep)
            out_cls.append(cls)
            out_mc.append(mc)
            out_prot.append(pid)
    return out_seq, np.array(out_cls, np.uint8), np.array(out_mc, np.uint8), np.array(out_prot, np.uint32)


def digest(fasta_paths, missed=2, min_len=7, max_len=50, semi=True, workers=None) -> dict:
    """Returns the vulcan_core.digest_fasta dict layout (keys, missed, prot_off, prot_idx, accessions, n_nonstandard)
    plus 'cls' (u8 per key)."""
    acc, seqs = read_fasta(fasta_paths)
    workers = workers or max(1, min(32, len(os.sched_getaffinity(0))))
    step = max(1, (len(seqs) + workers * 8 - 1) // (workers * 8))
    jobs = [(seqs[i:i + step], i, missed, min_len, max_len, semi) for i in range(0, len(seqs), step)]
    with mp.get_context("fork").Pool(workers) as pool:
        parts = pool.map(_digest_chunk, jobs)
    seq = [s for p in parts for s in p[0]]
    cls = np.concatenate([p[1] for p in parts])
    mc = np.concatenate([p[2] for p in parts])
    prot = np.concatenate([p[3] for p in parts])
    del parts
    codes, uniq = pd.factorize(pd.Series(seq), sort=True)
    del seq
    nk = len(uniq)
    kcls = np.full(nk, 255, np.uint8)
    np.minimum.at(kcls, codes, cls)
    # missed cleavages of the most specific occurrence (ties: min)
    best = cls == kcls[codes]
    kmc = np.full(nk, 255, np.uint8)
    np.minimum.at(kmc, codes[best], mc[best])
    o = np.lexsort((prot, codes))
    c_s, p_s = codes[o], prot[o]
    keep = np.r_[True, (c_s[1:] != c_s[:-1]) | (p_s[1:] != p_s[:-1])]
    c_s, p_s = c_s[keep], p_s[keep]
    off = np.zeros(nk + 1, np.int64)
    np.add.at(off, c_s + 1, 1)
    off = np.cumsum(off)
    return {"keys": list(uniq), "missed": kmc, "cls": kcls, "prot_off": off, "prot_idx": p_s.astype(np.uint32),
            "accessions": acc, "n_nonstandard": 0, "sequences": seqs}
