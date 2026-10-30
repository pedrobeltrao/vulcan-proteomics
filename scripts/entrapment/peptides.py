"""Shared peptide utilities for the M1e entrapment databases (DESIGN.md §4 L1, §6.1, §6.3).

* tryptic digest identical to FDRBench `-enzyme 1` / Sage `cleave_at=KR, restrict=P` / Comet enzyme 1
  (cleave after K/R unless followed by P; <= 2 missed cleavages; length 7-50; no N-term Met clipping)
* I/L collapse (I -> L), pseudo-reverse decoy with the protease C-term residue fixed,
* hash(peptide, seed) shuffles (C-term fixed) that do not depend on database order.
"""
from __future__ import annotations

import hashlib
import re

import numpy as np

STD_AA = set("ACDEFGHIKLMNPQRSTVWY")
# monoisotopic residue masses (Cys carbamidomethylated is applied at search time; here unmodified,
# matching Sage's peptide_min/max_mass filter which is applied to the unmodified+static-mod mass)
MONO = {
    "G": 57.02146372, "A": 71.03711379, "S": 87.03202841, "P": 97.05276385, "V": 99.06841391,
    "T": 101.0476785, "C": 103.0091845, "L": 113.0840640, "I": 113.0840640, "N": 114.0429275,
    "D": 115.0269431, "Q": 128.0585775, "K": 128.0949630, "E": 129.0425931, "M": 131.0404846,
    "H": 137.0589119, "F": 147.0684139, "R": 156.1011110, "Y": 163.0633285, "W": 186.0793130,
}
H2O = 18.0105646863
CAM = 57.021464

MIN_LEN, MAX_LEN, MISSED = 7, 50, 2
MIN_MASS, MAX_MASS = 500.0, 5000.0


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            elif line:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def accession(header: str) -> str:
    tok = header.split()[0]
    parts = tok.split("|")
    return parts[1] if len(parts) >= 3 else tok


_SITE = re.compile(r"(?<=[KR])(?!P)")


def digest(seq: str, missed: int = MISSED, min_len: int = MIN_LEN, max_len: int = MAX_LEN):
    """Fully tryptic peptides (K/R not before P) with up to `missed` missed cleavages."""
    frags = [f for f in _SITE.split(seq) if f]
    out = []
    n = len(frags)
    for i in range(n):
        pep = ""
        for j in range(i, min(n, i + missed + 1)):
            pep += frags[j]
            if len(pep) > max_len:
                break
            if len(pep) >= min_len:
                out.append(pep)
    return out


def il(s: str) -> str:
    return s.replace("I", "L")


def mass(pep: str, cam: bool = True) -> float:
    m = H2O + sum(MONO[a] for a in pep)
    if cam:
        m += CAM * pep.count("C")
    return m


def in_space(pep: str) -> bool:
    """Effective search space filter: standard residues, length and (CAM) mass range."""
    return (MIN_LEN <= len(pep) <= MAX_LEN and set(pep) <= STD_AA and MIN_MASS <= mass(pep) <= MAX_MASS)


def pseudo_reverse(pep: str) -> str:
    """DESIGN L1 default decoy: reverse with the protease C-term residue fixed."""
    return pep[-2::-1] + pep[-1]


def _rng(*keys) -> np.random.Generator:
    h = hashlib.blake2b("|".join(map(str, keys)).encode(), digest_size=8).digest()
    return np.random.default_rng(int.from_bytes(h, "little"))


def hash_shuffle(pep: str, seed: int, attempt: int = 0) -> str:
    """C-term-fixed shuffle seeded by hash(peptide, seed, attempt) (order independent)."""
    body = np.array(list(pep[:-1]))
    _rng(pep, seed, attempt).shuffle(body)
    return "".join(body) + pep[-1]


def isobaric_key(pep: str) -> str:
    """Canonical form for isobaric neighbours (DESIGN 6.3: N<->GG, Q<->AG): I->L, then N->GG and Q->AG/GA expanded,
    and every maximal run of G/A letters sorted. Two peptides share a key iff they differ only by N/GG and
    Q/AG/GA exchanges (plus G/A order inside a run, which is slightly broader than the spec and conservative: it only
    ever removes more foreign peptides). The previous greedy GG->N, AG->Q replacement was order dependent and missed
    e.g. AGG vs QG."""
    s = il(pep).replace("N", "GG").replace("Q", "AG")
    return _GA_RUN.sub(lambda m: "".join(sorted(m.group(0))), s)


_GA_RUN = re.compile(r"[GA]{2,}")

DECOY_SEEDS = (1001, 1002, 1003, 1004, 1005)   # DESIGN 6.1


def _shuf_chunk(args):
    seqs, seeds = args
    return [hash_shuffle(s, sd, 0) for s in seqs for sd in seeds]


def shuffle_decoy_set(seqs, seeds=DECOY_SEEDS, procs=None):
    """Set of C-term-fixed hash(peptide, seed) shuffled decoys (attempt 0) of every peptide for every decoy seed
    (the decoys of the 5-seed variability arm, DESIGN 6.2). Multiprocessed; deterministic."""
    import os
    from multiprocessing import Pool
    seqs = sorted(seqs)
    procs = procs or len(os.sched_getaffinity(0))
    chunks = [(seqs[i:i + 20000], seeds) for i in range(0, len(seqs), 20000)]
    out = set()
    with Pool(procs) as pool:
        for r in pool.imap_unordered(_shuf_chunk, chunks):
            out.update(r)
    return out


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def write_peptide_db(path: str, entries):
    """entries: iterable of (accession, sequence). One FASTA entry per peptide."""
    with open(path, "w") as fh:
        for acc, s in entries:
            fh.write(f">{acc}\n{s}\n")


def build_decoys(forward: dict[str, str], seed: int = 1001, max_tries: int = 50):
    """Pseudo-reverse decoys for every forward peptide (targets and entrapments alike).

    forward: {sequence (I/L collapsed): class}. A decoy that is palindromic or collides with any forward
    peptide or an already-made decoy is replaced by hash(peptide, DECOY_SEED) shuffles.
    Returns {forward seq: decoy seq} and the number of replacements.
    """
    fwd = set(forward)
    used = set()
    out = {}
    n_repl = 0
    for p in sorted(forward):  # sorted: deterministic, independent of DB order
        d = pseudo_reverse(p)
        k = 0
        while d == p or d in fwd or d in used:
            if k >= max_tries:
                d = None
                break
            d = hash_shuffle(p, seed, k)
            k += 1
        if k:
            n_repl += 1
        if d is not None:
            out[p] = d
            used.add(d)
    return out, n_repl


_LOOKUP = frozenset()


def _hit_chunk(args):
    seqs, seeds = args
    return [s for s in seqs if any(hash_shuffle(s, sd, 0) in _LOOKUP for sd in seeds)]


def shuffle_hits(seqs, lookup, seeds=DECOY_SEEDS, procs=None):
    """Peptides in `seqs` whose shuffled decoy (any decoy seed) is in `lookup` (fork-shared, multiprocessed)."""
    import os
    from multiprocessing import get_context
    global _LOOKUP
    _LOOKUP = lookup
    seqs = sorted(seqs)
    procs = procs or len(os.sched_getaffinity(0))
    chunks = [(seqs[i:i + 20000], seeds) for i in range(0, len(seqs), 20000)]
    out = []
    with get_context("fork").Pool(procs) as pool:
        for r in pool.imap_unordered(_hit_chunk, chunks):
            out.extend(r)
    return set(out)
