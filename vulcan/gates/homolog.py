"""A1.3 homolog filter for wrong-species G2 nulls (DESIGN §6.0 A1.3).

A null spectrum is removed if its matched peptide -- target or decoy, and (label-symmetric) the other member of the
matched target/decoy pair as well -- is within Hamming distance <= 2 of, or shares a 6-mer with, the null source's
proteome or the contaminants (I/L collapsed).

"Within Hamming distance <= 2 of the proteome" = equal, up to <= 2 substitutions, to some length-L window of any
protein (not only tryptic peptides; LiP / semi-specific peptides are covered).  Implementation:
  * 6-mer: set of every 6-residue window of every protein (I->L).
  * Hamming <= 2 for L <= 19: for each (L, i<j), a 64-bit polynomial hash of every proteome window of length L with
    positions i and j masked; a peptide is near iff for some (i, j) its masked hash is in that set (masking two
    positions covers distances 0, 1 and 2).  For L >= 20 a Hamming <= 2 match leaves <= 3 exact runs of total length
    >= L - 2 >= 18, so one run has length >= 6: the 6-mer test already implies it (pigeonhole), no extra work.
  * Windows that cross a protein boundary or contain a non-standard residue are not used.
Hash collisions (2^-64 per comparison) can only remove more spectra (conservative for the gate).
"""
from __future__ import annotations

import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"
_CODE = np.zeros(256, np.uint8)            # 0 = separator / non-standard
for _k, _a in enumerate(AA):
    _CODE[ord(_a)] = _k + 1
_CODE[ord("I")] = _CODE[ord("L")]          # I/L collapsed
_P = np.uint64(0x9E3779B97F4A7C15)
HAM_MAX_LEN = 19                           # >= 20: implied by the 6-mer test


def read_fasta_seqs(paths) -> list[str]:
    out = []
    for p in paths:
        cur = []
        with open(p) as fh:
            for line in fh:
                if line.startswith(">"):
                    if cur:
                        out.append("".join(cur))
                    cur = []
                else:
                    cur.append(line.strip())
        if cur:
            out.append("".join(cur))
    return out


def encode(seqs) -> np.ndarray:
    """Concatenate with a 0 separator between sequences (and at both ends)."""
    s = "\0" + "\0".join(seqs) + "\0"
    return _CODE[np.frombuffer(s.encode("latin-1"), np.uint8)]


def _pows(L):
    p = np.ones(L, np.uint64)
    with np.errstate(over="ignore"):
        for k in range(1, L):
            p[k] = p[k - 1] * _P
    return p


def _windows(codes: np.ndarray, L: int) -> np.ndarray:
    """(n_valid, L) uint64 matrix of every length-L window without a separator."""
    if codes.size < L:
        return np.zeros((0, L), np.uint64)
    W = np.lib.stride_tricks.sliding_window_view(codes, L)
    ok = (W != 0).all(axis=1)
    return W[ok].astype(np.uint64)


def _hash_rows(M: np.ndarray, pw: np.ndarray) -> np.ndarray:
    with np.errstate(over="ignore"):
        return (M * pw[None, :]).sum(axis=1, dtype=np.uint64)


class HomologFilter:
    """Near-homology test of peptides against one proteome (+ contaminants)."""

    def __init__(self, protein_seqs, k: int = 6, ham: int = 2, max_len_ham: int = HAM_MAX_LEN):
        if ham != 2:
            raise ValueError("only Hamming <= 2 is implemented (A1.3)")
        self.k = k
        self.codes = encode(protein_seqs)
        self.n_residues = int((self.codes != 0).sum())
        kw = _windows(self.codes, k)
        self.kmers = np.unique(_hash_rows(kw, _pows(k)))
        self.max_len_ham = max_len_ham

    @classmethod
    def from_fasta(cls, paths, extra_seqs=(), **kw):
        return cls(read_fasta_seqs(paths) + list(extra_seqs), **kw)

    def _ham_hits(self, L: int, Mq: np.ndarray) -> np.ndarray:
        """bool per query row (n, L uint64 codes): some proteome window of length L within Hamming <= 2."""
        W = _windows(self.codes, L)
        pw = _pows(L)
        hit = np.zeros(Mq.shape[0], bool)
        if W.shape[0] == 0:
            return hit
        with np.errstate(over="ignore"):
            Cw = W * pw[None, :]
            Hw = Cw.sum(axis=1, dtype=np.uint64)
            Cq = Mq * pw[None, :]
            Hq = Cq.sum(axis=1, dtype=np.uint64)
        for i in range(L):
            for j in range(i + 1, L):
                with np.errstate(over="ignore"):
                    q = Hq - Cq[:, i] - Cq[:, j]
                    w = Hw - Cw[:, i] - Cw[:, j]
                qs = np.unique(q)
                pos = np.minimum(np.searchsorted(qs, w), qs.size - 1)
                found = np.unique(w[qs[pos] == w])
                if found.size:
                    hit |= np.isin(q, found)
        return hit

    def near(self, peptides) -> np.ndarray:
        """bool per peptide: shares a k-mer with, or is within Hamming <= 2 of, the proteome (I/L collapsed)."""
        peps = list(peptides)
        out = np.zeros(len(peps), bool)
        if not peps:
            return out
        lens = np.fromiter((len(p) for p in peps), np.int64, len(peps))
        # k-mer test
        pk = _pows(self.k)
        for L in np.unique(lens):
            idx = np.flatnonzero(lens == L)
            if L < self.k:
                continue
            M = _CODE[np.frombuffer("".join(peps[i] for i in idx).encode("latin-1"), np.uint8)].reshape(idx.size, L)
            Wk = np.lib.stride_tricks.sliding_window_view(M, self.k, axis=1).astype(np.uint64)  # (n, L-k+1, k)
            with np.errstate(over="ignore"):
                hk = (Wk * pk[None, None, :]).sum(axis=2, dtype=np.uint64)
            out[idx] |= np.isin(hk, self.kmers).any(axis=1)
            # Hamming <= 2 (only L <= max_len_ham; longer ones are implied by the k-mer test)
            if L <= self.max_len_ham:
                todo = idx[~out[idx]]
                if todo.size == 0:
                    continue
                Mt = M[np.isin(idx, todo)].astype(np.uint64)
                out[todo] |= self._ham_hits(int(L), Mt)
        return out


def hamming_le2_bruteforce(pep: str, proteins) -> bool:
    """Reference implementation (tests only)."""
    p = pep.replace("I", "L")
    L = len(p)
    for s in proteins:
        s = s.replace("I", "L")
        for a in range(len(s) - L + 1):
            w = s[a:a + L]
            if sum(x != y for x, y in zip(p, w)) <= 2:
                return True
    return False


def shares_kmer_bruteforce(pep: str, proteins, k: int = 6) -> bool:
    p = pep.replace("I", "L")
    km = {p[i:i + k] for i in range(len(p) - k + 1)}
    for s in proteins:
        s = s.replace("I", "L")
        if any(s[i:i + k] in km for i in range(len(s) - k + 1)):
            return True
    return False
