"""Amendment A2.2(b) identity screen against the gate database T∪D, label-symmetric and outcome-independent.

The gate DB is handled ONLY as an unordered set of peptidoforms (``GateDB``): the label-free index side of an L1 index
(peptidoforms.arrow + residues.npy; labels.arrow is never opened) or, for tests and the P1 swap test, the union of a
"T" and a "D" peptidoform list, which is immediately canonicalised (sorted by (mass, sequence, mask), duplicates
removed) so the result cannot depend on which list carried which label or on input order.

Fixed implementation choices (A2.11 item 5 "f(c) details"; stricter reading where the text is open):
  * W = production precursor windows of the CERTIFIED precursor: neutral mass M_k = (mz_mono - PROTON) * z - k * NEUTRON
    for every k in the production isotope-error set (0, 1, 2), +-10 ppm of M_k (Sage f32 convention, as
    vulcan.gpu.host.make_hypotheses); candidates = DB peptidoforms with mass in any window (inclusive bounds).
  * Fragment masses: Sage residue masses (vulcan.gpu.engine.SAGE_MONO), static C +57.021464, variable M +15.9949 on the
    positions set in the peptidoform mask (bit k = residue k; bit 63 = variable N-term, value from the index meta),
    computed in float64; b_i = sum(res[:i]) + PROTON, y_i = sum(res[L-i:]) + H2O + PROTON (z = 1).
  * Ions of p: b2..b(L-1) and y2..y(L-1) (2(L-2) ions, as the text).  Ions of c: ALL b1..b(Lc-1) and y1..y(Lc-1)
    (stricter: more coinciding ions).  A p ion coincides if ANY c ion (b or y, cross-type allowed: stricter) lies
    within |m_p - m_c| <= 20e-6 * m_p.  No m/z range restriction (stricter than K0's 150-2000 window, which could only
    drop coinciding ions).
  * f(c) = #coinciding p ions / (2(L-2)); remove iff max_c f(c) >= 0.20.  The text's K2-kernel route is an
    implementation hint; it is NOT used because K0 (min_peaks 15, fragment m/z window, top-150 peaks) would drop ions of
    short peptides and lower f.  This exact numpy computation is the definition.
  * Hamming clause: p is near if some isobaric variant v of p (every subset of the rewrites N->GG, Q->AG, GG->N, AG->Q,
    applied to non-overlapping sites) has Hamming distance <= 2, I/L collapsed, to a T∪D PEPTIDE (stripped sequence) of
    the same length.  "T∪D peptide" = the DB's peptide entries (not arbitrary protein windows), as the text says.
    More than VARIANT_CAP = 4096 variants -> near (stricter).
"""
from __future__ import annotations

import os
from dataclasses import dataclass
import numpy as np

PROTON = 1.0072764
NEUTRON = 1.00335
H2O = 18.010565
# Sage MONOISOTOPIC masses by letter (vulcan.gpu.engine.SAGE_MONO; copied to stay torch-free; a test pins equality)
SAGE_MONO = np.array([71.03711, 0.0, 103.00919, 115.02694, 129.04259, 147.0684, 57.02146, 137.05891, 113.08406,
                      0.0, 128.09496, 113.08406, 131.0405, 114.04293, 237.14774, 97.05276, 128.05858, 156.1011,
                      87.03203, 101.04768, 150.95363, 99.06841, 186.07932, 0.0, 163.06332, 0.0], np.float32
                     ).astype(np.float64)
STATIC = {"C": 57.021464}
VAR = {"M": 15.9949}
NTERM_BIT = np.uint64(1) << np.uint64(63)
VARIANT_CAP = 4096
F_THRESHOLD = 0.20
FRAG_PPM = 20.0
PREC_PPM = 10.0
ISOTOPE_ERRORS = (0, 1, 2)


def residue_masses(seq: str, mask: int = 0, static=STATIC, var=VAR, nterm: float = 0.0) -> np.ndarray:
    m = np.empty(len(seq))
    for k, r in enumerate(seq):
        base = SAGE_MONO[ord(r) - 65]
        if (int(mask) >> k) & 1:
            m[k] = base + var.get(r, 0.0)
        else:
            m[k] = base + static.get(r, 0.0)
    if int(mask) & (1 << 63) and len(seq):
        m[0] += nterm
    return m


def peptide_mass(seq: str, mask: int = 0, **kw) -> float:
    return float(residue_masses(seq, mask, **kw).sum() + H2O)


def ion_ladder(seq: str, mask: int = 0, first: int = 1, **kw) -> np.ndarray:
    """z=1 b and y ion m/z: b_first..b_(L-1), y_first..y_(L-1)."""
    r = residue_masses(seq, mask, **kw)
    L = r.size
    cs = np.cumsum(r)
    b = cs[first - 1:L - 1] + PROTON
    rc = np.cumsum(r[::-1])
    y = rc[first - 1:L - 1] + H2O + PROTON
    return np.concatenate([b, y])


def f_overlap(p_ions: np.ndarray, c_ions: np.ndarray, ppm: float = FRAG_PPM) -> float:
    """Fraction of p's ions with a c ion within +-ppm (relative to the p ion)."""
    if p_ions.size == 0:
        return 0.0
    s = np.sort(c_ions)
    tol = p_ions * ppm * 1e-6
    lo = np.searchsorted(s, p_ions - tol, side="left")
    hi = np.searchsorted(s, p_ions + tol, side="right")
    return float(np.count_nonzero(hi > lo)) / p_ions.size


# ------------------------------------------------------------------------------------------------ gate DB
@dataclass
class GateDB:
    """Unordered T∪D peptidoform set, label-free.  Rows sorted by (mass, sequence, mask) for list input; for an L1
    index the index's own content-defined order (Sage f32 mass, xxhash id) is kept -- neither depends on labels.
    mass: float64 neutral mass used for windows (Sage f32 mass32 for an L1 index, as production); sequences live in a
    residue buffer (seq_off, length) so 1e7 peptidoforms fit in memory."""
    mass: np.ndarray
    mask: np.ndarray
    residues: np.ndarray
    seq_off: np.ndarray
    length: np.ndarray
    nterm: float = 0.0

    @classmethod
    def from_lists(cls, *lists, nterm: float = 0.0):
        """lists: any number of iterables of (sequence, mask), e.g. the T list and the D list.  The union is
        canonicalised (set, then sorted) so the result is invariant to which list is which and to order."""
        allp = sorted({(str(s), int(m)) for lst in lists for s, m in lst})
        mass = np.array([peptide_mass(s, m, nterm=nterm) for s, m in allp], float)
        o = sorted(range(len(allp)), key=lambda i: (mass[i], allp[i][0], allp[i][1]))
        allp = [allp[i] for i in o]
        mass = mass[o] if len(o) else mass
        seqs = [s for s, _ in allp]
        length = np.array([len(x) for x in seqs], np.int64)
        off = np.zeros(len(seqs), np.int64)
        if len(seqs):
            off[1:] = np.cumsum(length)[:-1]
        res = np.frombuffer("".join(seqs).encode(), np.uint8).copy()
        return cls(mass, np.array([m for _, m in allp], np.uint64), res, off, length, nterm)

    @classmethod
    def from_l1(cls, index_dir: str):
        """Label-free side of an L1 index ONLY (peptidoforms.arrow, residues.npy, meta.json); never labels.arrow."""
        import json

        import pyarrow as pa
        with pa.OSFile(os.path.join(index_dir, "peptidoforms.arrow"), "rb") as src:
            t = pa.ipc.open_file(src).read_all()
        meta = json.load(open(os.path.join(index_dir, "meta.json")))
        nterm = float(meta.get("params", {}).get("variable_mods", {}).get("^", 0.0))
        return cls(t.column("mass32").to_numpy().astype(np.float64), t.column("mask").to_numpy().astype(np.uint64),
                   np.load(os.path.join(index_dir, "residues.npy")), t.column("seq_off").to_numpy().astype(np.int64),
                   t.column("length").to_numpy().astype(np.int64), nterm)

    def __len__(self):
        return self.mass.size

    def seq(self, r: int) -> str:
        o = int(self.seq_off[r])
        return self.residues[o:o + int(self.length[r])].tobytes().decode()

    def window_rows(self, mz_mono: float, z: int, ppm: float = PREC_PPM, isotopes=ISOTOPE_ERRORS) -> np.ndarray:
        rows = []
        for k in isotopes:
            M = (mz_mono - PROTON) * z - k * NEUTRON
            lo, hi = M - M * ppm * 1e-6, M + M * ppm * 1e-6
            a = np.searchsorted(self.mass, lo, side="left")
            b = np.searchsorted(self.mass, hi, side="right")
            rows.append(np.arange(a, b))
        return np.unique(np.concatenate(rows)) if rows else np.zeros(0, np.int64)

    def peptides(self) -> list[str]:
        """Distinct stripped sequences (sorted)."""
        return sorted({self.seq(r) for r in range(len(self))})

    def fingerprint(self) -> str:
        """sha256 of the canonical content (sorted (sequence, mask) set): equal for T∪D and D∪T."""
        import hashlib
        h = hashlib.sha256()
        for s, m in sorted((self.seq(r), int(self.mask[r])) for r in range(len(self))):
            h.update(f"{s}:{m};".encode())
        return h.hexdigest()


def max_f_in_window(db: GateDB, p_seq: str, p_mask: int, mz_mono: float, z: int, cache: dict | None = None,
                    ppm_frag: float = FRAG_PPM, ppm_prec: float = PREC_PPM, isotopes=ISOTOPE_ERRORS,
                    static=STATIC, var=VAR) -> float:
    """max over c in (T∪D) ∩ W of f(c)  (0.0 if W is empty)."""
    p_ions = ion_ladder(p_seq, p_mask, first=2, static=static, var=var)
    best = 0.0
    for r in db.window_rows(mz_mono, z, ppm_prec, isotopes):
        key = (db.seq(r), int(db.mask[r]))
        if cache is not None and key in cache:
            c_ions = cache[key]
        else:
            c_ions = ion_ladder(key[0], key[1], first=1, static=static, var=var, nterm=db.nterm)
            if cache is not None:
                cache[key] = c_ions
        f = f_overlap(p_ions, c_ions, ppm_frag)
        if f > best:
            best = f
            if best >= 1.0:
                break
    return best


# ------------------------------------------------------------------------------------------------ Hamming (isobaric)
_ILC = str.maketrans("I", "L")
REWRITES = (("N", "GG"), ("Q", "AG"), ("GG", "N"), ("AG", "Q"))


def isobaric_variants(p: str, cap: int = VARIANT_CAP):
    """All strings obtained by applying any set of non-overlapping rewrites N<->GG, Q<->AG to p (I/L collapsed),
    including p itself.  Returns None if more than ``cap`` variants exist (caller treats p as near)."""
    p = p.translate(_ILC)
    sites = []
    for a, b in REWRITES:
        i = p.find(a)
        while i >= 0:
            sites.append((i, i + len(a), b))
            i = p.find(a, i + 1)
    sites.sort()
    out = set()

    def rec(k, chosen):
        if len(out) > cap:
            return
        if k == len(sites):
            s, last = [], 0
            for i, j, b in chosen:
                s += [p[last:i], b]
                last = j
            out.add("".join(s) + p[last:])
            return
        rec(k + 1, chosen)
        i, j, b = sites[k]
        if all(j0 <= i or i0 >= j for i0, j0, _ in chosen):
            rec(k + 1, chosen + [(i, j, b)])

    rec(0, [])
    return None if len(out) > cap else sorted(out)


_HP = np.uint64(0x9E3779B97F4A7C15)


def _block_hash(M: np.ndarray) -> np.ndarray:
    """64-bit polynomial hash of the rows of a uint8 matrix (collisions only add candidates; verified exactly)."""
    h = np.zeros(M.shape[0], np.uint64)
    with np.errstate(over="ignore"):
        for j in range(M.shape[1]):
            h = h * _HP + M[:, j].astype(np.uint64) + np.uint64(1)
    return h


class PeptideHamming:
    """Exact-length Hamming <= 2 (I/L collapsed) of query strings against a fixed peptide set.

    Pigeonhole: positions split into 3 contiguous blocks; distance <= 2 implies one identical block.  Candidates =
    peptides of the same length sharing a block (64-bit hash lookup); every candidate is verified exactly, so hash
    collisions cannot change the result."""

    def __init__(self, peptides):
        peps = sorted({p.translate(_ILC) for p in peptides})
        groups: dict[int, list[str]] = {}
        for p in peps:
            groups.setdefault(len(p), []).append(p)
        self.by_len: dict[int, np.ndarray] = {}
        self.index: dict[tuple[int, int], tuple[np.ndarray, np.ndarray]] = {}
        for L, ps in groups.items():
            M = np.frombuffer("".join(ps).encode(), np.uint8).reshape(len(ps), L)
            self.by_len[L] = M
            for bi, (a, b) in enumerate(self._bounds(L)):
                if b <= a:
                    continue
                h = _block_hash(M[:, a:b])
                o = np.argsort(h, kind="stable")
                self.index[(L, bi)] = (h[o], o)

    @staticmethod
    def _bounds(L):
        q, r = divmod(L, 3)
        out, s = [], 0
        for i in range(3):
            n = q + (1 if i < r else 0)
            out.append((s, s + n))
            s += n
        return out

    def near(self, q: str, maxd: int = 2) -> bool:
        q = q.translate(_ILC)
        L = len(q)
        M = self.by_len.get(L)
        if M is None:
            return False
        if L <= maxd:
            return True
        qa = np.frombuffer(q.encode(), np.uint8)
        cand = []
        for bi, (a, b) in enumerate(self._bounds(L)):
            if (L, bi) not in self.index:
                continue
            hs, o = self.index[(L, bi)]
            hq = _block_hash(qa[None, a:b])[0]
            lo = np.searchsorted(hs, hq, side="left")
            hi = np.searchsorted(hs, hq, side="right")
            if hi > lo:
                cand.append(o[lo:hi])
        if not cand:
            return False
        rows = np.unique(np.concatenate(cand))
        return bool(((M[rows] != qa[None, :]).sum(axis=1) <= maxd).any())


def hamming_isobaric_near(ph: PeptideHamming, p: str, cap: int = VARIANT_CAP) -> bool:
    vs = isobaric_variants(p, cap)
    if vs is None:
        return True
    return any(ph.near(v) for v in vs)


def hamming_isobaric_bruteforce(p: str, peptides, cap: int = VARIANT_CAP) -> bool:
    """Reference implementation (tests only)."""
    vs = isobaric_variants(p, cap)
    if vs is None:
        return True
    P = {x.translate(_ILC) for x in peptides}
    for v in vs:
        for q in P:
            if len(q) == len(v) and sum(a != b for a, b in zip(v, q)) <= 2:
                return True
    return False


def all_pairs_hamming_le(a: str, b: str, d: int = 2) -> bool:
    return len(a) == len(b) and sum(x != y for x, y in zip(a, b)) <= d


__all__ = ["GateDB", "PeptideHamming", "max_f_in_window", "hamming_isobaric_near", "isobaric_variants", "ion_ladder",
           "f_overlap", "peptide_mass"]
