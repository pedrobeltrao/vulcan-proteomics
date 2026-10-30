"""Label-blind input preparation for the split-then-joint encoder: spectrum peak tokens and peptidoform tokens.

Nothing here reads labels.arrow; peptidoforms come from the L1 index side (residues + var-mod mask), spectra from the
store.  Deterministic (stable sorts, no RNG)."""
from __future__ import annotations

import hashlib

import numpy as np

from .model import OX_TOKEN

N_PEAKS = 63
MAX_LEN = 50
MZ_MIN, MZ_MAX = 150.0, 2000.0


def spectra_topk(sp, n_peaks: int = N_PEAKS):
    """Spectra -> (mz (S, P) f32 ascending m/z, intensity (S, P) f32, n (S,) i16).  Top-n_peaks centroids by intensity
    (ties: lower m/z first) within [150, 2000] m/z, then sorted by m/z; zero padded."""
    S = len(sp)
    lens = np.diff(sp.offsets)
    seg = np.repeat(np.arange(S, dtype=np.int64), lens)
    mz = sp.mz.astype(np.float32)
    it = sp.intensity.astype(np.float32)
    ok = np.isfinite(mz) & (mz >= MZ_MIN) & (mz <= MZ_MAX) & (it > 0)
    seg, mz, it = seg[ok], mz[ok], it[ok]
    o = np.lexsort((mz, -it, seg))                     # by spectrum, intensity desc, m/z asc
    seg, mz, it = seg[o], mz[o], it[o]
    first = np.searchsorted(seg, seg, side="left")
    rank = np.arange(len(seg)) - first
    k = rank < n_peaks
    seg, mz, it, rank = seg[k], mz[k], it[k], rank[k]
    o = np.lexsort((mz, seg))
    seg, mz, it = seg[o], mz[o], it[o]
    first = np.searchsorted(seg, seg, side="left")
    pos = np.arange(len(seg)) - first
    M = np.zeros((S, n_peaks), np.float32)
    I = np.zeros((S, n_peaks), np.float32)
    M[seg, pos] = mz
    I[seg, pos] = it
    n = np.bincount(seg, minlength=S).astype(np.int16)
    return M, I, n


def spectrum_hash(M: np.ndarray, I: np.ndarray, pmz: np.ndarray) -> np.ndarray:
    """Exact-duplicate hash per spectrum: blake2b-64 of (precursor m/z rounded to 1e-3, top-peak m/z rounded to 1e-3,
    relative intensity rounded to 1e-3).  -> uint64 (S,)."""
    out = np.zeros(len(M), np.uint64)
    mzr = np.round(M.astype(np.float64) * 1000).astype(np.int64)
    rel = I / np.maximum(I.max(1, keepdims=True), 1e-12)
    itr = np.round(rel * 1000).astype(np.int32)
    pr = np.round(np.nan_to_num(pmz.astype(np.float64)) * 1000).astype(np.int64)
    for i in range(len(M)):
        h = hashlib.blake2b(digest_size=8)
        h.update(pr[i].tobytes()); h.update(mzr[i].tobytes()); h.update(itr[i].tobytes())
        out[i] = np.frombuffer(h.digest(), np.uint64)[0]
    return out


def peptide_tokens(l1, rows: np.ndarray, max_len: int = MAX_LEN) -> np.ndarray:
    """L1 rows -> (N, max_len) uint8 tokens: residue letter -> ord - 64; oxidised M (var-mod bit set on an M) -> OX."""
    rows = np.asarray(rows, np.int64)
    N = len(rows)
    off = l1.seq_off[rows].astype(np.int64)
    ln = l1.length[rows].astype(np.int64)
    mask = l1.mask[rows].astype(np.uint64)
    pos = np.arange(max_len)[None, :]
    valid = pos < ln[:, None]
    idx = np.where(valid, off[:, None] + pos, 0)
    res = l1.residues[idx]
    tok = np.where(valid, res.astype(np.int16) - 64, 0)
    bits = ((mask[:, None] >> pos.astype(np.uint64)) & np.uint64(1)).astype(bool) & valid
    tok = np.where(bits & (res == ord("M")), OX_TOKEN, tok)
    return tok.astype(np.uint8)


def stripped_il(l1, rows):
    """I/L-collapsed stripped sequences (bytes) of L1 rows (for negative filtering; training side only)."""
    out = []
    for r in np.asarray(rows, np.int64):
        o = int(l1.seq_off[r]); n = int(l1.length[r])
        out.append(l1.residues[o:o + n].tobytes().replace(b"I", b"L"))
    return out


def spectra_topk_gpu(sp, n_peaks: int = N_PEAKS, device="cuda"):
    """GPU version of spectra_topk with identical selection and ordering (stable sorts mirror the CPU lexsort);
    returns torch tensors (mz, intensity, n) on `device`.  Verified bit-identical to spectra_topk."""
    import torch
    S = len(sp)
    lens = torch.from_numpy(np.diff(sp.offsets)).to(device)
    seg = torch.repeat_interleave(torch.arange(S, device=device), lens)
    mz = torch.from_numpy(sp.mz.astype(np.float32)).to(device)
    it = torch.from_numpy(sp.intensity.astype(np.float32)).to(device)
    ok = torch.isfinite(mz) & (mz >= MZ_MIN) & (mz <= MZ_MAX) & (it > 0)
    seg, mz, it = seg[ok], mz[ok], it[ok]
    o = torch.argsort(mz, stable=True)
    o = o[torch.argsort(-it[o], stable=True)]
    o = o[torch.argsort(seg[o], stable=True)]
    seg, mz, it = seg[o], mz[o], it[o]
    first = torch.searchsorted(seg, seg, side="left")
    k = (torch.arange(len(seg), device=device) - first) < n_peaks
    seg, mz, it = seg[k], mz[k], it[k]
    o = torch.argsort(mz, stable=True)
    o = o[torch.argsort(seg[o], stable=True)]
    seg, mz, it = seg[o], mz[o], it[o]
    first = torch.searchsorted(seg, seg, side="left")
    pos = torch.arange(len(seg), device=device) - first
    M = torch.zeros((S, n_peaks), dtype=torch.float32, device=device)
    I = torch.zeros((S, n_peaks), dtype=torch.float32, device=device)
    M[seg, pos] = mz
    I[seg, pos] = it
    n = torch.bincount(seg, minlength=S)
    return M, I, n
