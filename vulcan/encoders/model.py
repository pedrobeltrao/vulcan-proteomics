"""Split-then-joint spectrum/peptide encoder (DESIGN §4 L3 cross-encoder arm; §10 E2b/E5/E6).

Split stage
  * spectrum encoder: top-`n_peaks` centroids (by intensity, m/z in [150, 2000]) -> tokens [Fourier(m/z), sqrt rel.
    intensity, intensity rank] + one charge token -> `n_spec` pre-LN Transformer layers.
  * peptide encoder: residue tokens [residue/mod embedding + Fourier(b1+ m/z of the prefix ending here, y1+ m/z of the
    suffix starting after it) + position] + charge token -> `n_pep` layers.
  * pooled outputs (masked mean -> linear -> L2 norm, `emb`-d) ARE the E5 bi-encoder embeddings.
Joint stage (cross-encoder head)
  * `n_joint` blocks: peptide token states cross-attend to the spectrum token states (+ FFN). The attention logits get
    an additive physics bias: a learned per-head mix of Gaussian fragment-match kernels (b1, y1, b2, y2 theoretical
    m/z of the token vs each peak m/z, sigma = `tol_ppm`). Pooled (mean, max) + bi-encoder cosine -> MLP -> logit.

Label blindness: inputs are spectrum peaks/charge and the peptidoform (residues + var-mod mask) only. No target/decoy
or entrapment flag, no precursor mass (candidates are already precursor-matched; withholding it keeps shifted-
precursor negatives a fragment-level task and removes a trivial cue).  The model is a function of (spectrum,
peptidoform) and is identical for a target and its paired decoy (same code path), as required by G0.
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

PROTON = 1.00727646688
H2O = 18.0105646837
# residue masses (monoisotopic) incl. fixed carbamidomethyl C; index by ord(c) - 65
_AA = {"A": 71.03711, "C": 103.00919 + 57.021464, "D": 115.02694, "E": 129.04259, "F": 147.06841, "G": 57.02146,
       "H": 137.05891, "I": 113.08406, "K": 128.09496, "L": 113.08406, "M": 131.04049, "N": 114.04293, "O": 237.14774,
       "P": 97.05276, "Q": 128.05858, "R": 156.10111, "S": 87.03203, "T": 101.04768, "U": 150.95363, "V": 99.06841,
       "W": 186.07931, "Y": 163.06333}
AA_MASS = np.zeros(32, np.float64)
for _c, _m in _AA.items():
    AA_MASS[ord(_c) - 64] = _m           # token ids 1..26 (A..Z); 0 = pad
OX_TOKEN = 27                            # oxidised M (var mod +15.9949)
AA_MASS[OX_TOKEN] = _AA["M"] + 15.9949
N_TOKENS = 32


@dataclass
class EncoderConfig:
    d: int = 192
    heads: int = 6
    ff: int = 384
    n_spec: int = 2
    n_pep: int = 2
    n_joint: int = 1
    n_peaks: int = 63          # 1 charge token + 63 peaks = 64 keys (SDPA fast path needs key length % 8 == 0)
    emb: int = 128
    n_freq: int = 24
    max_len: int = 63
    tol_ppm: float = 20.0
    dropout: float = 0.0

    def to_dict(self):
        return asdict(self)


SIZES = {  # E2b candidates (parameter counts measured in E2b)
    "xs": EncoderConfig(d=96, heads=4, ff=192, n_spec=1, n_pep=1, n_joint=1),
    "s": EncoderConfig(d=192, heads=6, ff=384, n_spec=2, n_pep=2, n_joint=1),
    "m": EncoderConfig(d=256, heads=8, ff=768, n_spec=3, n_pep=3, n_joint=2),
    "l": EncoderConfig(d=448, heads=8, ff=1536, n_spec=4, n_pep=4, n_joint=2),
    "xl": EncoderConfig(d=640, heads=10, ff=2560, n_spec=4, n_pep=4, n_joint=2),
}


def fourier(x: torch.Tensor, n_freq: int, lo: float = 0.02, hi: float = 4000.0) -> torch.Tensor:
    """x (...,) masses in Da -> (..., 2*n_freq) sin/cos at geometric wavelengths lo..hi (fp32 math)."""
    lam = torch.logspace(math.log10(lo), math.log10(hi), n_freq, device=x.device, dtype=torch.float32)
    a = (2 * math.pi) * x.float().unsqueeze(-1) / lam
    return torch.cat([torch.sin(a), torch.cos(a)], -1)


class Block(nn.Module):
    def __init__(self, c: EncoderConfig, cross: bool = False):
        super().__init__()
        self.h, self.cross = c.heads, cross
        self.ln1 = nn.LayerNorm(c.d)
        self.q = nn.Linear(c.d, c.d)
        self.kv = nn.Linear(c.d, 2 * c.d)
        self.o = nn.Linear(c.d, c.d)
        self.ln2 = nn.LayerNorm(c.d)
        self.ff = nn.Sequential(nn.Linear(c.d, c.ff), nn.GELU(), nn.Linear(c.ff, c.d))
        if cross:
            self.lnk = nn.LayerNorm(c.d)

    def forward(self, x, mem=None, mask=None):
        """x (N, L, d); mem (N, P, d) for cross; mask: additive float (N, H|1, L, P|L) or None."""
        N, L, d = x.shape
        y = self.ln1(x)
        src = self.lnk(mem) if self.cross else y
        q = self.q(y).view(N, L, self.h, -1).transpose(1, 2)
        k, v = self.kv(src).view(N, src.shape[1], 2, self.h, -1).permute(2, 0, 3, 1, 4)
        a = F.scaled_dot_product_attention(q, k, v, attn_mask=mask)
        x = x + self.o(a.transpose(1, 2).reshape(N, L, d))
        return x + self.ff(self.ln2(x))


class SplitJointEncoder(nn.Module):
    def __init__(self, c: EncoderConfig):
        super().__init__()
        self.c = c
        nf = 2 * c.n_freq
        self.peak_in = nn.Linear(nf + 2, c.d)
        self.charge_emb = nn.Embedding(8, c.d)
        self.tok = nn.Embedding(N_TOKENS, c.d, padding_idx=0)
        self.pep_in = nn.Linear(2 * nf, c.d)
        self.pos = nn.Embedding(c.max_len + 1, c.d)
        self.spec_layers = nn.ModuleList([Block(c) for _ in range(c.n_spec)])
        self.pep_layers = nn.ModuleList([Block(c) for _ in range(c.n_pep)])
        self.spec_ln = nn.LayerNorm(c.d)
        self.pep_ln = nn.LayerNorm(c.d)
        self.spec_proj = nn.Linear(c.d, c.emb)
        self.pep_proj = nn.Linear(c.d, c.emb)
        self.joint = nn.ModuleList([Block(c, cross=True) for _ in range(c.n_joint)])
        self.bias_w = nn.Parameter(torch.full((c.n_joint, c.heads, 4), 4.0))
        self.head = nn.Sequential(nn.LayerNorm(2 * c.d + 1), nn.Linear(2 * c.d + 1, c.d), nn.GELU(), nn.Linear(c.d, 1))
        self.logit_scale = nn.Parameter(torch.tensor(math.log(20.0)))

    # ------------------------------------------------------------------ split stage
    def encode_spectra(self, mz, inten, n, charge):
        """mz/inten (S, P) padded with 0, n (S,) peak counts, charge (S,) -> states (S, 1+P, d), key mask, emb (S, emb)."""
        c = self.c
        S, P = mz.shape
        valid = torch.arange(P, device=mz.device)[None] < n[:, None]
        rel = torch.sqrt(inten / inten.amax(1, keepdim=True).clamp_min(1e-12))
        rank = torch.arange(P, device=mz.device, dtype=torch.float32)[None].expand(S, P) / P
        x = self.peak_in(torch.cat([fourier(mz, c.n_freq), rel.unsqueeze(-1), rank.unsqueeze(-1)], -1).to(self.peak_in.weight.dtype))
        x = torch.cat([self.charge_emb(charge.clamp(0, 7)).unsqueeze(1), x], 1)
        kvalid = torch.cat([torch.ones(S, 1, dtype=torch.bool, device=mz.device), valid], 1)
        m = torch.zeros(S, 1, 1, P + 1, device=mz.device, dtype=x.dtype).masked_fill(~kvalid[:, None, None], float("-inf"))
        for b in self.spec_layers:
            x = b(x, mask=m)
        x = self.spec_ln(x)
        w = kvalid.unsqueeze(-1).to(x.dtype)
        e = F.normalize(self.spec_proj((x * w).sum(1) / w.sum(1)).float(), dim=-1)
        return x, kvalid, e

    @staticmethod
    def ion_mz(tokens: torch.Tensor) -> torch.Tensor:
        """tokens (N, L) -> (N, L, 4) theoretical m/z of b1, y1, b2, y2 for the cleavage after each residue
        (the last position gets the full precursor-like ions; they rarely match, harmless)."""
        am = torch.as_tensor(AA_MASS, device=tokens.device)[tokens.long()]       # f64
        pre = torch.cumsum(am, 1)
        tot = pre[:, -1:]
        b = pre
        y = tot - pre + H2O
        out = torch.stack([b + PROTON, y + PROTON, (b + 2 * PROTON) / 2, (y + 2 * PROTON) / 2], -1)
        return out.float()

    def encode_peptides(self, tokens, charge):
        """tokens (N, L) int (0 pad), charge (N,) -> states (N, 1+L, d), valid, emb (N, emb), ion m/z (N, L, 4)."""
        c = self.c
        N, L = tokens.shape
        if (L + 1) % 8:                                        # SDPA fast path: key length multiple of 8
            tokens = F.pad(tokens, (0, 8 - (L + 1) % 8))
            L = tokens.shape[1]
        valid = tokens > 0
        ions = self.ion_mz(tokens)
        x = self.tok(tokens) + self.pep_in(torch.cat([fourier(ions[..., 0] - PROTON, c.n_freq),
                                                      fourier(ions[..., 1] - PROTON, c.n_freq)], -1).to(self.pep_in.weight.dtype))
        x = x + self.pos(torch.arange(L, device=tokens.device))[None]
        x = torch.cat([self.charge_emb(charge.clamp(0, 7)).unsqueeze(1), x], 1)
        kvalid = torch.cat([torch.ones(N, 1, dtype=torch.bool, device=tokens.device), valid], 1)
        m = torch.zeros(N, 1, 1, L + 1, device=tokens.device, dtype=x.dtype).masked_fill(~kvalid[:, None, None], float("-inf"))
        for b in self.pep_layers:
            x = b(x, mask=m)
        x = self.pep_ln(x)
        w = kvalid.unsqueeze(-1).to(x.dtype)
        e = F.normalize(self.pep_proj((x * w).sum(1) / w.sum(1)).float(), dim=-1)
        return x, kvalid, e, ions

    # ------------------------------------------------------------------ joint stage
    def joint_score(self, spec_x, spec_valid, spec_mz, spec_e, pep_x, pep_valid, pep_ions, pep_e, sidx):
        """Score candidates. sidx (N,) maps each candidate to its spectrum row. Returns logits (N,) fp32."""
        c = self.c
        mem = spec_x[sidx]                                   # (N, 1+P, d)
        mz = spec_mz[sidx]                                   # (N, P)
        N, P = mz.shape
        ions = torch.cat([torch.zeros_like(pep_ions[:, :1]), pep_ions], 1)    # charge token: no ions
        d = (ions.unsqueeze(-1) - mz[:, None, None, :]) / (mz[:, None, None, :].clamp_min(1.0) * (c.tol_ppm * 1e-6))
        ker = torch.exp(-0.5 * d * d).to(pep_x.dtype)        # (N, 1+L, 4, P)
        kmask = torch.zeros(N, 1, 1, P + 1, device=mz.device, dtype=pep_x.dtype).masked_fill(
            ~spec_valid[sidx][:, None, None], float("-inf"))
        x = pep_x
        for j, b in enumerate(self.joint):
            bias = torch.einsum("nlip,hi->nhlp", ker, self.bias_w[j].to(ker.dtype))
            bias = torch.cat([torch.zeros_like(bias[..., :1]), bias], -1) + kmask
            x = b(x, mem=mem, mask=bias)
        w = pep_valid.unsqueeze(-1).to(x.dtype)
        mean = (x * w).sum(1) / w.sum(1)
        mx = x.masked_fill(~pep_valid.unsqueeze(-1), float("-inf")).amax(1)
        cos = (spec_e[sidx] * pep_e).sum(-1, keepdim=True).to(x.dtype)
        return self.head(torch.cat([mean, mx, cos], -1)).squeeze(-1).float()

    def forward(self, spec, pep, sidx):
        sx, sv, se = self.encode_spectra(spec["mz"], spec["int"], spec["n"], spec["charge"])
        px, pv, pe, ions = self.encode_peptides(pep["tokens"], pep["charge"])
        return self.joint_score(sx, sv, spec["mz"], se, px, pv, ions, pe, sidx), se, pe


def n_params(m: nn.Module) -> int:
    return sum(p.numel() for p in m.parameters())
