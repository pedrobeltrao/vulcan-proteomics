"""Frozen cross-encoder scoring of PSM tables (L3, label-blind) and split-stage embeddings (E5).

Deterministic: rows are processed in table order in fixed-size row blocks, so identical PSM tables (G0) give
bit-identical scores.  Inputs: the PSM table (pf_row, charge, spec_idx), the L1 index side and the run's spectra."""
from __future__ import annotations

import hashlib
import time

import numpy as np
import torch

from .data import peptide_tokens, spectra_topk_gpu
from .model import EncoderConfig, SplitJointEncoder

ROW_BLOCK = 8192


def load_model(path: str, dtype=torch.bfloat16, device="cuda"):
    blob = open(path, "rb").read()
    import io
    ck = torch.load(io.BytesIO(blob), map_location="cpu", weights_only=False)
    m = SplitJointEncoder(EncoderConfig(**ck["config"]))
    m.load_state_dict(ck["state_dict"])
    m = m.to(device).to(dtype).eval()
    m.sha256 = hashlib.sha256(blob).hexdigest()
    return m


class SpectraTensors:
    def __init__(self, sp, device="cuda"):
        # GPU peak selection, bit-identical to data.spectra_topk (verified; the CPU version was 4.7 core-s on D0)
        self.mz, self.int, self.n = spectra_topk_gpu(sp, device=device)


@torch.inference_mode()
def score_table(model, l1, psm: dict, sp, device="cuda", block: int = ROW_BLOCK):
    """psm: pf_row, charge, spec_idx (numpy).  Returns (logits float32 (n,), info)."""
    t0 = time.perf_counter()
    st = SpectraTensors(sp, device)
    n = len(psm["pf_row"])
    out = np.zeros(n, np.float32)
    tok_all = peptide_tokens(l1, psm["pf_row"])
    t_prep = time.perf_counter() - t0
    torch.cuda.synchronize()
    ev0, ev1 = torch.cuda.Event(enable_timing=True), torch.cuda.Event(enable_timing=True)
    ev0.record()
    for b in range(0, n, block):
        sl = slice(b, min(n, b + block))
        si = psm["spec_idx"][sl].astype(np.int64)
        u, inv = np.unique(si, return_inverse=True)
        ut = torch.from_numpy(u).to(device)
        ch = torch.from_numpy(psm["charge"][sl].astype(np.int64)).to(device)
        # spectrum charge token: instrument charge, else the charge of the spectrum's first row in the block (label-free)
        first = np.zeros(len(u), np.int64)
        first[inv[::-1]] = np.arange(len(inv))[::-1]
        zc = np.where(sp.charge[u] > 0, sp.charge[u].astype(np.int64), psm["charge"][sl][first].astype(np.int64))
        spec = {"mz": st.mz[ut], "int": st.int[ut], "n": st.n[ut], "charge": torch.from_numpy(zc).to(device)}
        tok = tok_all[sl]
        L = int((tok > 0).sum(1).max())
        pep = {"tokens": torch.from_numpy(tok[:, :L].astype(np.int64)).to(device), "charge": ch}
        logit, _, _ = model(spec, pep, torch.from_numpy(inv.astype(np.int64)).to(device))
        out[sl] = logit.float().cpu().numpy()
    ev1.record()
    torch.cuda.synchronize()
    return out, {"n_rows": n, "n_spectra": len(sp), "prep_s": t_prep, "gpu_s": ev0.elapsed_time(ev1) / 1e3,
                 "wall_s": time.perf_counter() - t0}


@torch.inference_mode()
def embed_peptides(model, l1, rows: np.ndarray, charge: int = 2, block: int = 65536, device="cuda"):
    """Split-stage peptide embeddings (fp16, L2-normalised) for L1 rows."""
    out = torch.empty((len(rows), model.c.emb), dtype=torch.float16, device=device)
    for b in range(0, len(rows), block):
        tok = peptide_tokens(l1, rows[b:b + block])
        L = int((tok > 0).sum(1).max())
        t = torch.from_numpy(tok[:, :L].astype(np.int64)).to(device)
        _, _, e, _ = model.encode_peptides(t, torch.full((len(t),), charge, device=device))
        out[b:b + len(t)] = e.half()
    return out


@torch.inference_mode()
def embed_spectra(model, sp, charge: np.ndarray, block: int = 16384, device="cuda"):
    st = SpectraTensors(sp, device)
    S = len(sp)
    out = torch.empty((S, model.c.emb), dtype=torch.float16, device=device)
    ch = torch.from_numpy(np.asarray(charge, np.int64)).to(device)
    for b in range(0, S, block):
        sl = slice(b, min(S, b + block))
        _, _, e = model.encode_spectra(st.mz[sl], st.int[sl], st.n[sl], ch[sl])
        out[sl] = e.half()
    return out
