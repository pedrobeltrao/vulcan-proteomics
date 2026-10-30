"""AlphaPeptDeep (peptdeep, Apache-2.0) MS2 + RT inference on the GPU with a vectorized input path.

peptdeep's own `predict()` builds every mini-batch with pandas (groupby nAA, str.split for mods), which caps the
A100 at ~11k precursors/s (MS2) and ~3.4k peptides/s (RT).  Here the same models (generic pretrained_models_v3)
are fed directly: residues (u8 ASCII) + variable-mod bitmask (vulcan L1 index layout) -> aa indices and mod
features built on the GPU, length-bucketed exactly like peptdeep (no padding inside a batch, so outputs are
identical to `ModelManager.*.predict`; checked by tests/predict).

Decision D4: AlphaPeptDeep features are NEVER used on D2/T2 (PXD009449/PXD000138 are in its training data);
`check_allowed(dataset_id)` enforces it.  Environment: conda_envs/pred (pe + peptdeep 1.5.1 --no-deps + transformers).
Model files: data/cache/peptdeep/pretrained_models/pretrained_models_v3.zip (sha256 345cb241...).
"""
from __future__ import annotations

import hashlib
import os

import numpy as np
import torch

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PEPTDEEP_HOME = os.path.join(ROOT, "data", "cache", "peptdeep")
MODEL_ZIP = os.path.join(PEPTDEEP_HOME, "pretrained_models", "pretrained_models_v3.zip")
FORBIDDEN_SETS = {"D2", "T2"}                 # decision D4
FORBIDDEN_ACCESSIONS = {"PXD009449", "PXD000138"}
FRAG_TYPES = ("b_z1", "b_z2", "y_z1", "y_z2")


def check_allowed(dataset_id: str = "", accession: str = ""):
    if dataset_id in FORBIDDEN_SETS or accession in FORBIDDEN_ACCESSIONS:
        raise PermissionError(f"decision D4: AlphaPeptDeep features are not allowed on {dataset_id or accession}")


def model_sha256() -> str:
    h = hashlib.sha256()
    with open(MODEL_ZIP, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


class APD:
    """Generic AlphaPeptDeep MS2 (b/y z1/z2) and RT models on one device."""

    def __init__(self, device="cuda", static_mods=None, var_mods=None):
        from peptdeep.settings import global_settings as g
        g["PEPTDEEP_HOME"] = PEPTDEEP_HOME
        import peptdeep
        from peptdeep.pretrained_models import ModelManager
        from peptdeep.model.featurize import MOD_TO_FEATURE, parse_instrument_indices
        self.version = {"peptdeep": peptdeep.__version__, "models": "generic (pretrained_models_v3.zip)",
                        "model_zip_sha256": model_sha256()}
        mm = ModelManager(mask_modloss=True, device="gpu" if device == "cuda" else "cpu")
        mm.load_installed_models("generic")
        self.mm = mm
        self.dev = torch.device(device)
        self.ms2 = mm.ms2_model.model.to(self.dev).eval()
        self.rt = mm.rt_model.model.to(self.dev).eval()
        self.charge_factor = mm.ms2_model.charge_factor
        self.nce_factor = mm.ms2_model.NCE_factor
        self.min_inten = mm.ms2_model.min_inten
        sup = list(self.ms2.supported_charged_frag_types) if hasattr(self.ms2, "supported_charged_frag_types") \
            else list(mm.ms2_model.charged_frag_types)
        self.col = [sup.index(f) for f in FRAG_TYPES]
        self.parse_inst = parse_instrument_indices
        static_mods = {"C": 57.021464} if static_mods is None else static_mods
        var_mods = {"M": 15.9949} if var_mods is None else var_mods
        names = {("C", 57.021464): "Carbamidomethyl@C", ("M", 15.9949): "Oxidation@M"}
        F = len(next(iter(MOD_TO_FEATURE.values())))
        self.F = F
        st = np.zeros((128, F), np.float32)
        for aa, m in static_mods.items():
            st[ord(aa)] = MOD_TO_FEATURE[names[(aa, m)]]
        vr = np.zeros((128, F), np.float32)
        for aa, m in var_mods.items():
            vr[ord(aa)] = MOD_TO_FEATURE[names[(aa, m)]]
        self.static_feat = torch.from_numpy(st).to(self.dev)
        self.var_feat = torch.from_numpy(vr).to(self.dev)

    # ---------------------------------------------------------------- inputs
    def _inputs(self, res: torch.Tensor, mask: torch.Tensor):
        """res: [B, L] u8 ASCII (uniform L); mask: [B] int64 var-mod bitmask -> (aa_idx [B,L+2], mod_x [B,L+2,F])."""
        B, L = res.shape
        r = res.long()
        aa = torch.zeros(B, L + 2, dtype=torch.long, device=self.dev)
        aa[:, 1:L + 1] = r - ord("A") + 1
        bits = ((mask[:, None] >> torch.arange(L, device=self.dev)[None, :]) & 1).bool()
        mx = torch.zeros(B, L + 2, self.F, device=self.dev)
        mx[:, 1:L + 1] = self.static_feat[r] + torch.where(bits[..., None], self.var_feat[r], 0.0)
        return aa, mx

    @torch.inference_mode()
    def predict(self, residues: np.ndarray, seq_off: np.ndarray, length: np.ndarray, mask: np.ndarray,
                charge: np.ndarray | None, nce: float, instrument: str = "QE", batch: int = 16384,
                want_ms2: bool = True, want_rt: bool = True):
        """Predict for n peptidoforms (rows given by offsets into `residues`).

        Returns dict with 'rt' [n] float32 (peptdeep normalized RT) and 'ms2' (flat float16 [sum(L-1)*4]) +
        'ms2_off' [n+1] (offsets in fragments; fragment f of row i = cleavage after residue f, columns
        b_z1,b_z2,y_z1,y_z2 with y counted from the C-term as in alphabase).
        """
        n = len(length)
        length = np.asarray(length, np.int64)
        out_rt = np.zeros(n, np.float32)
        fo = np.zeros(n + 1, np.int64)
        np.cumsum(np.maximum(length - 1, 0) * 4, out=fo[1:])
        out_ms2 = np.zeros(fo[-1], np.float16) if want_ms2 else None
        res_t = torch.from_numpy(np.ascontiguousarray(residues)).to(self.dev)
        inst = int(self.parse_inst([instrument])[0])
        order = np.argsort(length, kind="stable")
        ls = length[order]
        bounds = np.flatnonzero(np.diff(ls)) + 1
        starts = np.r_[0, bounds]
        ends = np.r_[bounds, n]
        for s0, e0 in zip(starts, ends):
            L = int(ls[s0])
            for b0 in range(s0, e0, batch):
                b1 = min(e0, b0 + batch)
                rows = order[b0:b1]
                rows_t = torch.from_numpy(rows).to(self.dev)
                off = torch.from_numpy(np.asarray(seq_off, np.int64)[rows]).to(self.dev)
                res = res_t[off[:, None] + torch.arange(L, device=self.dev)[None, :]]
                aa, mx = self._inputs(res, torch.from_numpy(np.asarray(mask, np.int64)[rows]).to(self.dev))
                if want_rt:
                    out_rt[rows] = self.rt(aa, mx).float().cpu().numpy().reshape(-1)
                if want_ms2:
                    B = len(rows)
                    ch = torch.from_numpy(np.asarray(charge)[rows].astype(np.float32)).to(self.dev)[:, None] \
                        * self.charge_factor
                    nc = torch.full((B, 1), float(nce) * self.nce_factor, device=self.dev)
                    ins = torch.full((B,), inst, dtype=torch.long, device=self.dev)
                    p = self.ms2(aa, mx, ch, nc, ins)                       # [B, L-1, n_types]
                    p = p[:, :, self.col]
                    apex = p.reshape(B, -1).amax(1)
                    apex = torch.where(apex <= 0, torch.ones_like(apex), apex)
                    p = p / apex[:, None, None]
                    p = torch.where(p < self.min_inten, torch.zeros_like(p), p)
                    p = p.reshape(B, -1).half().cpu().numpy()
                    # scatter rows of equal length into the flat output
                    idx = fo[rows][:, None] + np.arange((L - 1) * 4)[None, :]
                    out_ms2[idx.ravel()] = p.ravel()
        return {"rt": out_rt, "ms2": out_ms2, "ms2_off": fo}
