"""Prosit (Koina) MS2 + iRT predictor for the IT-CID instrument class (tree 3c).

AlphaPeptDeep's pretrained MS2 models cover HCD only (instruments QE/Lumos/timsTOF/SciexTOF/ThermoTOF; checked
2026-10-27), so ion-trap CID fragment intensities come from the CID-trained Prosit models served by Koina
(https://koina.wilhelmlab.org, KServe v2 REST):
  label-free   Prosit_2020_intensity_CID   + Prosit_2019_irt
  TMT          Prosit_2020_intensity_TMT (fragmentation_types = "CID", collision energy = 35) + Prosit_2020_irt_TMT
Training data of these models (data/manifests/predictor_training_accessions.tsv): PXD004732, PXD010595, PXD021013,
PXD023119, PXD023120 (ProteomeTools).  None is a dev, sealed or 3c-training accession.  Decision D4 allows Prosit
everywhere.

Interface = vulcan.predict.apd.APD.predict (so vulcan.predict.cache.PredCache caches it unchanged): returns
'ms2' flat f16 in the APD layout (per peptide L-1 rows x [b_z1, b_z2, y_z1, y_z2]; row f = cleavage after residue f,
y counted from the C-term) normalised to max 1, and 'rt' (Prosit iRT units).
Peptides Prosit cannot take (length > 30, residues outside the 20 standard amino acids) get all-zero MS2 and
rt = NaN; the caller treats them as "no prediction" (label-blind: decoys are pseudo-reverse, same length/composition).
"""
from __future__ import annotations

import json
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import numpy as np

URL = "https://koina.wilhelmlab.org/v2/models/%s/infer"
MODELS = {"lf": ("Prosit_2020_intensity_CID", "Prosit_2019_irt"),
          "tmt": ("Prosit_2020_intensity_TMT", "Prosit_2020_irt_TMT")}
TMT_MASS = 229.162932
MAX_LEN = 30
STD = set(b"ACDEFGHIKLMNPQRSTVWY")
BATCH = 1000          # Koina limit per request for these models


def _post(model, inputs, tries=8):
    body = json.dumps({"id": "0", "inputs": inputs}).encode()
    for k in range(tries):
        try:
            req = urllib.request.Request(URL % model, data=body, headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=300) as r:
                out = json.load(r)
            return {o["name"]: o for o in out["outputs"]}
        except Exception:  # noqa: BLE001
            if k == tries - 1:
                raise
            time.sleep(2 * (k + 1))


class KoinaProsit:
    def __init__(self, kind: str, static_mods: dict, var_mods: dict, threads: int = 8, ce: float = 35.0):
        assert kind in MODELS
        self.kind = kind
        self.ms2_model, self.rt_model = MODELS[kind]
        self.static, self.var = dict(static_mods), dict(var_mods)
        self.threads, self.ce = threads, float(ce)
        # PredCache keys on version['model_zip_sha256']; Koina models are identified by name + served version
        self.version = {"model_zip_sha256": f"koina:{self.ms2_model}/1+{self.rt_model}/1", "source": "Koina",
                        "ms2": self.ms2_model, "rt": self.rt_model, "ce": self.ce,
                        "fragmentation": "CID"}
        self.n_requests = 0

    # ------------------------------------------------------------------ peptidoform -> ProForma-like Koina string
    def _strings(self, residues, seq_off, length, mask):
        out, ok = [], np.zeros(len(length), bool)
        for i in range(len(length)):
            L = int(length[i])
            s = bytes(residues[int(seq_off[i]):int(seq_off[i]) + L])
            m = int(mask[i]) & ((1 << 64) - 1)
            if L > MAX_LEN or any(c not in STD for c in s):
                out.append(None)
                continue
            parts = []
            if self.kind == "tmt":
                parts.append("[UNIMOD:737]-")
            for k, c in enumerate(s.decode()):
                tok = c
                if (m >> k) & 1:
                    if c == "M" and abs(self.var.get("M", 0) - 15.9949) < 1e-3:
                        tok = "M[UNIMOD:35]"
                    else:
                        out.append(None); break
                elif c == "C" and abs(self.static.get("C", 0) - 57.021464) < 1e-3:
                    tok = "C[UNIMOD:4]"
                elif c == "K" and self.kind == "tmt":
                    tok = "K[UNIMOD:737]"
                parts.append(tok)
            else:
                out.append("".join(parts))
                ok[i] = True
                continue
        return out, ok

    def _run(self, model, strs, charges=None):
        idx = [i for i, s in enumerate(strs) if s is not None]
        batches = [idx[b:b + BATCH] for b in range(0, len(idx), BATCH)]

        def one(bi):
            P = [strs[i] for i in bi]
            n = len(P)
            inp = [{"name": "peptide_sequences", "shape": [n, 1], "datatype": "BYTES", "data": P}]
            if charges is not None:
                inp.append({"name": "precursor_charges", "shape": [n, 1], "datatype": "INT32",
                            "data": [int(min(max(charges[i], 1), 6)) for i in bi]})
                if self.kind == "tmt":
                    inp += [{"name": "collision_energies", "shape": [n, 1], "datatype": "FP32", "data": [self.ce] * n},
                            {"name": "fragmentation_types", "shape": [n, 1], "datatype": "BYTES", "data": ["CID"] * n}]
            return bi, _post(model, inp)

        res = []
        if batches:
            with ThreadPoolExecutor(self.threads) as ex:
                res = list(ex.map(one, batches))
        self.n_requests += len(batches)
        return res

    def predict(self, residues, seq_off, length, mask, charge, nce, instrument="IT", batch=None,
                want_ms2=True, want_rt=True):
        n = len(length)
        length = np.asarray(length, np.int64)
        strs, ok = self._strings(residues, seq_off, length, np.asarray(mask).view(np.uint64))
        out_rt = np.full(n, np.nan, np.float32)
        fo = np.zeros(n + 1, np.int64)
        np.cumsum(np.maximum(length - 1, 0) * 4, out=fo[1:])
        out_ms2 = np.zeros(fo[-1], np.float16) if want_ms2 else None
        if want_rt:
            for bi, o in self._run(self.rt_model, strs):
                out_rt[bi] = np.asarray(o["irt"]["data"], np.float32).reshape(-1)
        if want_ms2:
            ch = np.asarray(charge, np.int64)
            for bi, o in self._run(self.ms2_model, strs, ch):
                inten = np.asarray(o["intensities"]["data"], np.float32).reshape(len(bi), 29, 2, 3)  # pos, (y,b), z
                inten = np.where(inten < 0, 0.0, inten)
                for j, i in enumerate(bi):
                    L = int(length[i])
                    rows = np.zeros((L - 1, 4), np.float32)
                    k = np.arange(1, L)                      # ion number 1..L-1
                    # b_k -> row k-1 ; y_k -> row L-1-k ; z1/z2 only (z3 dropped)
                    rows[k - 1, 0] = inten[j, k - 1, 1, 0]
                    rows[k - 1, 1] = inten[j, k - 1, 1, 1]
                    rows[L - 1 - k, 2] = inten[j, k - 1, 0, 0]
                    rows[L - 1 - k, 3] = inten[j, k - 1, 0, 1]
                    mx = rows.max()
                    if mx > 0:
                        rows /= mx
                    out_ms2[fo[i]:fo[i + 1]] = rows.reshape(-1).astype(np.float16)
        return {"rt": out_rt, "ms2": out_ms2, "ms2_off": fo}
