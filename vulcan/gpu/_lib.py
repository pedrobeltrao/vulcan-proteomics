"""ctypes binding of tools/vulcan_kernels/libvulcan_kernels.so (built from vulcan/gpu/kernels.cu).

The .so is (re)built by scripts/m2/build_kernels.sh whenever the recorded source hash differs."""
from __future__ import annotations

import ctypes as C
import hashlib
import os
import subprocess

import torch

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SRC = os.path.join(ROOT, "vulcan", "gpu", "kernels.cu")
OUT = os.path.join(ROOT, "tools", "vulcan_kernels")
SO = os.path.join(OUT, "libvulcan_kernels.so")
BUILD = os.path.join(ROOT, "scripts", "m2", "build_kernels.sh")

P = C.c_void_p
I32, I64, F32, F64 = C.c_int, C.c_int64, C.c_float, C.c_double

_SIGS = {
    "vk_set_params": [P, P, P],
    "vk_frag_count": [I64, P, P, P, P, P, I32, F32, F32, P, P],
    "vk_frag_emit": [I64, P, P, P, P, P, I32, F32, F32, F64, I32, I64, P, P, P],
    "vk_k1": [I32, P, P, P, P, P, P, P, P, P, P, I32, I32, I32, I32, F32, F32, F32, F64, I32, I64, P, P, P, P,
              P, P, P],
    "vk_k1_brute": [I64, P, P, P, P, P, P, P, P, P, P, P, I32, F32, F32, F32, F32, F64, P, P],
    "vk_k2": [I32, P, P, P, P, P, P, P, P, P, P, P, P, F32, F32, P, P, P, P, P, P, P, P, I32, P],
    # M4: variable peptide-N-term mod (mask bit 63) and K2 with per-candidate precursor delta (open/offset modes)
    "vk_set_nterm": [F32],
    "vk_k2_shift": [I32, P, P, P, P, P, P, P, P, P, P, P, P, F32, F32, P, P, P, P, P, P, P, P, P, I32, P],
    # 3c: fragment tolerance unit (0 = ppm, 1 = Da; IT-CID)
    "vk_set_tol_da": [I32],
}

_lib = None


def _src_hash() -> str:
    with open(SRC, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def lib():
    global _lib
    if _lib is not None:
        return _lib
    stamp = os.path.join(OUT, "source.sha256")
    fresh = os.path.exists(SO) and os.path.exists(stamp) and open(stamp).read().strip() == _src_hash()
    if not fresh:
        subprocess.run(["bash", BUILD], check=True)
    torch.cuda.init()
    L = C.CDLL(SO)
    for name, args in _SIGS.items():
        f = getattr(L, name)
        f.argtypes = args
        f.restype = C.c_int
    _lib = L
    return L


def ptr(t: torch.Tensor | None):
    if t is None:
        return None
    assert t.is_cuda and t.is_contiguous(), "kernels need contiguous CUDA tensors"
    return t.data_ptr()


def stream():
    return torch.cuda.current_stream().cuda_stream


def check(code: int, what: str):
    if code != 0:
        raise RuntimeError(f"{what}: CUDA error {code}")
