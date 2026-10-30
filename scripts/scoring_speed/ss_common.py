"""Shared settings for the scoring-speed study (results/2026-10-27_scoring_speed).  Dev data only.

Runs: D0 (2 runs), D1 (1 run), D8 holdout (first 10 runs of data/manifests/e9_d8.tsv, M3e stores).  All searched
against the M3 production indexes (human_contam_prod / pfu_contam_prod), report_psms 5, as in M3/M3e.
"""
from __future__ import annotations

import glob
import os
import resource
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts", "m3"))
sys.path.insert(0, ROOT)
import m3_common as C  # noqa: E402

W = os.environ.get("VULCAN_SS_WORK", "/tmp/vulcan_ss")
RES = os.path.join(ROOT, "results", "2026-10-27_scoring_speed")
M3E_W = "/tmp/vulcan_m3e"
REF_CACHE = os.path.join(C.WORK, "predcache")          # the M3 cache root the reference features were made from
LIB29_KEY = "a9d5af1427616250"
MODEL = os.path.join(C.RES, "models", "static_gbm_pred.pkl")


def d0_stores():
    return [C.E2.store_path(r) for r in C.DEV["D0"]["raw"]]


def d1_stores():
    return [C.E2.store_path(r) for r in C.DEV["D1"]["raw"]]


def d8_stores():
    return sorted(glob.glob(os.path.join(M3E_W, "store_d8", "*.vspec.arrow")))


def run_name(store):
    return os.path.basename(store)[:-len(".vspec.arrow")]


def ref_psm_dir(ds):
    return {"D0": os.path.join(C.WORK, "psms", "D0_prod"), "D1": os.path.join(C.WORK, "psms", "D1_prod"),
            "D8": os.path.join(M3E_W, "psms", "D8_prod")}[ds]


def family(ds):
    return "pfu_contam" if ds == "D1" else "human_contam"


def plan():
    """[(tag, dataset, store)] in the M3e bench order: D0 run 1 warm-up, 10 D8 runs, D0 run 2, then D1."""
    d0 = d0_stores()
    return ([("warmup", "D0", d0[0])] + [("D8", "D8", f) for f in d8_stores()] + [("D0", "D0", d0[1])]
            + [("D1", "D1", f) for f in d1_stores()])


def cpu():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


class Acc:
    """Accumulating wall/core-s timers keyed by name (cuda-synchronized at both ends)."""

    def __init__(self, torch=None):
        self.torch, self.rec = torch, {}

    def reset(self):
        self.rec = {}

    def __call__(self, name):
        acc = self

        class _C:
            def __enter__(s):
                if acc.torch is not None:
                    acc.torch.cuda.synchronize()
                s.w, s.c = time.perf_counter(), cpu()

            def __exit__(s, *a):
                if acc.torch is not None:
                    acc.torch.cuda.synchronize()
                r = acc.rec.setdefault(name, {"wall_s": 0.0, "cpu_core_s": 0.0, "n": 0})
                r["wall_s"] += time.perf_counter() - s.w
                r["cpu_core_s"] += cpu() - s.c
                r["n"] += 1
        return _C()
