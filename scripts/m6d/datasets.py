"""M6d dev quant sets: paths, runs, conditions and expected log2 ratios (A/B).  Dev only."""
from __future__ import annotations

import os
from dataclasses import dataclass, field

WORK = os.environ.get("VULCAN_M6_WORK", "/tmp/vulcan_m6")
W6D = os.environ.get("VULCAN_M6D_WORK", "/tmp/vulcan_m6d")


@dataclass
class DS:
    name: str
    A: list
    B: list
    expected: dict
    ids: str
    ms1: str
    store_dir: str
    sel: tuple = ()
    extra: dict = field(default_factory=dict)

    def store(self, run):
        return os.path.join(self.store_dir, run + ".vspec.arrow")

    def subset(self, which):
        if which == "all":
            runs = self.A + self.B
        elif which == "sel":
            runs = list(self.sel)
        else:
            runs = [r for r in self.A + self.B if r not in self.sel]
        return runs, [r for r in runs if r in self.A], [r for r in runs if r in self.B]


_D9A = [f"LFQ_Astral_DDA_15min_50ng_Condition_A_REP{i}" for i in (1, 2, 3)]
_D9B = [f"LFQ_Astral_DDA_15min_50ng_Condition_B_REP{i}" for i in (1, 2, 3)]
REG = {
    "d9": DS("d9", _D9A, _D9B, {"HUMAN": 0.0, "YEAST": 1.0, "ECOLI": -2.0}, f"{WORK}/ids/d9_hye", f"{WORK}/ms1",
             f"{WORK}/store", sel=(_D9A[0], _D9A[1], _D9B[0])),
}


def get(name) -> DS:
    if name not in REG:
        import json
        p = os.path.join(W6D, name, "dataset.json")      # written by scripts/m6d/run_set.py for D11 / D12
        d = json.load(open(p))
        REG[name] = DS(**d)
    return REG[name]
