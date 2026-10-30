"""Semi-specific L1 index (DESIGN §4 L1): the standard paired-decoy builder (vulcan.index.build.build) fed with the
semi-specific digest of vulcan.lip.digest instead of vulcan_core's specific digest.

    python -m vulcan.lip.index --fasta a.fasta --out DIR [--missed 2] [--entrap-seeds ...]

Decoys.  Every forward peptide (specific or semi) gets its paired decoy by the standard recipe (pseudo-reverse with
the C-terminal residue fixed; collisions -> hash shuffles), and the specificity class is a property of the PAIR
(taken from the forward member), written to `lip_class.parquet` keyed by target sequence.  Hence #decoys = #targets
in every class, which is what DESIGN L1 asks of "decoys from decoy proteins".  Literal semi-specific digestion of
reversed proteins was not used because (i) it breaks the 1:1 target-decoy pairing that paired double competition,
G0 and family-max entrapment rely on, and (ii) for C-tryptic half-tryptic targets (ending in K/R) a reversed-protein
decoy starts with K/R instead of ending with it, so its y-ion series is systematically different (not exchangeable);
the C-term-fixed pseudo-reverse keeps the terminal residue of every class.

Index-side output is the standard layout (peptidoforms.arrow, residues.npy, labels.arrow, proteins.arrow, ...), so
the GPU engine, scorer and L4 code run unchanged.  `missed` = internal cleavage sites of the most specific occurrence.
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from vulcan.index import build as B
from vulcan.index._core import core as _real_core

from .digest import digest

SEMI_VERSION = "lip-semi-0.1"


class _CoreProxy:
    """vulcan_core with digest_fasta replaced by a precomputed digest (everything else delegated)."""

    def __init__(self, dg):
        self._dg = dg

    def digest_fasta(self, *a, **k):
        return self._dg

    def __getattr__(self, name):
        return getattr(_real_core, name)


def build_semi(fasta_paths, out_dir, missed=2, min_len=7, max_len=50, entrap_seeds=(), semi=True, workers=None,
               var_mods=None, max_var=2, shuffle_rows_seed=None, swap_labels=False) -> dict:
    """shuffle_rows_seed / swap_labels: G0 test knobs passed to vulcan.index.build.build (index side must not change)."""
    t0 = time.perf_counter()
    dg = digest(fasta_paths, missed=missed, min_len=min_len, max_len=max_len, semi=semi, workers=workers)
    t_dig = time.perf_counter() - t0
    params = {"enzyme": {"cleave_at": "KR", "restrict": "P", "c_terminal": True, "missed_cleavages": missed,
                         "min_len": min_len, "max_len": max_len, "semi": SEMI_VERSION if semi else False},
              "entrap_seeds": list(entrap_seeds)}
    if var_mods is not None:
        params["variable_mods"] = var_mods
        params["max_variable_mods"] = max_var
    old = B.core
    B.core = _CoreProxy(dg)
    try:
        meta = B.build(fasta_paths, out_dir, params, workers=workers, shuffle_rows_seed=shuffle_rows_seed,
                       swap_labels=swap_labels)
    finally:
        B.core = old
    pq.write_table(pa.table({"target_seq": pa.array(dg["keys"], pa.string()), "cls": pa.array(dg["cls"], pa.uint8())}),
                   os.path.join(out_dir, "lip_class.parquet"))
    cnt = np.bincount(dg["cls"], minlength=3)
    meta["semi"] = {"version": SEMI_VERSION if semi else None, "digest_s": t_dig,
                    "n_targets_by_class": {"specific": int(cnt[0]), "n_semi": int(cnt[1]), "c_semi": int(cnt[2])}}
    meta["index_version"] = B.INDEX_VERSION + ("+" + SEMI_VERSION if semi else "+pydigest")
    with open(os.path.join(out_dir, "meta.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    return meta


def load_class(index_dir) -> dict:
    t = pq.read_table(os.path.join(index_dir, "lip_class.parquet"))
    return dict(zip(t.column("target_seq").to_pylist(), t.column("cls").to_numpy().tolist()))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--missed", type=int, default=2)
    ap.add_argument("--specific", action="store_true", help="fully specific digest (python digester) for checks")
    ap.add_argument("--entrap-seeds", nargs="*", type=int, default=[])
    ap.add_argument("--workers", type=int, default=None)
    a = ap.parse_args(argv)
    meta = build_semi(a.fasta, a.out, missed=a.missed, entrap_seeds=a.entrap_seeds, semi=not a.specific,
                      workers=a.workers)
    print(json.dumps({k: v for k, v in meta.items() if k not in ("params", "fasta")}, indent=1))


if __name__ == "__main__":
    main()
