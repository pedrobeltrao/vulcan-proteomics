"""Loading PSM + L3 feature tables and (L4-only) labels for M3 scripts."""
from __future__ import annotations

import glob
import os

import numpy as np
import pyarrow.ipc as ipc
import pyarrow.parquet as pq

PSM_COLS = ["run", "spec_idx", "rank", "pf_row", "pf_id", "charge", "isotope_error", "hyperscore", "delta_next",
            "delta_best", "matched_b", "matched_y", "matched_peaks", "matched_intensity_pct", "average_ppm",
            "delta_mass_ppm", "longest_b", "longest_y", "longest_y_pct", "peptide_len", "missed_cleavages", "n_var",
            "scored_candidates", "poisson", "rt"]


def load_runs(d, with_feat=True):
    """-> list of (run, psm dict, feat dict|None) for every <run>.psms.parquet in d."""
    out = []
    for f in sorted(glob.glob(os.path.join(d, "*.psms.parquet"))):
        run = os.path.basename(f)[:-len(".psms.parquet")]
        t = pq.read_table(f)
        psm = {c: t.column(c).to_numpy(zero_copy_only=False) for c in PSM_COLS}
        feat = None
        if with_feat:
            ff = os.path.join(d, run + ".feat.parquet")
            ft = pq.read_table(ff)
            feat = {c: ft.column(c).to_numpy() for c in ft.column_names}
        out.append((run, psm, feat))
    return out


def labels(index_dir):
    """L4 ONLY."""
    with ipc.open_file(os.path.join(index_dir, "labels.arrow")) as r:
        t = r.read_all()
    cols = ("pf_id", "family", "entrap_seed", "is_decoy", "is_entrapment", "tdc_pair", "entrap_pair", "target_seq")
    return {k: t.column(k).to_numpy(zero_copy_only=False) for k in cols}


def concat(dicts):
    return {k: np.concatenate([d[k] for d in dicts]) for k in dicts[0]}
