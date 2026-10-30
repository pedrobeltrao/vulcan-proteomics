"""Static-first scorers (DESIGN §4 L3 "Default scorer: static-first"): a FROZEN low-capacity model per instrument
class over label-blind + predicted features.  Trained once on project-disjoint training data (scripts/m3/
train_static.py); applying it never looks at labels, so plain paired TDC / picked protein remain valid on the
evaluated data (no RESET needed, §6.4).

Feature matrix (`feature_matrix`) is a pure function of the PSM table + L3 feature table of ONE run.
"""
from __future__ import annotations

import hashlib
import json
import pickle

import numpy as np

BASE_FEATURES = ["hyperscore", "hs_runz", "delta_next", "delta_best", "matched_peaks", "matched_b", "matched_y",
                 "matched_intensity_pct", "abs_average_ppm", "abs_delta_mass_ppm", "longest_b", "longest_y",
                 "longest_y_pct", "peptide_len", "missed_cleavages", "n_var", "z2", "z3", "z4p", "isotope_error",
                 "log_scored_candidates", "poisson"]
PRED_FEATURES = ["pred_sa", "pred_sa_sqrt", "pred_pearson", "pred_entropy", "pred_frac_obs", "pred_top6",
                 "pred_obs_tic", "rt_abs_resid", "rt_resid_z"]
# NOT used as a model feature: rt_pred (a pure sequence property; target/decoy AUC 0.71 on D1 -> G1-type leak risk)


def feature_matrix(psm: dict, feat: dict | None, names) -> np.ndarray:
    """psm/feat: dicts of numpy columns for one run (all ranks)."""
    r1 = psm["rank"] == 1
    hs = psm["hyperscore"].astype(np.float64)
    med = np.median(hs[r1]) if r1.any() else 0.0
    iqr = (np.quantile(hs[r1], 0.9) - med) if r1.any() else 1.0
    z = psm["charge"].astype(np.int64)
    d = {
        "hs_runz": (hs - med) / max(iqr, 1e-6),
        "abs_average_ppm": np.abs(psm["average_ppm"]),
        "abs_delta_mass_ppm": np.abs(psm["delta_mass_ppm"]),
        "z2": (z == 2).astype(float), "z3": (z == 3).astype(float), "z4p": (z >= 4).astype(float),
        "log_scored_candidates": np.log1p(psm["scored_candidates"].astype(np.float64)),
    }
    cols = []
    for n in names:
        if n in d:
            v = d[n]
        elif n in psm:
            v = psm[n]
        else:
            v = feat[n]
        cols.append(np.asarray(v, dtype=np.float64))
    X = np.column_stack(cols)
    return np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)


class StaticModel:
    """Frozen scorer: kind 'linear' (standardized logistic regression) or 'gbm' (sklearn HistGradientBoosting)."""

    def __init__(self, kind: str, names, model, mu=None, sd=None, meta=None):
        self.kind, self.names, self.model, self.mu, self.sd = kind, list(names or []), model, mu, sd
        self.meta = meta or {}

    def score(self, X: np.ndarray) -> np.ndarray:
        if self.kind == "linear":
            Z = (X - self.mu) / self.sd
            return Z @ self.model["w"] + self.model["b"]
        return self.model.decision_function(X)

    def apply(self, psm, feat):
        return self.score(feature_matrix(psm, feat, self.names))

    def save(self, path):
        blob = pickle.dumps(self)
        with open(path, "wb") as fh:
            fh.write(blob)
        meta = dict(self.meta, sha256=hashlib.sha256(blob).hexdigest(), kind=self.kind, features=self.names)
        if self.kind == "linear":
            meta["weights"] = dict(zip(self.names, map(float, self.model["w"])))
            meta["bias"] = float(self.model["b"])
            meta["mu"] = list(map(float, self.mu))
            meta["sd"] = list(map(float, self.sd))
        with open(path + ".json", "w") as fh:
            json.dump(meta, fh, indent=1, default=float)
        return meta["sha256"]

    @staticmethod
    def load(path):
        with open(path, "rb") as fh:
            return pickle.load(fh)
