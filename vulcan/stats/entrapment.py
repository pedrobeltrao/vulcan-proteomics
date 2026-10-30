"""Entrapment FDP estimators and the pre-registered decision rule (DESIGN.md §6.3, v0.3).

Estimators (Wen, Freestone, Riffle, MacCoss, Noble, Keich, Nat Methods 22:1454 (2025);
code: github.com/Noble-Lab/FDRBench, ``FDPCalc.java``):

* lower bound   LB      = N_E / (N_T + N_E)                       (failure detector only)
* combined      FDP_C   = N_E (1 + 1/r) / (N_T + N_E)             (estimated upper bound)
* paired (r=1)  FDP_P   = (N_E + N_{E>=s>T} + 2 N_{E>T>=s}) / (N_T + N_E)

Here N_T = original-target discoveries, N_E = entrapment discoveries (decoys are
never counted), s = the score cutoff of the reported list, N_{E>=s>T} = discovered
entrapments whose paired target is *not* discovered, N_{E>T>=s} = discovered
entrapments whose paired target is also discovered with a lower score. This
matches FDRBench ``paired_fdp = (n_entrapment_targets + 2 n_p_t_s + n_p_s_t) /
(n_entrapment_targets + n_targets)`` (verified against the source,
2026-10-26). Exact score ties between an entrapment and its paired target are
counted as E > T (the conservative choice, weight 2).

Decision rule (v0.3): replicate = (entrapment seed, decoy seed); statistic = mean
FDP over replicates; 95% CI from a Poisson bootstrap of target-entrapment pairs
(the same pair weights for all replicates, keyed by the pair id) combined with
between-seed variance:

    Var(mean) = Var_boot(mean over replicates) + s^2_between / J,

    CI = mean +- 1.96 sqrt(Var(mean)).

| Verdict          | Condition                                            |
|------------------|------------------------------------------------------|
| not-evaluable    | alpha * N_T < 20   (N_T = mean over replicates)      |
| violated         | lower 95% CI of the mean LB > alpha                  |
| supported        | mean FDP <= alpha and upper 95% CI <= 1.5 alpha      |
| supported-noisy  | mean FDP <= alpha and upper 95% CI > 1.5 alpha       |
| inconclusive     | otherwise                                            |

Gated levels: alpha in {1%, 5%}; other levels are reported, not gated.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np

from ._util import as_codes, make_rng

__all__ = [
    "lower_bound", "combined", "paired",
    "EntrapmentCounts", "paired_counts", "combined_counts",
    "protein_group_is_target", "protein_paired_counts",
    "Replicate", "evaluate_arm", "combine_arms", "verdict",
    "GATED_ALPHAS", "VERDICTS",
]

GATED_ALPHAS = (0.01, 0.05)
VERDICTS = ("supported", "supported-noisy", "violated", "not-evaluable", "inconclusive")


def lower_bound(n_t, n_e):
    n_t, n_e = np.asarray(n_t, float), np.asarray(n_e, float)
    return np.where(n_t + n_e > 0, n_e / np.maximum(n_t + n_e, 1), 0.0)


def combined(n_t, n_e, r: float):
    n_t, n_e = np.asarray(n_t, float), np.asarray(n_e, float)
    return np.where(n_t + n_e > 0, n_e * (1.0 + 1.0 / r) / np.maximum(n_t + n_e, 1), 0.0)


def paired(n_t, n_e, n_e_only, n_e_gt_t):
    """(N_E + N_{E>=s>T} + 2 N_{E>T>=s}) / (N_T + N_E)."""
    n_t, n_e = np.asarray(n_t, float), np.asarray(n_e, float)
    num = n_e + np.asarray(n_e_only, float) + 2.0 * np.asarray(n_e_gt_t, float)
    return np.where(n_t + n_e > 0, num / np.maximum(n_t + n_e, 1), 0.0)


@dataclass
class EntrapmentCounts:
    """Per-pair contributions (so that pairs can be bootstrapped)."""

    key: np.ndarray       # pair key (one row per pair touching the list)
    t: np.ndarray         # original target of the pair accepted (0/1)
    e: np.ndarray         # entrapment of the pair accepted (0/1)
    e_only: np.ndarray    # E accepted, paired T not accepted
    e_gt_t: np.ndarray    # both accepted and E scored >= T

    def totals(self, w=None):
        if w is None:
            return self.t.sum(), self.e.sum(), self.e_only.sum(), self.e_gt_t.sum()
        return w @ self.t, w @ self.e, w @ self.e_only, w @ self.e_gt_t


def paired_counts(pair_id, is_entrapment, score, accepted) -> EntrapmentCounts:
    """Per-pair counts for the paired estimator.

    Rows: the *target-side* entries of the tool's result (original targets and
    entrapment peptides; decoys removed), each with the key of its
    target-entrapment pair, its score, and whether it was accepted (q <= alpha).
    At most one row per (pair, side) is expected (take the best otherwise).
    """
    pair = as_codes(pair_id)
    ent = np.asarray(is_entrapment, dtype=bool)
    s = np.asarray(score, dtype=float)
    acc = np.asarray(accepted, dtype=bool)
    keys, inv = np.unique(pair, return_inverse=True)
    m = keys.shape[0]
    t_acc = np.zeros(m, dtype=np.int64)
    e_acc = np.zeros(m, dtype=np.int64)
    t_s = np.full(m, -np.inf)
    e_s = np.full(m, -np.inf)
    np.maximum.at(t_acc, inv[~ent], acc[~ent].astype(np.int64))
    np.maximum.at(e_acc, inv[ent], acc[ent].astype(np.int64))
    np.maximum.at(t_s, inv[~ent & acc], s[~ent & acc])
    np.maximum.at(e_s, inv[ent & acc], s[ent & acc])
    e_only = (e_acc == 1) & (t_acc == 0)
    e_gt_t = (e_acc == 1) & (t_acc == 1) & (e_s >= t_s)
    return EntrapmentCounts(key=keys, t=t_acc, e=e_acc, e_only=e_only.astype(np.int64),
                            e_gt_t=e_gt_t.astype(np.int64))


def combined_counts(item_id, is_entrapment, accepted) -> EntrapmentCounts:
    """Counts for the combined / lower-bound estimators; bootstrap unit = item."""
    ent = np.asarray(is_entrapment, dtype=bool)
    acc = np.asarray(accepted, dtype=bool)
    key = as_codes(item_id)
    z = np.zeros(key.shape[0], dtype=np.int64)
    return EntrapmentCounts(key=key, t=(acc & ~ent).astype(np.int64), e=(acc & ent).astype(np.int64),
                            e_only=z, e_gt_t=z)


# ------------------------------------------------------------------- protein level

def protein_group_is_target(groups: Sequence[Sequence], is_entrapment_protein: Mapping) -> np.ndarray:
    """Protein-level counting rule: a group with >= 1 original target protein is a
    target; otherwise (all members entrapment) it is an entrapment group."""
    return np.array([any(not is_entrapment_protein[p] for p in g) for g in groups], dtype=bool)


def protein_paired_counts(groups: Sequence[Sequence], score, accepted, is_entrapment_protein: Mapping,
                          partner: Mapping) -> EntrapmentCounts:
    """Paired estimator at protein-group level (FDRBench convention).

    An entrapment group is represented by its first member (as FDRBench); its
    paired original target protein is ``partner[first]``. The pair counts as
    "both accepted" when any accepted target group contains that partner; then
    the comparison is between the two group scores. Pair key = target protein.
    """
    score = np.asarray(score, float)
    accepted = np.asarray(accepted, bool)
    is_t = protein_group_is_target(groups, is_entrapment_protein)
    best_t = {}
    for g, s, a, it in zip(groups, score, accepted, is_t):
        if a and it:
            for p in g:
                if not is_entrapment_protein[p]:
                    best_t[p] = max(best_t.get(p, -np.inf), s)
    keys, t, e, e_only, e_gt = [], [], [], [], []
    seen_t = set()
    for g, s, a, it in zip(groups, score, accepted, is_t):
        if not a:
            continue
        if it:
            rep = next(p for p in g if not is_entrapment_protein[p])
            if rep in seen_t:
                continue
            seen_t.add(rep)
            keys.append(("T", rep)); t.append(1); e.append(0); e_only.append(0); e_gt.append(0)
        else:
            tp = partner[g[0]]
            keys.append(("E", tp)); t.append(0); e.append(1)
            if tp in best_t:
                e_only.append(0); e_gt.append(int(s >= best_t[tp]))
            else:
                e_only.append(1); e_gt.append(0)
    # bootstrap unit = the target protein of the pair
    k = as_codes(np.array([x[1] for x in keys], dtype=object).astype(str)) if keys else np.zeros(0, np.int64)
    return EntrapmentCounts(key=k, t=np.asarray(t, np.int64), e=np.asarray(e, np.int64),
                            e_only=np.asarray(e_only, np.int64), e_gt_t=np.asarray(e_gt, np.int64))


# --------------------------------------------------------------------- decision rule

@dataclass
class Replicate:
    """One (entrapment seed, decoy seed) replicate at one alpha."""

    counts: EntrapmentCounts
    entrap_seed: int | None = None
    decoy_seed: int | None = None


def _fdp(estimator: str, n_t, n_e, n_eo, n_egt, r):
    if estimator == "paired":
        return paired(n_t, n_e, n_eo, n_egt)
    if estimator == "combined":
        return combined(n_t, n_e, r)
    if estimator == "lower_bound":
        return lower_bound(n_t, n_e)
    raise ValueError(estimator)


def verdict(alpha: float, mean_fdp: float, ci_hi: float, lb_ci_lo: float, n_t: float) -> str:
    """Apply the v0.3 verdict table."""
    if alpha * n_t < 20:
        return "not-evaluable"
    if lb_ci_lo > alpha:
        return "violated"
    if mean_fdp <= alpha and ci_hi <= 1.5 * alpha:
        return "supported"
    if mean_fdp <= alpha:
        return "supported-noisy"
    return "inconclusive"


CI_CONVENTIONS = ("a", "b")


def evaluate_arm(replicates: Sequence[Replicate], alpha: float, estimator: str = "paired", r: float = 1.0,
                 n_boot: int = 2000, seed=11, z: float = 1.959963984540054, ci_convention: str = "a") -> dict:
    """Mean FDP over replicates, 95% CI, and the verdict.

    The Poisson(1) bootstrap draws one weight per pair key; replicates sharing a
    key (same original target across entrapment/decoy seeds) share the weight.

    ``ci_convention`` (DESIGN §6.0 A1.6): "a" = bootstrap variance + between-seed
    variance / J (the v0.3 rule, default); "b" = bootstrap variance only.  Both CIs
    are always returned (``ci95_a``, ``ci95_b``); ``ci95`` / the verdict use the
    requested one.
    """
    if ci_convention not in CI_CONVENTIONS:
        raise ValueError(ci_convention)
    J = len(replicates)
    if J == 0:
        raise ValueError("no replicates")
    allkeys = np.unique(np.concatenate([rep.counts.key for rep in replicates]))
    rng = make_rng(seed)
    fdp = np.zeros(J)
    lb = np.zeros(J)
    nT = np.zeros(J)
    boot_fdp = np.zeros((n_boot, J))
    boot_lb = np.zeros((n_boot, J))
    idxs = [np.searchsorted(allkeys, rep.counts.key) for rep in replicates]
    for j, rep in enumerate(replicates):
        c = rep.counts
        n_t, n_e, n_eo, n_egt = c.totals()
        nT[j] = n_t
        fdp[j] = _fdp(estimator, n_t, n_e, n_eo, n_egt, r)
        lb[j] = lower_bound(n_t, n_e)
    chunk = max(1, int(2e7 // max(allkeys.shape[0], 1)))
    for b0 in range(0, n_boot, chunk):
        b1 = min(n_boot, b0 + chunk)
        W = rng.poisson(1.0, size=(b1 - b0, allkeys.shape[0])).astype(np.float32)
        for j, rep in enumerate(replicates):
            c = rep.counts
            Wj = W[:, idxs[j]]
            bt, be = Wj @ c.t.astype(np.float32), Wj @ c.e.astype(np.float32)
            beo, begt = Wj @ c.e_only.astype(np.float32), Wj @ c.e_gt_t.astype(np.float32)
            boot_fdp[b0:b1, j] = _fdp(estimator, bt, be, beo, begt, r)
            boot_lb[b0:b1, j] = lower_bound(bt, be)
    mean = float(fdp.mean())
    mean_lb = float(lb.mean())
    between = float(fdp.var(ddof=1)) / J if J > 1 else 0.0
    between_lb = float(lb.var(ddof=1)) / J if J > 1 else 0.0
    var_boot = float(boot_fdp.mean(axis=1).var(ddof=1))
    var_boot_lb = float(boot_lb.mean(axis=1).var(ddof=1))
    se_a, se_b = float(np.sqrt(var_boot + between)), float(np.sqrt(var_boot))
    se_lb_a, se_lb_b = float(np.sqrt(var_boot_lb + between_lb)), float(np.sqrt(var_boot_lb))
    se, se_lb = (se_a, se_lb_a) if ci_convention == "a" else (se_b, se_lb_b)
    ci = (mean - z * se, mean + z * se)
    lb_ci = (mean_lb - z * se_lb, mean_lb + z * se_lb)
    v = verdict(alpha, mean, ci[1], lb_ci[0], float(nT.mean()))
    return {
        "alpha": alpha, "estimator": estimator, "r": r, "n_replicates": J,
        "mean_fdp": mean, "ci95": ci, "se": se, "fdp_per_replicate": fdp.tolist(),
        "mean_lb": mean_lb, "lb_ci95": lb_ci, "mean_n_t": float(nT.mean()),
        "verdict": v, "gated": any(abs(alpha - a) < 1e-12 for a in GATED_ALPHAS),
        "ci_convention": ci_convention,
        "ci95_a": (mean - z * se_a, mean + z * se_a), "ci95_b": (mean - z * se_b, mean + z * se_b),
        "lb_ci95_a": (mean_lb - z * se_lb_a, mean_lb + z * se_lb_a),
        "lb_ci95_b": (mean_lb - z * se_lb_b, mean_lb + z * se_lb_b),
        "var_boot": var_boot, "var_between_over_J": between, "lb_per_replicate": lb.tolist(),
    }


REQUIRED_ARMS = {
    "fixed": None,                      # P or C
    "static_low_capacity": ("P", "C"),
    "sequence_capacity": ("P", "C", "M"),
}


def combine_arms(arm_verdicts: Mapping[str, str], model_class: str) -> str:
    """Overall verdict across entrapment arms (DESIGN v0.3 §6.3 "Required arms").

    * any arm violated -> violated
    * fixed hand-crafted scores: supported if P or C supported (noisy likewise)
    * otherwise every required arm must be present; a missing/not-evaluable
      required arm -> not-evaluable; all supported -> supported; all in
      {supported, supported-noisy} -> supported-noisy; else inconclusive.
    """
    if any(v == "violated" for v in arm_verdicts.values()):
        return "violated"
    if model_class not in REQUIRED_ARMS:
        raise ValueError(model_class)
    req = REQUIRED_ARMS[model_class]
    if req is None:
        vs = [arm_verdicts.get(a) for a in ("P", "C") if a in arm_verdicts]
        if "supported" in vs:
            return "supported"
        if "supported-noisy" in vs:
            return "supported-noisy"
        if vs and all(v == "not-evaluable" for v in vs):
            return "not-evaluable"
        return "inconclusive" if vs else "not-evaluable"
    vs = [arm_verdicts.get(a) for a in req]
    if any(v is None or v == "not-evaluable" for v in vs):
        return "not-evaluable"
    if all(v == "supported" for v in vs):
        return "supported"
    if all(v in ("supported", "supported-noisy") for v in vs):
        return "supported-noisy"
    return "inconclusive"
