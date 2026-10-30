"""CF-MS transfers (DESIGN §4 L5 CF-MS; E10).

Policy (the default is mode='local', k=1):
  * identifications: run-specific (per-fraction) q-values; transfers are never identifications;
  * every transferred value carries provenance 'transfer_local' / 'transfer_global' and a transfer q-value;
  * fraction-local: a precursor can be transferred into fraction f only from donors identified in fractions
    f-k..f+k (k=1: adjacent fractions only) of the same series; mode='global' allows donors from any fraction
    (unconstrained MBR, the E10 comparator);
  * profile-agnostic: the transfer score is the M6 a-priori score (mass, RT, isotope and trace cosine) WITHOUT the
    intensity-consistency term, and nothing in score, window or q-value looks at the donor's (or any neighbour's)
    intensity; a profile-continuity prior would manufacture exactly the co-elution CF-MS scores;
  * no normalization: intensities are not centred or scaled per run anywhere (the M6 score's per-run median
    centring is not used either).
Transfer q-values are the M6 PIP-ECHO-style estimate (RT-shifted and m/z-shifted decoy peaks, decoy-peptide donors),
computed per acceptor run (vulcan.quant.lfq.transfer_qvalues).  Features already claimed by an MS2 ID of the
acceptor run are excluded from every window (target and decoy alike), as in M6.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from vulcan.quant import lfq as Q


@dataclass
class CFTransferParams:
    mode: str = "local"         # 'local' | 'global'
    k: int = 1                  # fraction window half-width for mode 'local'
    rt_k: float = 4.0
    rt_min_w: float = 0.05
    ppm_k: float = 3.0
    ppm_min: float = 2.0
    w_iso: float = 2.0
    w_tr: float = 1.0
    rt_shift_lo: float = 4.0
    rt_shift_hi: float = 8.0
    rt_shift_min: float = 0.5
    seed: int = 11
    exclude_claimed: bool = True


_NO_INTENSITY_SD = 1e30   # makes the intensity term of lfq._best_in_window exactly 0 (profile-agnostic)


def cf_transfers(tables: dict, ms1: dict, align: Q.Alignment, decoy_donor: np.ndarray, fraction: dict,
                 params: CFTransferParams | None = None, donor_keys=None) -> tuple[pd.DataFrame, dict]:
    """tables: run -> lfq.precursor_table (MS2-identified precursors of that run incl. decoy-peptide donors);
    fraction: run -> integer fraction index (one series).  Returns (candidates, info) with lfq.match_between_runs
    columns plus 'n_donor_fractions' and 'provenance'."""
    P = params or CFTransferParams()
    rng = np.random.default_rng(P.seed)
    allp = pd.concat([t.assign(run=run) for run, t in tables.items()])
    q = allp[allp["feature_id"] >= 0].copy()
    if donor_keys is not None:
        q = q[[k in donor_keys for k in zip(q["pf_id"], q["charge"])]]
    q["rt_cons"] = np.nan
    for run in tables:
        s = q["run"] == run
        q.loc[s, "rt_cons"] = align.to_consensus(run, q.loc[s, "rt_apex"].to_numpy())
    q["frac"] = np.array([fraction[r] for r in q["run"].astype(str)], np.int64)
    present = set(zip(allp["run"], allp["pf_id"], allp["charge"]))
    out, info = [], {"mode": P.mode, "k": P.k, "runs": {}}
    for run in sorted(tables, key=lambda r: fraction[r]):
        f = fraction[run]
        src = q if P.mode == "global" else q[(q["frac"] - f).abs().between(1, P.k)]
        if not len(src):
            continue
        D = src.groupby(["pf_id", "charge"]).agg(rt_cons=("rt_cons", "median"), calcmass=("calcmass", "first"),
                                                 n_donor_fractions=("run", "nunique")).reset_index()
        D = D[[(run, a, b) not in present for a, b in zip(D["pf_id"], D["charge"])]].reset_index(drop=True)
        if not len(D):
            continue
        D["is_decoy_donor"] = np.isin(D["pf_id"].to_numpy(), decoy_donor)
        r = ms1[run]
        F = r.feat
        fkey_raw = F["charge"].to_numpy().astype(np.float64) * 1e5 + F["neutral_mass"].to_numpy()
        o = np.argsort(fkey_raw, kind="stable")
        fok = np.ones(len(F), np.bool_)
        if P.exclude_claimed:
            cl = tables[run]["feature_id"].to_numpy()
            fok[cl[cl >= 0]] = False
        args = (fkey_raw[o], F["rt_apex"].to_numpy().astype(np.float64), F["neutral_mass"].to_numpy(),
                F["iso_cos"].to_numpy().astype(np.float64), F["trace_cos"].to_numpy().astype(np.float64),
                np.zeros(len(F)), fok, o)
        t_hat = align.to_run(run, D["rt_cons"].to_numpy())
        sd_rt = np.maximum(align.local_sd(run, D["rt_cons"].to_numpy()), 1e-3)
        W = np.maximum(P.rt_k * sd_rt, P.rt_min_w)
        sd_eff = W / P.rt_k
        mass = D["calcmass"].to_numpy() * (1 + r.ppm_offset * 1e-6)
        tol = np.full(len(D), max(P.ppm_k * r.ppm_sd, P.ppm_min))
        sdp = np.full(len(D), max(r.ppm_sd, 0.5))
        z = D["charge"].to_numpy().astype(np.float64)
        mu = np.zeros(len(D))
        St, At = Q._best_in_window(*args, z, mass, t_hat, W, sd_eff, tol, sdp, mu, _NO_INTENSITY_SD, P.w_iso, P.w_tr)
        mag = np.maximum(rng.uniform(P.rt_shift_lo, P.rt_shift_hi, len(D)) * W, P.rt_shift_min)
        sign = np.where(rng.random(len(D)) < 0.5, -1.0, 1.0)
        lo_rt, hi_rt = r.rt_range
        t_d = t_hat + sign * mag
        bad = (t_d < lo_rt + W) | (t_d > hi_rt - W)
        t_d = np.where(bad, t_hat - sign * mag, t_d)
        bad2 = (t_d < lo_rt + W) | (t_d > hi_rt - W)
        Sr, _ = Q._best_in_window(*args, z, mass, t_d, W, sd_eff, tol, sdp, mu, _NO_INTENSITY_SD, P.w_iso, P.w_tr)
        Sr[bad2] = -np.inf
        Sm, _ = Q._best_in_window(*args, z, mass + Q.MZ_DECOY_SHIFT_DA, t_hat, W, sd_eff, tol, sdp, mu,
                                  _NO_INTENSITY_SD, P.w_iso, P.w_tr)
        T = pd.DataFrame({"run": run, "fraction": f, "pf_id": D["pf_id"], "charge": D["charge"],
                          "calcmass": D["calcmass"], "n_donor_fractions": D["n_donor_fractions"],
                          "is_decoy_donor": D["is_decoy_donor"], "rt_pred": t_hat, "rt_window": W,
                          "St": St, "Sd_rt": Sr, "Sd_mz": Sm, "feature_id": At})
        ok = At >= 0
        for col, s_ in (("intensity", "intensity_area"), ("rt_apex", "rt_apex"), ("feat_mass", "neutral_mass")):
            v = np.full(len(T), np.nan)
            v[ok] = F[s_].to_numpy()[At[ok]]
            T[col] = v
        qa, qb = Q.transfer_qvalues(T, P.seed)
        T["transfer_q"], T["transfer_q_peak"] = qa, qb
        out.append(T)
        info["runs"][run] = {"fraction": f, "n_candidates": int(len(T)), "n_with_target_peak": int(ok.sum())}
    T = pd.concat(out, ignore_index=True)
    T["provenance"] = "transfer_" + P.mode
    return T, info
