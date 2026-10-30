"""L3 localization operator (DESIGN.md §4 L3): exhaustive site configurations, decoy-amino-acid candidates,
site-determining ions with neutral-loss awareness, scoring blind to residue identity given masses.

For each PSM (a peptidoform found by the identification search) the operator re-places its k localizable
modifications (phospho by default) over a set of candidate positions and scores **every** configuration
(C(n, k) of them; no pruning). Candidate positions are supplied by the caller: S/T/Y plus the decoy residue
(Ala primary, Leu/Ile secondary; DESIGN §6.5). The decoy residue enters only here, never in identification.

Residue blindness
-----------------
The scoring kernel `score_configs` sees only numbers: per-residue masses (static and fixed variable mods
included, the relocalized mod excluded), an N-terminal delta, candidate positions, k and the spectrum. It never
receives residue letters, so a configuration placing the mod on an Ala is scored with exactly the same ions and
the same neutral loss as one placing it on a Ser at the same position with the same masses.
`tests/ptm/test_localize.py` checks this (decoy-AA G0 cases).

Ions (HCD): b and y, fragment charges 1..min(z-1, max_frag_charge); for every fragment that contains >= 1 of the
relocalized mods also the neutral-loss ion (-H3PO4 by default; one loss), for every acceptor alike.

Score: peaks are reduced to the top `peaks_per_window` per `window_th` Th window; each kept peak weighs
w = ln(1 + I / median I of the kept peaks) (default, `weight="intensity"`) or, with `weight="rank"`, a peak of rank r in its
window weighs w = ln((N + 1) / r). S(config) = sum of w over the distinct peaks matched (within `frag_ppm`) by any
of the config's ions. Ions common to all configurations add the same constant, so differences between
configurations come only from site-determining ions.

Outputs per PSM: best configuration (ties: lowest configuration index, i.e. lexicographic positions, which depends
on positions only), its score, and per site of the best configuration the site score
    delta_site = S(best) - max{ S(c) : c does not contain that site }   (+inf if every configuration contains it).
Higher delta = more confident. FLR̂ (vulcan.stats.flr) is computed over sites ranked by delta_site.
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

import numpy as np

PROTON = 1.007276467
H2O = 18.010564684
H3PO4 = 97.976895573
PHOSPHO = 79.966331

# monoisotopic residue masses (Unimod), index = ASCII - 'A'; used only by `residue_masses` (the caller side)
RESIDUE = {"A": 71.037113805, "C": 103.009184505, "D": 115.026943065, "E": 129.042593135, "F": 147.068413945,
           "G": 57.021463721, "H": 137.058911875, "I": 113.084064015, "K": 128.094963050, "L": 113.084064015,
           "M": 131.040484645, "N": 114.042927470, "P": 97.052763875, "Q": 128.058577540, "R": 156.101111050,
           "S": 87.032028435, "T": 101.047678505, "V": 99.068413945, "W": 186.079312980, "Y": 163.063328575}


@dataclass
class LocParams:
    mod_mass: float = PHOSPHO
    nl_mass: float = H3PO4            # neutral loss of the relocalized mod (0 = none)
    frag_ppm: float = 20.0
    peaks_per_window: int = 12         # dev grid (NOTES §E8): {4,8,12,20,150} x {rank,intensity}
    window_th: float = 100.0
    max_frag_charge: int = 2
    min_mz: float = 100.0
    max_mz: float = 2000.0
    weight: str = "intensity"         # "rank": ln((N+1)/r) of the peak's intensity rank r in its window;
                                      # "intensity": ln(1 + I / median I of the kept peaks)
    max_configs: int = 20000          # per PSM; above this the PSM is flagged and skipped (never silently pruned)
    batch_ions: int = 60_000_000


def residue_masses(seq: str, static: dict, fixed_var: dict[int, float] | None = None,
                   table: dict | None = None) -> np.ndarray:
    """Caller-side helper: residue masses incl. static mods and position-fixed variable mods (e.g. M-ox)."""
    t = RESIDUE if table is None else table
    m = np.array([t[c] + static.get(c, 0.0) for c in seq], np.float64)
    for i, d in (fixed_var or {}).items():
        m[i] += d
    return m


def candidate_positions(seq: str, acceptors: str) -> np.ndarray:
    """Caller-side helper: positions whose residue is in `acceptors` (I is read as L when 'L' is given)."""
    s = seq.replace("I", "L") if "L" in acceptors else seq
    return np.array([i for i, c in enumerate(s) if c in acceptors], np.int64)


# ------------------------------------------------------------------------------------------------ peaks
def prepare_peaks(mz_list, int_list, p: LocParams):
    """Top-N per window, rank weights. Returns (offsets, mz f64, weight f64) concatenated over spectra."""
    offs = [0]
    mzs, ws = [], []
    for mz, it in zip(mz_list, int_list):
        mz = np.asarray(mz, np.float64)
        it = np.asarray(it, np.float64)
        ok = (mz >= p.min_mz) & (mz <= p.max_mz) & np.isfinite(mz) & (it > 0)
        mz, it = mz[ok], it[ok]
        win = np.floor(mz / p.window_th).astype(np.int64)
        o = np.lexsort((-mz, -it, win))          # window, intensity desc, m/z desc (deterministic)
        w_s = win[o]
        first = np.searchsorted(w_s, w_s, side="left")
        rank = np.arange(len(o)) - first + 1
        keep = rank <= p.peaks_per_window
        sel = o[keep]
        if p.weight == "rank":
            wt = np.log((p.peaks_per_window + 1.0) / rank[keep])
        else:
            ik = it[sel]
            wt = np.log1p(ik / max(np.median(ik), 1e-12)) if len(ik) else np.zeros(0)
        srt = np.argsort(mz[sel], kind="stable")
        mzs.append(mz[sel][srt])
        ws.append(wt[srt])
        offs.append(offs[-1] + len(sel))
    return (np.array(offs, np.int64), np.concatenate(mzs) if mzs else np.zeros(0),
            np.concatenate(ws) if ws else np.zeros(0))


# ------------------------------------------------------------------------------------------------ configurations
def enumerate_configs(n_cand: int, k: int) -> np.ndarray:
    """All k-subsets of range(n_cand) in lexicographic order, as an (n_cfg, k) int array."""
    if k == 0 or k > n_cand:
        return np.zeros((1 if k == 0 else 0, k), np.int64)
    return np.array(list(itertools.combinations(range(n_cand), k)), np.int64).reshape(-1, k)


# ------------------------------------------------------------------------------------------------ scoring
def score_configs(res_mass: list[np.ndarray], nterm: np.ndarray, cand_pos: list[np.ndarray], k: np.ndarray,
                  charge: np.ndarray, spec_slot: np.ndarray, peaks, p: LocParams, device: str = "cuda"):
    """Score every configuration of every PSM. Letter-free by construction (see module docstring).

    res_mass[i]: residue masses of PSM i (relocalized mod excluded); nterm[i]: N-terminal delta; cand_pos[i]:
    candidate positions; k[i]: number of relocalized mods; charge[i]: precursor charge; spec_slot[i]: index of the
    spectrum in `peaks` (from prepare_peaks).
    Returns per PSM: dict(cfg = (n_cfg, k) positions, score = (n_cfg,) f64) or None if skipped (too many configs).
    """
    import torch
    dev = torch.device(device if torch.cuda.is_available() else "cpu")
    p_off, p_mz, p_w = peaks
    # global sorted key: slot * 1e4 + m/z  (m/z < 1e4; f64 keeps ~1e-12 relative)
    slot_of_peak = np.repeat(np.arange(len(p_off) - 1), np.diff(p_off))
    key = torch.from_numpy(slot_of_peak * 1.0e4 + p_mz).to(dev)
    wt = torch.from_numpy(p_w).to(dev)
    n_pk = len(p_mz)
    out: list = [None] * len(res_mass)

    # build per-PSM configuration lists, then process PSMs in batches bounded by the ion count
    cfgs = []
    for i in range(len(res_mass)):
        n = len(cand_pos[i])
        kk = int(k[i])
        if kk == 0 or kk > n:
            cfgs.append(None)
            continue
        nc = math.comb(n, kk)
        if nc > p.max_configs:
            cfgs.append(None)
            continue
        cfgs.append(np.asarray(cand_pos[i], np.int64)[enumerate_configs(n, kk)])

    def ions_per(i):
        L = len(res_mass[i])
        zmax = max(1, min(int(charge[i]) - 1, p.max_frag_charge))
        return len(cfgs[i]) * 2 * (L - 1) * zmax * (2 if p.nl_mass else 1)

    order = [i for i in range(len(res_mass)) if cfgs[i] is not None]
    b0 = 0
    while b0 < len(order):
        tot, b1 = 0, b0
        while b1 < len(order) and (tot == 0 or tot + ions_per(order[b1]) <= p.batch_ions):
            tot += ions_per(order[b1])
            b1 += 1
        batch = order[b0:b1]
        b0 = b1
        # ---- flatten (psm, config) rows and the per-row modified residue masses (host, vectorized per PSM)
        rows_psm, prefix_rows, nmod_rows, tot_rows, zmax_rows, Ls = [], [], [], [], [], []
        Lmax = max(len(res_mass[i]) for i in batch)
        for i in batch:
            m = np.asarray(res_mass[i], np.float64)
            L = len(m)
            cf = cfgs[i]
            n_cfg = len(cf)
            onehot = np.zeros((n_cfg, L), np.float64)
            np.put_along_axis(onehot, cf, 1.0, axis=1)
            mm = m[None, :] + onehot * p.mod_mass
            pre = np.zeros((n_cfg, Lmax), np.float64)
            pre[:, :L - 1] = np.cumsum(mm, axis=1)[:, :L - 1] + nterm[i]        # b-ion neutral masses (b1..b_{L-1})
            cnt = np.zeros((n_cfg, Lmax), np.float64)
            cnt[:, :L - 1] = np.cumsum(onehot, axis=1)[:, :L - 1]             # mods inside each b fragment
            rows_psm.append(np.full(n_cfg, i))
            prefix_rows.append(pre)
            nmod_rows.append(cnt)
            tot_rows.append(np.full(n_cfg, m.sum() + p.mod_mass * cf.shape[1] + nterm[i] + H2O))
            zmax_rows.append(np.full(n_cfg, max(1, min(int(charge[i]) - 1, p.max_frag_charge))))
            Ls.append(np.full(n_cfg, L))
        rp = np.concatenate(rows_psm)
        B = torch.from_numpy(np.concatenate(prefix_rows)).to(dev)            # (R, Lmax) b neutral
        Cn = torch.from_numpy(np.concatenate(nmod_rows)).to(dev)             # (R, Lmax) mods in b fragment
        T = torch.from_numpy(np.concatenate(tot_rows)).to(dev)[:, None]
        Kc = torch.from_numpy(np.array([cfgs[i].shape[1] for i in rp], np.float64)).to(dev)[:, None]
        Lr = torch.from_numpy(np.concatenate(Ls)).to(dev)[:, None]
        Zr = torch.from_numpy(np.concatenate(zmax_rows)).to(dev)[:, None]
        R = B.shape[0]
        j = torch.arange(Lmax, device=dev)[None, :]
        valid = j < (Lr - 1)
        Y = T - B                       # y neutral complementary to b_j (total includes H2O); mods inside y = K - in b
        Cy = Kc - Cn
        ions, has = [], []
        for base, c in ((B, Cn), (Y, Cy)):
            ions.append(base)
            has.append(valid)
            if p.nl_mass:
                ions.append(base - p.nl_mass)
                has.append(valid & (c > 0))
        ion = torch.stack(ions, 0)                                            # (V, R, Lmax) neutral masses
        ok = torch.stack(has, 0)
        slots = torch.from_numpy(np.asarray(spec_slot)[rp].astype(np.float64)).to(dev)[None, :, None]
        row_id = torch.arange(R, device=dev)[None, :, None].expand_as(ion)
        keys_all, rows_all = [], []
        for z in range(1, p.max_frag_charge + 1):
            okz = ok & (Zr[None] >= z)
            mz = (ion + z * PROTON) / z
            okz = okz & (mz >= p.min_mz) & (mz <= p.max_mz)
            kq = (slots.expand_as(mz) * 1.0e4 + mz)[okz]
            keys_all.append(kq)
            rows_all.append(row_id[okz])
        q = torch.cat(keys_all)
        qr = torch.cat(rows_all)
        mzq = q - torch.floor(q / 1.0e4) * 1.0e4
        tol = mzq * p.frag_ppm * 1e-6
        lo = torch.searchsorted(key, q - tol, side="left")
        hit = lo < n_pk
        lo_c = lo.clamp(max=max(n_pk - 1, 0))
        hit = hit & (key[lo_c] <= q + tol)
        # nearest peak within tolerance: also look at lo+1 (peaks are sparse after top-N; pick max weight)
        lo2 = (lo_c + 1).clamp(max=max(n_pk - 1, 0))
        hit2 = (lo2 != lo_c) & (key[lo2] <= q + tol) & (key[lo2] >= q - tol)
        pk = torch.where(hit2 & (~hit | (wt[lo2] > wt[lo_c])), lo2, lo_c)
        anyhit = hit | hit2
        pk, qr = pk[anyhit], qr[anyhit]
        # distinct (row, peak) pairs
        u = torch.unique(qr.to(torch.int64) * (n_pk + 1) + pk)
        ur = u // (n_pk + 1)
        up = u % (n_pk + 1)
        sc = torch.zeros(R, dtype=torch.float64, device=dev)
        sc.index_add_(0, ur, wt[up])
        sc = sc.cpu().numpy()
        # split back per PSM (rows of one PSM are contiguous, in config order)
        st = 0
        for i in batch:
            n_cfg = len(cfgs[i])
            out[i] = {"cfg": cfgs[i], "score": sc[st:st + n_cfg]}
            st += n_cfg
    return out


def summarize(scored, cand_pos):
    """Best configuration and per-site delta scores. Returns list of dicts (None for skipped PSMs):
    best_cfg (positions), best_score, second_score (best config with a different site set), site_delta (per site
    of best_cfg)."""
    res = []
    for s in scored:
        if s is None:
            res.append(None)
            continue
        sc, cf = s["score"], s["cfg"]
        b = int(np.argmax(sc))                         # first max = lowest config index (position-only tie rule)
        best = cf[b]
        others = np.delete(sc, b)
        second = float(others.max()) if len(others) else -np.inf
        deltas = []
        for site in best:
            without = ~np.any(cf == site, axis=1)
            deltas.append(float(sc[b] - sc[without].max()) if without.any() else float("inf"))
        res.append({"best_cfg": best, "best_score": float(sc[b]), "second_score": second,
                    "site_delta": np.array(deltas), "n_cfg": len(cf)})
    return res
