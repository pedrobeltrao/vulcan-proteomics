"""Per-run mass calibration as offsets (DESIGN §4 L0.5): ppm error = f(rt) + g(m/z), never baked into the index.

The model is additive and piecewise linear: f and g are given by knots (equal-count bins, robust medians, 3-point
running-median smoothing) and fitted by backfitting. Outside the knot range the end values hold. A component is only
kept if it lowers the robust spread (MAD) of the anchors by >= 2% (otherwise it is set to zero), which keeps
the model from fitting noise on small anchor sets.

    corrected m/z = observed m/z / (1 + ppm(rt, mz) * 1e-6)

Empirical tolerances (both reported; the search uses the "core" rule):
  * core (default): iterated sigma clipping (median, 1.4826*MAD, clip at 5 SD, 5 rounds) gives the location mu and
    spread sd of the error core; tolerance = mu +- k*sd (k = 4), floored at +-floor. The fraction of anchors inside is
    reported. Anchor errors on D0 have a Gaussian-like core plus a flat tail of wide-window matches and a deamidation-
    vs-isotope shoulder, so quantile rules over all anchors measure the contaminants, not the instrument.
  * coverage: matched errors inside a wide window contain random matches. The background is estimated as a uniform
    density from the outer band (0.6..1 of the window) and subtracted; the tolerance is the narrowest window around
    the median that holds `coverage` (99%) of the background-subtracted signal plateau (median over the outer 70%).
    Noisy when the background is a sizeable share of the matches; reported as a diagnostic.
  * quantile (legacy, precursor only): [q0.5%, q99.5%] after trimming > 6 robust SD, padded 10%.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field

import numpy as np


def _mad_sd(x):
    x = np.asarray(x, float)
    if len(x) == 0:
        return float("nan")
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


@dataclass
class Curve:
    """Piecewise-linear 1-D offset curve (ppm), constant beyond the knots. Empty knots = zero."""
    knots: list = field(default_factory=list)
    values: list = field(default_factory=list)

    def __call__(self, x):
        x = np.asarray(x, float)
        if not self.knots:
            return np.zeros_like(x)
        return np.interp(x, np.asarray(self.knots), np.asarray(self.values))

    @staticmethod
    def fit(x, y, n_min_per_bin: int = 40, max_bins: int = 24, min_bins: int = 2) -> "Curve":
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        n = len(x)
        nb = int(np.clip(n // n_min_per_bin, 0, max_bins))
        if nb < min_bins:
            return Curve([float(np.median(x))] if n else [], [float(np.median(y))] if n else [])
        o = np.argsort(x, kind="stable")
        xs, ys = x[o], y[o]
        edges = np.linspace(0, n, nb + 1).astype(int)
        kx = np.array([np.median(xs[a:b]) for a, b in zip(edges[:-1], edges[1:])])
        ky = np.array([np.median(ys[a:b]) for a, b in zip(edges[:-1], edges[1:])])
        if nb >= 3:   # 3-point running median (end points kept)
            sm = ky.copy()
            sm[1:-1] = np.median(np.stack([ky[:-2], ky[1:-1], ky[2:]]), axis=0)
            ky = sm
        return Curve(kx.tolist(), ky.tolist())


@dataclass
class OffsetModel:
    """ppm(rt, mz) = rt_curve(rt) + mz_curve(mz)."""
    rt_curve: Curve = field(default_factory=Curve)
    mz_curve: Curve = field(default_factory=Curve)
    n_anchors: int = 0
    mad_before: float = float("nan")
    mad_after: float = float("nan")
    median_before: float = float("nan")
    used: dict = field(default_factory=lambda: {"rt": False, "mz": False})

    def __call__(self, rt, mz):
        return self.rt_curve(rt) + self.mz_curve(mz)

    @staticmethod
    def fit(rt, mz, err_ppm, n_iter: int = 3, min_gain: float = 0.02, **kw) -> "OffsetModel":
        rt, mz, e = (np.asarray(a, float) for a in (rt, mz, err_ppm))
        m = OffsetModel(n_anchors=len(e), mad_before=_mad_sd(e), median_before=float(np.median(e)) if len(e) else 0.0)
        if len(e) == 0:
            return m
        # global offset first (always used), then backfit rt and mz
        c0 = float(np.median(e))
        f_rt = Curve([float(np.median(rt))], [c0])
        f_mz = Curve()
        for _ in range(n_iter):
            f_rt = Curve.fit(rt, e - f_mz(mz), **kw)
            r = e - f_rt(rt)
            f_mz = Curve.fit(mz, r, **kw)
            # center the mz curve so the rt curve carries the global offset
            if f_mz.values:
                mid = float(np.median(f_mz(mz)))
                f_mz = Curve(f_mz.knots, (np.asarray(f_mz.values) - mid).tolist())
        const = Curve([float(np.median(rt))], [c0])
        mad_c = _mad_sd(e - c0)
        mad_rt = _mad_sd(e - f_rt(rt))
        use_rt = mad_rt < (1 - min_gain) * mad_c
        base = f_rt if use_rt else const
        mad_b = _mad_sd(e - base(rt))
        mad_full = _mad_sd(e - base(rt) - f_mz(mz))
        use_mz = bool(f_mz.values) and mad_full < (1 - min_gain) * mad_b
        m.rt_curve = base
        m.mz_curve = f_mz if use_mz else Curve()
        m.used = {"rt": bool(use_rt), "mz": bool(use_mz)}
        m.mad_after = _mad_sd(e - m(rt, mz))
        return m

    def to_dict(self):
        return {"rt_knots": self.rt_curve.knots, "rt_ppm": self.rt_curve.values, "mz_knots": self.mz_curve.knots,
                "mz_ppm": self.mz_curve.values, "n_anchors": self.n_anchors, "mad_sd_before_ppm": self.mad_before,
                "mad_sd_after_ppm": self.mad_after, "median_before_ppm": self.median_before, "used": self.used}

    @staticmethod
    def from_dict(d) -> "OffsetModel":
        return OffsetModel(Curve(list(d["rt_knots"]), list(d["rt_ppm"])), Curve(list(d["mz_knots"]), list(d["mz_ppm"])),
                           d.get("n_anchors", 0), d.get("mad_sd_before_ppm", float("nan")),
                           d.get("mad_sd_after_ppm", float("nan")), d.get("median_before_ppm", float("nan")),
                           d.get("used", {}))


def precursor_tolerance(err_ppm, q=(0.005, 0.995), pad: float = 0.10, floor_ppm: float = 2.0, trim_sd: float = 6.0):
    e = np.asarray(err_ppm, float)
    if len(e) < 20:
        return None
    med, sd = np.median(e), _mad_sd(e)
    e = e[np.abs(e - med) <= trim_sd * max(sd, 1e-3)]
    lo, hi = np.quantile(e, q)
    w = hi - lo
    lo, hi = lo - pad * w / 2, hi + pad * w / 2
    return [float(min(lo, -floor_ppm)), float(max(hi, floor_ppm))]


def core_tolerance(err, k: float = 4.0, clip: float = 5.0, iters: int = 5, floor: float = 0.0):
    """mu +- k*sd of the sigma-clipped error core. Returns ([lo, hi], info) or (None, info) if < 20 points."""
    e = np.asarray(err, float)
    e = e[np.isfinite(e)]
    if len(e) < 20:
        return None, {"n": int(len(e))}
    m = np.ones(len(e), bool)
    mu, sd = float(np.median(e)), _mad_sd(e)
    for _ in range(iters):
        m = np.abs(e - mu) <= clip * max(sd, 1e-6)
        mu, sd = float(np.median(e[m])), _mad_sd(e[m])
    t = max(k * sd, floor)
    inside = float(np.mean(np.abs(e - mu) <= t))
    return [mu - t, mu + t], {"n": int(len(e)), "mu": mu, "sd": sd, "k": k, "frac_inside": inside,
                              "frac_core": float(m.mean())}


def fragment_tolerance(err, window: float, coverage: float = 0.99, band=(0.6, 1.0), floor: float = 0.0, bins: int = 400):
    """Background-subtracted tolerance (same unit as err/window). Returns ([lo, hi], info)."""
    e = np.asarray(err, float)
    e = e[np.abs(e) <= window]
    if len(e) < 100:
        return None, {"n": int(len(e))}
    med = float(np.median(e[np.abs(e) <= 0.5 * window]))
    outer = (np.abs(e) >= band[0] * window) & (np.abs(e) <= band[1] * window)
    dens = outer.sum() / (2 * (band[1] - band[0]) * window)          # background matches per unit
    d = np.sort(np.abs(e - med))
    grid = np.linspace(0, window - abs(med), bins + 1)[1:]
    cnt = np.searchsorted(d, grid, side="right")
    sig = cnt - dens * 2 * grid
    total = float(np.median(sig[grid >= 0.3 * window]))     # plateau (a max would chase background noise)
    k = int(np.argmax(sig >= coverage * total))
    t = float(max(grid[k], floor))
    info = {"n": int(len(e)), "median": med, "background_per_unit": float(dens),
            "signal_total": float(total), "background_frac_in_window": float(dens * 2 * window / len(e))}
    return [med - t, med + t], info


@dataclass
class Calibration:
    """Per-run calibration (JSON). Offsets in ppm; apply with correct_precursor / correct_fragment."""
    run_id: int = 0
    run_name: str = ""
    precursor: OffsetModel = field(default_factory=OffsetModel)
    fragment: OffsetModel = field(default_factory=OffsetModel)
    fragment_unit: str = "ppm"
    tolerances: dict = field(default_factory=dict)
    extra: dict = field(default_factory=dict)

    def correct_precursor(self, mz, rt):
        mz = np.asarray(mz, float)
        return mz / (1.0 + self.precursor(rt, mz) * 1e-6)

    def correct_fragment(self, mz, rt):
        """rt broadcast to mz (pass the per-peak rt)."""
        mz = np.asarray(mz, float)
        if self.fragment_unit != "ppm":
            return mz
        return mz / (1.0 + self.fragment(rt, mz) * 1e-6)

    def summary(self) -> dict:
        return {"run_id": self.run_id, "precursor": self.precursor.to_dict(), "fragment": self.fragment.to_dict(),
                "fragment_unit": self.fragment_unit, "tolerances": self.tolerances}

    def to_json(self, path: str):
        d = self.summary()
        d.update(run_name=self.run_name, extra=self.extra, format="vulcan.calib 1")
        with open(path, "w") as fh:
            json.dump(d, fh, indent=1, default=float)

    @staticmethod
    def from_json(path: str) -> "Calibration":
        d = json.load(open(path))
        return Calibration(d["run_id"], d.get("run_name", ""), OffsetModel.from_dict(d["precursor"]),
                           OffsetModel.from_dict(d["fragment"]), d.get("fragment_unit", "ppm"), d.get("tolerances", {}),
                           d.get("extra", {}))
