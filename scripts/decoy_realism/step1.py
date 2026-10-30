#!/usr/bin/env python3
"""3g step 1 + H3 strata (PLAN.md): homology dose-response of the arm-C E-vs-dec(E) coin on the E9 tail tables, the
within-pair ladder-sharing model, and the coin by length / charge / missed / mass decile / n_var / decoy kind.
    conda_envs/pe/bin/python scripts/decoy_realism/step1.py
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import chi2_contingency  # noqa: E402
import common as C  # noqa: E402


def il(s):
    return s.replace("I", "L") if isinstance(s, str) else s


def levels(H):
    l1 = H.ham <= 2
    l2 = H.k6
    l3 = H.k5
    l4 = H.sw_evalue <= 10
    alone = {"L1_ham2": l1, "L2_k6": l2, "L3_k5": l3, "L4_sw_e10": l4}
    l1, l2, l3, l4 = (np.asarray(x) for x in (l1, l2, l3, l4))
    alone = {"L1_ham2": l1, "L2_k6": l2, "L3_k5": l3, "L4_sw_e10": l4}
    cum = {"L0": np.zeros(len(H), bool), "L1": l1, "L2": l1 | l2, "L3": l1 | l2 | l3, "L4": l1 | l2 | l3 | l4}
    return alone, cum


def strata_table(tab, col, bins=None):
    v = tab[col].values
    if bins is not None:
        v = np.digitize(v, bins)
    out, ct = {}, []
    tot_exc = tab.E_won.sum() - len(tab) / 2
    for k in np.unique(v):
        m = v == k
        w = tab.E_won.values[m]
        c = C.coin(w)
        c["share_pairs"] = float(m.mean())
        c["share_excess"] = float((w.sum() - m.sum() / 2) / tot_exc) if tot_exc else None
        out[str(k)] = c
        ct.append([w.sum(), m.sum() - w.sum()])
    ct = np.array(ct)
    ok = ct.sum(1) >= 20
    p = float(chi2_contingency(ct[ok])[1]) if ok.sum() >= 2 else None
    flag = any(o["share_excess"] is not None and o["share_excess"] > 0.5 and o["share_pairs"] < 0.25 and o["n"] >= 100
               for o in out.values())
    return {"strata": out, "heterogeneity_p": p, "H3_flag": bool(flag)}


def logit_fit(y, X, n_boot=500):
    from sklearn.linear_model import LogisticRegression
    def fit(yy, XX):
        m = LogisticRegression(penalty=None, max_iter=1000).fit(XX, yy)
        return np.concatenate([[m.intercept_[0]], m.coef_[0]])
    b = fit(y, X)
    rng = np.random.default_rng(C.SEED)
    bs = np.array([fit(y[i], X[i]) for i in (rng.integers(0, len(y), len(y)) for _ in range(n_boot))])
    sig = lambda z: 1 / (1 + np.exp(-z))  # noqa: E731
    return {"coef": b.tolist(), "coef_ci95": np.quantile(bs, [.025, .975], axis=0).T.tolist(),
            "coin_at_delta0": float(sig(b[0])), "coin_at_delta0_ci95": [float(sig(np.quantile(bs[:, 0], .025))),
                                                                       float(sig(np.quantile(bs[:, 0], .975)))]}


def main():
    H = pd.read_parquet(os.path.join(C.SCR, "homology_base.parquet")).set_index("seq")
    Dd = pd.read_parquet(os.path.join(C.SCR, "homology_dec.parquet")).set_index("seq")
    rnd = pd.read_parquet(os.path.join(C.SCR, "random_base.parquet")).seq.map(il)
    res = {}
    Hr = H.loc[rnd.values]
    alone_r, cum_r = levels(Hr)
    res["background_removal_random20k"] = {**{k: float(np.mean(v)) for k, v in alone_r.items()},
                                           **{f"cum_{k}": float(np.mean(v)) for k, v in cum_r.items()},
                                           "sw_evalue_le_1000": float((Hr.sw_evalue <= 1000).mean()),
                                           "hy_mean": float(Hr.hy.mean()), "hb_mean": float(Hr.hb.mean())}
    for col in ("score", "hyperscore"):
        tab = pd.read_parquet(os.path.join(C.SCR, f"e9tail_C_{col}.parquet"))
        tab["b"] = tab.base_seq.map(il)
        Ht = H.loc[tab.b.values].reset_index(drop=True)
        alone, cum = levels(Ht)
        tab["lstr"] = C.len_stratum(tab.length.values)
        r = {}
        for thr, msk in (("0.05", np.ones(len(tab), bool)), ("0.01", tab.above1.values)):
            w0 = tab.E_won.values[msk]
            ls0 = tab.lstr.values[msk]
            rr = {"cumulative": {}, "alone": {}}
            for kind, dct in (("cumulative", cum), ("alone", alone)):
                for k, rem in dct.items():
                    keep = ~np.asarray(rem)[msk]
                    w = w0[keep]
                    e = C.coin(w)
                    e["frac_pairs_kept"] = float(keep.mean())
                    e["evaluable"] = bool(len(w) >= 400)
                    if len(w):
                        e["length_std"] = C.std_coin(w, ls0[keep], ls0)
                        e["retained_excess"] = C.boot_ratio(w0, w)
                        # length-standardised retained excess (point estimate)
                        e["retained_excess_std"] = float((e["length_std"]["coin_std"] - .5) / (w0.mean() - .5))
                    rr[kind][k] = e
            # continuous: SW normalised-score quintiles, ladder sharing
            q = np.quantile(Ht.sw_norm.values[msk], [.2, .4, .6, .8])
            rr["by_sw_norm_quintile"] = strata_table(pd.DataFrame({"E_won": w0, "x": Ht.sw_norm.values[msk]}), "x", q)
            rr["by_sw_norm_quintile"]["edges"] = q.tolist()
            # within-pair ladder model
            hyE, hbE = Ht.hy.values[msk], Ht.hb.values[msk]
            ds = tab.dec_seq.map(il).values[msk]
            okd = pd.notna(ds)
            hyD = np.array([Dd.hy.get(s, np.nan) if s is not None else np.nan for s in ds], float)
            hbD = np.array([Dd.hb.get(s, np.nan) if s is not None else np.nan for s in ds], float)
            ok = okd & np.isfinite(hyD)
            dy, db = hyE[ok] - hyD[ok], hbE[ok] - hbD[ok]
            rr["ladder"] = {"n": int(ok.sum()), "mean_hy_E": float(hyE[ok].mean()), "mean_hy_dec": float(hyD[ok].mean()),
                            "mean_hb_E": float(hbE[ok].mean()), "mean_hb_dec": float(hbD[ok].mean()),
                            "model_dy_db": logit_fit(w0[ok].astype(int), np.column_stack([dy, db])),
                            "coin_delta_both_0": C.coin(w0[ok][(dy == 0) & (db == 0)]),
                            "coin_by_dy": {str(k): C.coin(w0[ok][np.clip(dy, -3, 3) == k]) for k in range(-3, 4)}}
            # H3 strata
            sub = tab[msk]
            rr["H3"] = {"length": strata_table(sub, "lstr"), "charge": strata_table(sub, "charge"),
                        "missed": strata_table(sub, "missed"), "n_var": strata_table(sub, "n_var"),
                        "mass_decile": strata_table(sub, "mass", np.quantile(sub.mass, np.linspace(.1, .9, 9))),
                        "dec_is_reverse": strata_table(sub, "dec_is_reverse")}
            r[thr] = rr
        res[col] = r
        c5 = r["0.05"]["cumulative"]
        print(col, {k: (v["n"], round(v["coin"], 4) if v["coin"] else None, v.get("length_std", {}).get("coin_std"))
                    for k, v in c5.items()}, flush=True)
        print(col, "ladder", r["0.05"]["ladder"]["model_dy_db"], flush=True)
    json.dump(res, open(os.path.join(C.SCR, "step1.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
