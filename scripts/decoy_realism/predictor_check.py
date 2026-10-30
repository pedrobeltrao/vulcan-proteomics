#!/usr/bin/env python3
"""3g step 2e (PLAN.md): spectrum-free predictor check.  For random samples (seed 11) of 5,000 pairs each of
  C: Arabidopsis E vs its decoy        P: shuffled-human E_2001 vs its decoy        T: human target vs its decoy
(unmodified peptidoforms) predict with AlphaPeptDeep (generic models, NCE 29, charge 2, QE) and compare paired:
  ms2 entropy (normalised over b/y z1 ions), top-6 intensity share, y-ion share, predicted RT, and the expected RT
  residual against the human-target predicted-RT distribution (E|X - rt_pred|, X ~ 20k human target rt_pred: the
  spectrum-free analogue of rt_abs_resid for a null spectrum eluting like a human peptide), and the M3e G1
  sequence-classifier logit (/tmp/vulcan_m3e/seqclf.pkl).
    conda_envs/pred/bin/python scripts/decoy_realism/predictor_check.py
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "m3e"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "m3"))
import numpy as np  # noqa: E402
import common as C  # noqa: E402


def sample_pairs(lab, fam_fwd, n, rng):
    fam = lab["family"]
    cand = np.flatnonzero((fam == fam_fwd) & (lab["n_var"] == 0))
    fw = rng.choice(cand, n, replace=False)
    dec_fam = fam_fwd + 2
    dmask = (fam == dec_fam) & (lab["n_var"] == 0)
    drow = np.flatnonzero(dmask)
    pmap = dict(zip(lab["tdc_pair"][drow].tolist(), drow.tolist()))
    dw = np.array([pmap.get(int(lab["tdc_pair"][r]), -1) for r in fw])
    ok = dw >= 0
    return fw[ok], dw[ok]


def stats(apd, lab, rows):
    L = lab["length"][rows].astype(np.int64)
    res = np.asarray(lab["residues"])
    o = apd.predict(res, lab["seq_off"][rows], L, np.zeros(len(rows), np.int64), np.full(len(rows), 2), 29.0)
    ent, top6, yshare = [], [], []
    for i in range(len(rows)):
        p = o["ms2"][o["ms2_off"][i]:o["ms2_off"][i + 1]].astype(float).reshape(-1, 4)[:, [0, 2]]   # b_z1, y_z1
        s = p.sum()
        if s <= 0:
            ent.append(np.nan); top6.append(np.nan); yshare.append(np.nan)
            continue
        q = (p / s).ravel()
        q = q[q > 0]
        ent.append(float(-(q * np.log(q)).sum()))
        top6.append(float(np.sort(p.ravel())[::-1][:6].sum() / s))
        yshare.append(float(p[:, 1].sum() / s))
    return {"entropy": np.array(ent), "top6": np.array(top6), "yshare": np.array(yshare), "rt": o["rt"].astype(float)}


def paired(a, b):
    d = a - b
    d = d[np.isfinite(d)]
    nz = d[d != 0]
    return {"mean_fwd": float(np.nanmean(a)), "mean_dec": float(np.nanmean(b)), "mean_diff": float(d.mean()),
            "frac_fwd_greater": float((nz > 0).mean()), "ci95": C.wilson(int((nz > 0).sum()), len(nz)), "n": int(len(d))}


def main():
    from vulcan.predict.apd import APD
    apd = APD()
    rng = np.random.default_rng(C.SEED)
    out = {}
    labC = C.labels(C.IDX["C"])
    labP = C.labels(C.IDX["P"])
    tr = np.flatnonzero((labC["family"] == 0) & (labC["n_var"] == 0))
    hum = stats(apd, labC, rng.choice(tr, 20000, replace=False))["rt"]
    hs = np.sort(hum)
    cs = np.concatenate([[0], np.cumsum(hs)])

    def exp_resid(x):            # E|X - x| for X ~ empirical human target rt_pred
        k = np.searchsorted(hs, x)
        n = len(hs)
        return (x * k - cs[k] + (cs[n] - cs[k]) - x * (n - k)) / n
    try:
        import seqclf
        clf = seqclf.load()
    except Exception as e:  # noqa: BLE001
        clf = None
        out["seqclf_error"] = repr(e)
    for name, lab, ff in (("C_arabidopsis", labC, 1), ("P_shuffled", labP, 1), ("T_human", labC, 0)):
        fw, dw = sample_pairs(lab, ff, 5000, rng)
        a, b = stats(apd, lab, fw), stats(apd, lab, dw)
        r = {k: paired(a[k], b[k]) for k in a}
        ra, rb = exp_resid(a["rt"]), exp_resid(b["rt"])
        r["exp_rt_resid_vs_human"] = paired(-ra, -rb)       # 'fwd greater' = fwd has the SMALLER expected residual
        r["exp_rt_resid_vs_human"]["note"] = "sign flipped: frac_fwd_greater = P(fwd expected |RT resid| < dec)"
        if clf is not None:
            sa = [C.row_seq(lab, x) for x in fw]
            sb = [C.row_seq(lab, x) for x in dw]
            r["seqclf_logit"] = paired(clf.logit(sa, a["rt"]), clf.logit(sb, b["rt"]))
        out[name] = r
        print(name, json.dumps({k: (round(v["mean_diff"], 4), round(v["frac_fwd_greater"], 3)) for k, v in r.items()}),
              flush=True)
    json.dump(out, open(os.path.join(C.SCR, "predictor_check.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
