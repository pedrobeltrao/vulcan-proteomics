#!/usr/bin/env python3
"""3g step 3 (PLAN.md): synthetic controls on the human_contam_3g search of D8 runs 0-99 (frozen static_gbm_pred and
hyperscore; L6 global competition, N = all scored runs).  Families of foreign pair base sequences: A (Arabidopsis),
B (distant prokaryotes), S (shuffled A).  Coins c_A, c_B, c_S at the global 5% (primary) / 1% cuts; c_A - c_B;
odds ratios R_AB, R_AS (= decoy-normalised natural-vs-natural / natural-vs-shuffled share ratios); L1 (Hamming<=2)
and L2 (6-mer) homology filters per family; H3 strata (length) per family.
    conda_envs/pe/bin/python scripts/decoy_realism/step3.py
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import common as C  # noqa: E402
import homology as HM  # noqa: E402

WORK = "/tmp/vulcan_3g"
IDX3 = os.path.join(WORK, "index", "human_contam_3g")


def odds(w):
    p = np.mean(w)
    return p / (1 - p)


def boot_or(w1, w2, n_boot=2000):
    rng = np.random.default_rng(C.SEED)
    rs = [odds(w1[rng.integers(0, len(w1), len(w1))]) / odds(w2[rng.integers(0, len(w2), len(w2))])
          for _ in range(n_boot)]
    return {"OR": float(odds(w1) / odds(w2)), "ci95": [float(np.quantile(rs, .025)), float(np.quantile(rs, .975))]}


def main():
    man = C.manifest()[:100]
    runs = [r["run"] for r in man]
    present = [r for r in runs if os.path.exists(os.path.join(WORK, "out", "P", r + ".scored.parquet"))]
    lab = C.labels(IDX3)
    famtab = pd.read_csv(os.path.join(WORK, "family_3g.tsv"), sep="\t")
    fmap = dict(zip(famtab.seq, famtab.family))
    prots = HM.read_fasta(HM.HUMAN)
    k6 = {p[i:i + 6] for p in prots for i in range(len(p) - 5)}
    hum = np.unique(lab["target_seq"][lab["family"] == 0])
    byL = {}
    for s in hum:
        byL.setdefault(len(s), []).append(HM.il(s))
    res = {"n_runs_present": len(present), "n_runs_requested": len(runs)}
    for col in ("score", "hyperscore"):
        d = C.load_winners(WORK, "P", runs, col)
        g, cut = C.compete(d, lab)
        rows = d["pf_row"][g.row]
        f = lab["family"][rows]
        m5 = (g.score >= cut[0.05]) & np.isin(f, (1, 3))
        tab = pd.DataFrame({"E_won": f[m5] == 1, "above1": g.score[m5] >= cut[0.01],
                            "base": lab["target_seq"][rows[m5]], "length": lab["length"][rows[m5]]})
        tab["fam"] = tab.base.map(fmap)
        assert tab.fam.notna().all()
        # homology flags for base seqs (L1 Hamming<=2, L2 6-mer)
        u = sorted(set(tab.base.map(HM.il)))
        maxL = 51
        ref_off, ref_n, blocks, pos = np.zeros(maxL, np.int64), np.zeros(maxL, np.int64), [], 0
        for L in range(maxL):
            ss = byL.get(L, [])
            ref_off[L], ref_n[L] = pos, len(ss)
            if ss:
                blocks.append(HM.CODE[np.frombuffer("".join(ss).encode(), np.uint8)])
                pos += L * len(ss)
        qs, qo, ql = HM.encode(u)
        ham = dict(zip(u, HM.ham_many(qs, qo, ql, np.concatenate(blocks), ref_off, ref_n, 2)))
        tab["L1"] = tab.base.map(lambda s: ham[HM.il(s)] <= 2)
        tab["k6"] = tab.base.map(lambda s: any(HM.il(s)[i:i + 6] in k6 for i in range(len(s) - 5)))
        tab["L2"] = tab.L1 | tab.k6
        r = {"cut": cut}
        for thr in ("0.05", "0.01"):
            t = tab if thr == "0.05" else tab[tab.above1]
            rr = {}
            for fam in ("A", "B", "S"):
                tf = t[t.fam == fam]
                rr[fam] = {"all": C.coin(tf.E_won.values),
                           "L1_kept": C.coin(tf.E_won.values[~tf.L1.values]),
                           "L2_kept": C.coin(tf.E_won.values[~tf.L2.values]),
                           "frac_removed_L2": float(tf.L2.mean()) if len(tf) else None,
                           "by_length": {str(k): C.coin(tf.E_won.values[C.len_stratum(tf.length.values) == k])
                                         for k in range(5)}}
            wA = t[t.fam == "A"].E_won.values.astype(float)
            wB = t[t.fam == "B"].E_won.values.astype(float)
            wS = t[t.fam == "S"].E_won.values.astype(float)
            rr["cA_minus_cB"] = C.boot_diff(wA, wB)
            rr["cA_minus_cS"] = C.boot_diff(wA, wS)
            rr["R_AB_oddsratio"] = boot_or(wA, wB)
            rr["R_AS_oddsratio"] = boot_or(wA, wS)
            rr["cB_excess_over_cA_excess"] = C.boot_ratio(wA, wB)
            # decoy-normalised shares (natural vs natural): A share among A+B winners vs decA share among decA+decB
            nA, nB = int(wA.sum()), int(wB.sum())
            dA, dB = int(len(wA) - wA.sum()), int(len(wB) - wB.sum())
            rr["shares"] = {"A_of_AB_fwd": nA / max(nA + nB, 1), "decA_of_AB_dec": dA / max(dA + dB, 1),
                            "A_of_AS_fwd": nA / max(nA + int(wS.sum()), 1),
                            "decA_of_AS_dec": dA / max(dA + int(len(wS) - wS.sum()), 1)}
            r[thr] = rr
            print(col, thr, {k: (rr[k]["all"]["n"], round(rr[k]["all"]["coin"] or 0, 4)) for k in "ABS"},
                  rr["cA_minus_cB"], flush=True)
        res[col] = r
    json.dump(res, open(os.path.join(C.SCR, "step3.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
