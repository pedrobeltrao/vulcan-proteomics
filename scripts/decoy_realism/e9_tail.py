#!/usr/bin/env python3
"""3g steps 1-2 input: E9 arm C / arm P global competitions (N = 1000 prefix, frozen score and hyperscore; exactly the
E9 L6 procedure) -> coins at the global 1%/5% thresholds (re-derivation check vs E9 metrics) and a table of every
entrapment-family pair winner (E or dec(E)) at/above the 5% cut, with its best-evidence PSM (run, spec_idx, pf_row),
base sequence, actual decoy sequence, length, charge, missed cleavages, mass, n_var, collision-shuffle flag, project.
    conda_envs/pe/bin/python scripts/decoy_realism/e9_tail.py [--n 1000]
"""
import argparse, json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import common as C  # noqa: E402
from vulcan.index.decoys import pseudo_reverse  # noqa: E402

E_FAM, DE_FAM = 1, 3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--work", default="/tmp/vulcan_e9")
    ap.add_argument("--arms", nargs="+", default=["C", "P"])
    a = ap.parse_args()
    os.makedirs(C.SCR, exist_ok=True)
    man = C.manifest()[:a.n]
    runs = [r["run"] for r in man]
    proj = np.array([r["accession"] for r in man])
    res = {}
    for arm in a.arms:
        lab = C.labels(C.IDX[arm])
        fam = lab["family"]
        for col in ("score", "hyperscore"):
            t0 = time.time()
            d = C.load_winners(a.work, arm, runs, col)
            g, cut = C.compete(d, lab)
            rows = d["pf_row"][g.row]
            f = fam[rows]
            out = {"cut": cut, "n_spectrum_winners": int(len(d["pf_row"]))}
            for al, c in cut.items():
                m = (g.score >= c) & np.isin(f, (E_FAM, DE_FAM))
                out[str(al)] = C.coin(f[m] == E_FAM)
            # tail table at the 5% cut
            m = (g.score >= cut[0.05]) & np.isin(f, (E_FAM, DE_FAM))
            gi = g.row[m]
            r = rows[m]
            base = lab["target_seq"][r]
            wseq = np.array([C.row_seq(lab, x) for x in r])
            is_dec = lab["is_decoy"][r]
            dseq = np.where(is_dec, wseq, None)
            # decoy sequence for E winners: the pair partner's sequence (row of the decoy with the same tdc_pair)
            need = np.flatnonzero(~is_dec)
            if len(need):
                pair_all = lab["tdc_pair"]
                dec_rows = np.flatnonzero(lab["is_decoy"] & np.isin(pair_all, lab["tdc_pair"][r[need]]))
                pmap = {}
                for x in dec_rows:
                    pmap.setdefault(int(pair_all[x]), x)
                for j in need:
                    x = pmap.get(int(lab["tdc_pair"][r[j]]))
                    dseq[j] = C.row_seq(lab, x) if x is not None else None
            tab = pd.DataFrame({
                "run": d["run"][gi], "spec_idx": d["spec_idx"][gi], "pf_row": r, "pf_id": d["pf_id"][gi],
                "score": g.score[m], "E_won": ~is_dec, "above1": g.score[m] >= cut[0.01],
                "base_seq": base, "win_seq": wseq, "dec_seq": dseq, "length": lab["length"][r],
                "charge": d["charge"][gi], "missed": lab["missed"][r], "n_var": lab["n_var"][r],
                "mass": lab["mass"][r], "project": proj[d["run"][gi]]})
            # strip mod-free sequence comparisons: base_seq is the unmodified forward sequence
            tab["dec_is_reverse"] = [ds is not None and ds == pseudo_reverse(b) for b, ds in zip(tab.base_seq, tab.dec_seq)]
            tab.to_parquet(os.path.join(C.SCR, f"e9tail_{arm}_{col}.parquet"))
            # project-cluster CI for the headline coins
            for al in (0.01, 0.05):
                mm = tab.score >= cut[al]
                out[str(al)]["cluster_ci95_project"] = C.cluster_boot(tab.E_won[mm].values, tab.project[mm].values)
            out["wall_s"] = time.time() - t0
            res[f"{arm}_{col}"] = out
            print(arm, col, json.dumps({k: out[k] for k in ("0.01", "0.05")}), flush=True)
            del d, g
    json.dump(res, open(os.path.join(C.SCR, "e9_coins.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
