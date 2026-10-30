#!/usr/bin/env python3
"""E9 arm-C diagnostic: are accepted foreign (Arabidopsis) entrapment peptides enriched for sequence similarity to human
targets (a 6-mer shared with any human target peptide, I/L collapsed; DESIGN §6.0 A1.3 homolog criterion)?  Compares
accepted E at the global 1% / 5% threshold vs a random sample of all foreign peptides in the index.  If enriched, the
E-over-dec(E) excess in arm C reflects residual homology (near-true matches), not a scorer tilt on true nulls.

    python scripts/e9/c_homology.py --n 1000 --out results/2026-10-27_E9_scale/c_homology.json
"""
import argparse, csv, json, os, re, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
import analyze as A  # noqa: E402
from vulcan.stats import global_fdr as G  # noqa: E402
from vulcan.stats._util import hash_tiebreak  # noqa: E402


def kmers(s, k=6):
    s = s.replace("I", "L")
    return {s[i:i + k] for i in range(len(s) - k + 1)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--work", default="/tmp/vulcan_e9")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    man = list(csv.DictReader(open(A.MAN), delimiter="\t"))
    want = [f"{r['id']}__{re.sub(r'[^A-Za-z0-9._-]', '_', r['filename'])[:-4]}" for r in man[:a.n]]
    lab = A.labels(A.IDX["C"])
    fam, seq = lab["family"], lab["target_seq"]
    human = set()
    for s in seq[fam == 0]:
        human |= kmers(s)
    rng = np.random.default_rng(11)
    allE = np.flatnonzero(fam == 1)
    samp = rng.choice(allE, 20000, replace=False)
    base = float(np.mean([bool(kmers(seq[i]) & human) for i in samp]))
    out = {"n_human_6mers": len(human), "random_foreign_frac_shared_6mer": base}
    for col in ("score", "hyperscore"):
        d, _, _ = A.load_arm(a.work, "C", want, col=col)
        rows = d["pf_row"]
        dec = lab["is_decoy"][rows]
        g = G.global_competition(d["run"], d["pf_id"], lab["tdc_pair"][rows], d["score_raw"], dec,
                                 hash_tiebreak(d["pf_id"], 11))
        for al in (0.01, 0.05):
            acc = g.accepted(al)
            r = rows[g.row[acc]]
            e = r[fam[r] == 1]
            # decoys of foreign peptides above the same threshold: their base (foreign) sequence
            s_cut = g.score[acc].min()
            dr = rows[g.row[(g.score >= s_cut) & g.is_decoy]]
            de = dr[fam[dr] == 3]
            fe = float(np.mean([bool(kmers(seq[i]) & human) for i in e])) if len(e) else None
            fd = float(np.mean([bool(kmers(seq[i]) & human) for i in de])) if len(de) else None
            out[f"{col}_{al}"] = {"n_E": int(len(e)), "E_frac_shared_6mer": fe, "n_decE": int(len(de)),
                                  "decE_base_frac_shared_6mer": fd}
            print(col, al, out[f"{col}_{al}"], flush=True)
    json.dump(out, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
