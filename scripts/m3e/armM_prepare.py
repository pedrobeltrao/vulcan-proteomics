#!/usr/bin/env python3
"""Arm M on dev (DESIGN §6.3 (M), A1.5): D7 yeast (PXD031616) searched against yeast + contaminants + the spiked human
alpha-synuclein P37840 (true sample components) as targets, with human Swiss-Prot peptides as entrapment.
Human peptides are neighbour-filtered against the targets exactly as arm C (scripts/entrapment/build_dbs.py
arm_foreign: I/L-identical, isobaric N<->GG / Q<->AG, single-substitution, decoy collisions), then restricted to the
effective search space (in_space) and seeded-subsampled to r = 1 (|E_eff| = |T_eff|), one replicate per
ENTRAP_SEEDS_DEV seed (2001, 2002, 2003; numpy default_rng(seed)).  One L1 index per seed (vulcan.index.foreign).
    python scripts/m3e/armM_prepare.py"""
import json, os, sys
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "entrapment"))
from peptides import in_space, mass  # noqa: E402
W = os.environ.get("VULCAN_M3E_WORK", "/tmp/vulcan_m3e")
DB = os.path.join(W, "armM", "db")
SEEDS = (2001, 2002, 2003)


def main():
    T, E = [], []
    with open(os.path.join(DB, "yeast_contam.C.peptides.tsv")) as fh:
        next(fh)
        for line in fh:
            seq, kind = line.split("\t", 2)[:2]
            (T if kind == "T" else E if kind == "E" else []).append(seq)
    Teff = [t for t in T if in_space(t)]
    Eeff = sorted(e for e in E if in_space(e))
    info = {"n_target": len(T), "n_target_eff": len(Teff), "n_foreign_filtered": len(E), "n_foreign_eff": len(Eeff),
            "seeds": {}}
    from vulcan.index.foreign import build_foreign
    tm = np.array([mass(t) for t in Teff])
    edges = np.quantile(tm, np.linspace(0, 1, 11)); edges[0], edges[-1] = -np.inf, np.inf
    ct = np.histogram(tm, edges)[0]
    for s in SEEDS:
        pick = np.sort(np.random.default_rng(s).choice(len(Eeff), size=len(Teff), replace=False))
        sub = [Eeff[i] for i in pick]
        fl = os.path.join(W, "armM", f"foreign_{s}.txt")
        with open(fl, "w") as fh:
            fh.write("\n".join(sub) + "\n")
        em = np.array([mass(e) for e in sub])
        ce = np.histogram(em, edges)[0]
        d = os.path.join(W, "index", f"yeast_M{s}")
        if not os.path.exists(os.path.join(d, "meta.json")):
            build_foreign([os.path.join(W, "armM", "yeast_contam_snca.fasta")], d, fl)
        meta = json.load(open(os.path.join(d, "meta.json")))
        info["seeds"][s] = {"r_peptide_eff": len(sub) / len(Teff), "r_per_mass_decile": (ce / ct).tolist(),
                            "index": d, "n_peptidoforms_by_family": meta.get("n_peptidoforms_by_family")}
        print(s, json.dumps(info["seeds"][s]), flush=True)
    info["decile_mass_edges_da"] = edges[1:-1].tolist()
    json.dump(info, open(os.path.join(W, "armM", "armM_db.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
