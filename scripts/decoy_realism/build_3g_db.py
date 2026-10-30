#!/usr/bin/env python3
"""3g step 3 (PLAN.md): foreign families for the synthetic-control index `human_contam_3g`.

  A  the existing arm-C Arabidopsis list (/tmp/vulcan_m3/foreign/human_contam.C.foreign.txt), unchanged
  B  distant natural prokaryotic union (E. coli K-12, P. furiosus, B. cereus, S. aureus USA300, Synechococcus PCC 7002;
     FASTA sequence space only), arm-C neighbour filters (scripts/entrapment/build_dbs.arm_foreign rules) against
     human targets AND against A (I/L-identical, isobaric N<->GG/Q<->AG, single substitution, = pseudo-reverse,
     = shuffled decoys 1001-1005, shuffled-decoy collisions)
  S  shuffled A: production paired-entrapment shuffle (vulcan.index.decoys.make_entrapment, seed 2001, C-term fixed)
     of every A peptide, unique against T, A, B (and their pseudo-reverses)

Writes WORK/foreign_3g.txt (A|B|S union, one per line), WORK/family_3g.tsv (seq, family), WORK/db_3g.json, then builds
the L1 index with the production builder vulcan.index.foreign.build_foreign (decoy recipe unchanged, seed 1001).
    conda_envs/pe/bin/python scripts/decoy_realism/build_3g_db.py [--work /tmp/vulcan_3g]
"""
import argparse, hashlib, json, os, sys, time
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "entrapment"))
from build_dbs import target_peptides  # noqa: E402
from peptides import isobaric_key, pseudo_reverse, shuffle_decoy_set, shuffle_hits  # noqa: E402

HUMAN = os.path.join(ROOT, "data", "cache", "M1e", "entrapment", "inputs", "human_contam.fasta")
A_LIST = "/tmp/vulcan_m3/foreign/human_contam.C.foreign.txt"
B_FASTAS = ["ecoli_K12_UP000000625.fasta", "pfu_UP000001013.fasta", "dn_UP000001417.fasta", "dn_UP000001939.fasta",
            "dn_UP000001688.fasta"]


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def neighbour_filter(F, ref, removed):
    F = F - ref
    removed["IL_identical"] = removed.get("IL_identical", 0)
    riso = {isobaric_key(t) for t in ref}
    F2 = {f for f in F if isobaric_key(f) not in riso}
    removed["isobaric"] = len(F) - len(F2)
    F = F2
    keys = set()
    for t in ref:
        for i in range(len(t)):
            keys.add(hash((t[:i], t[i + 1:], len(t))))
    F2 = {f for f in F if not any(hash((f[:i], f[i + 1:], len(f))) in keys for i in range(len(f)))}
    removed["single_substitution"] = len(F) - len(F2)
    del keys
    F = F2
    rev = {pseudo_reverse(t) for t in ref}
    F2 = {f for f in F if f not in rev}
    removed["equal_to_ref_pseudo_reverse"] = len(F) - len(F2)
    F = F2
    ST = shuffle_decoy_set(ref)
    F2 = {f for f in F if f not in ST}
    removed["equal_to_ref_shuffled_decoy"] = len(F) - len(F2)
    return F2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="/tmp/vulcan_3g")
    ap.add_argument("--no-index", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    t0 = time.time()
    T = set(target_peptides(HUMAN)[0])
    A = {x.strip() for x in open(A_LIST) if x.strip() and not x.startswith("#")}
    Braw = set()
    for f in B_FASTAS:
        Braw |= set(target_peptides(os.path.join(ROOT, "data", "fasta", f))[0])
    meta = {"human_fasta": HUMAN, "human_sha256": sha(HUMAN), "A_list": A_LIST, "A_sha256": sha(A_LIST),
            "B_fastas": {f: sha(os.path.join(ROOT, "data", "fasta", f)) for f in B_FASTAS},
            "n_T": len(T), "n_A": len(A), "n_B_raw": len(Braw)}
    rem = {}
    n0 = len(Braw)
    B = Braw - T - A
    rem["IL_identical_to_T_or_A"] = n0 - len(B)
    B = neighbour_filter(B, T | A, rem)
    bad = shuffle_hits(B, T | A | B)
    B = B - bad
    rem["shuffled_decoy_collides"] = len(bad)
    meta["B_removed"] = rem
    meta["n_B"] = len(B)
    print("B", len(B), rem, flush=True)
    from vulcan.index.decoys import make_entrapment
    taken = set(T) | A | B
    S_map, nr, nd = make_entrapment(sorted(A), 2001, taken, workers=16)
    S = {e for e, _ in S_map.values()}
    meta.update({"n_S": len(S), "S_redrawn": nr, "S_dropped": nd, "S_seed": 2001})
    assert not (S & (T | A | B)) and not (B & (T | A))
    fam = {**{s: "A" for s in A}, **{s: "B" for s in B}, **{s: "S" for s in S}}
    fl = os.path.join(a.work, "foreign_3g.txt")
    with open(fl, "w") as fh:
        fh.write("\n".join(sorted(fam)) + "\n")
    with open(os.path.join(a.work, "family_3g.tsv"), "w") as fh:
        fh.write("seq\tfamily\tsource\n")
        src = {e: t for t, (e, _) in S_map.items()}
        for s in sorted(fam):
            fh.write(f"{s}\t{fam[s]}\t{src.get(s, '')}\n")
    meta["foreign_list_sha256"] = sha(fl)
    meta["prep_s"] = time.time() - t0
    json.dump(meta, open(os.path.join(a.work, "db_3g.json"), "w"), indent=1)
    print(json.dumps(meta, indent=1), flush=True)
    if not a.no_index:
        from vulcan.index.foreign import build_foreign
        out = os.path.join(a.work, "index", "human_contam_3g")
        info = build_foreign([HUMAN], out, fl)
        meta["index"] = {k: info.get(k) for k in ("content_hash", "n_peptidoforms", "n_peptidoforms_by_family",
                                                    "n_decoy_shuffled", "timing_s")} if isinstance(info, dict) else str(info)
        json.dump(meta, open(os.path.join(a.work, "db_3g.json"), "w"), indent=1, default=str)
        print("index", meta["index"], flush=True)


if __name__ == "__main__":
    main()
