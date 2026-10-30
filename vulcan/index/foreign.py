"""L1 index with arm-C (natural foreign proteome) entrapment (DESIGN §6.3 arm C) for M3 dev entrapment.

Identical to `vulcan.index.build.build` (copied from commit 166e6de) except that a list of neighbour-filtered foreign
peptides (I/L collapsed; from scripts/entrapment/build_dbs.py arm_foreign, the M1e `<name>.C.peptides.tsv` kind=E rows)
is added as entrapment family E (entrap_seed = -1).  Kept separate so that the shared builder (concurrently extended
by M2b) is not touched; to be folded into build.py as a `foreign_peptides` option.

    python -m vulcan.index.foreign --fasta a.fasta --foreign peptides.txt --out DIR
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time

import numpy as np
import pyarrow as pa

from . import decoys as D
from ._core import core
from .build import FAMILY, INDEX_VERSION, _write, content_hash, default_params


def build_foreign(fasta_paths, out_dir: str, foreign_peptides: str, params: dict | None = None,
                  workers: int | None = None, shuffle_rows_seed: int | None = None, swap_labels: bool = False) -> dict:
    """Build the index. shuffle_rows_seed / swap_labels exist ONLY for the G0 invariance test: they permute the
    pair list handed to expand_pairs and swap the forward/decoy roles of every pair; the index side must come
    out bit-identical (labels.arrow is the only output allowed to differ)."""
    p = default_params()
    p.update(params or {})
    workers = workers or D.default_workers()
    t0 = time.perf_counter()
    timing = {}
    text = "".join(open(f).read() for f in fasta_paths)
    e = p["enzyme"]
    dg = core.digest_fasta(text, e["missed_cleavages"], e["min_len"], e["max_len"], e["cleave_at"],
                           e["restrict"], e["c_terminal"])
    targets = list(dg["keys"])
    t_missed = dict(zip(targets, dg["missed"].tolist()))
    timing["digest_s"] = time.perf_counter() - t0

    # ---- families
    t1 = time.perf_counter()
    taken = set(targets)
    entrap = {}   # seed -> {target: (E, perm)}
    ent_stats = {}
    for s in p["entrap_seeds"]:
        m, nr, nd = D.make_entrapment(targets, s, taken, workers=workers)
        entrap[s] = m
        ent_stats[str(s)] = {"n": len(m), "redrawn": nr, "dropped_unpairable": nd}
    fwd_fam = {t: ("T", 0, t) for t in targets}           # seq -> (family, seed, target seq)
    for s, m in entrap.items():
        for t, (ent, _) in m.items():
            fwd_fam[ent] = ("E", s, t)
    # ---- arm C: natural foreign-proteome entrapment peptides (already neighbour-filtered against the targets by
    # scripts/entrapment/build_dbs.py arm_foreign).  Family E, entrap_seed -1, each foreign peptide its own pair key.
    with open(foreign_peptides) as fh:
        foreign = [x.strip() for x in fh if x.strip() and not x.startswith("#")]
    n_foreign = 0
    for f in foreign:
        if f in fwd_fam:
            continue
        fwd_fam[f] = ("E", -1, f)
        t_missed[f] = sum(1 for i in range(len(f) - 1) if f[i] in e["cleave_at"] and f[i + 1] not in e["restrict"])
        n_foreign += 1
    ent_stats["foreign"] = {"n": n_foreign, "file": os.path.abspath(foreign_peptides),
                            "sha256": hashlib.sha256(open(foreign_peptides, "rb").read()).hexdigest()}
    forward = sorted(fwd_fam)
    dec, n_dec_shuffled = D.make_decoys(forward, p["decoy_seed"])
    forward = [f for f in forward if f in dec]
    timing["families_s"] = time.perf_counter() - t1

    # ---- pairs -> peptidoforms (vulcan_core / sage-core)
    t2 = time.perf_counter()
    fw = forward
    dc = [dec[f][0] for f in fw]
    perms = [dec[f][1] for f in fw]
    order = np.arange(len(fw))
    if shuffle_rows_seed is not None:
        order = np.random.default_rng(shuffle_rows_seed).permutation(len(fw))
    a_seq = [fw[i] for i in order]
    b_seq = [dc[i] for i in order]
    pr = []
    for i in order:
        q = perms[i]
        if q is None:
            n = len(fw[i])
            q = np.array([n - 2 - j for j in range(n - 1)] + [n - 1], np.uint8)
        pr.append(q)
    if swap_labels:   # the decoy becomes the "forward" side: perm' = inverse(perm)
        a_seq, b_seq = b_seq, a_seq
        pr = [np.argsort(q).astype(np.uint8) for q in pr]
    poff = np.zeros(len(pr) + 1, np.int64)
    np.cumsum([len(q) for q in pr], out=poff[1:])
    pflat = np.concatenate(pr).astype(np.uint8) if pr else np.zeros(0, np.uint8)
    smods = sorted(p["static_mods"].items())
    vmods = sorted(p["variable_mods"].items())
    r = core.expand_pairs(a_seq, b_seq, pflat, poff, smods, vmods, p["max_variable_mods"],
                          p["peptide_min_mass"], p["peptide_max_mass"], p["id_seed"])
    timing["expand_s"] = time.perf_counter() - t2

    # ---- peptidoform arrays (both sides), content-defined order
    t3 = time.perf_counter()
    n_pairs = len(r["seq"])
    seq_i = r["seq"].astype(np.int64)
    # residue table: a-side sequences then b-side sequences
    a_bytes = np.frombuffer("".join(a_seq).encode(), np.uint8)
    b_bytes = np.frombuffer("".join(b_seq).encode(), np.uint8)
    lens = np.array([len(s) for s in a_seq], np.int64)
    a_off = np.zeros(len(a_seq) + 1, np.int64)
    np.cumsum(lens, out=a_off[1:])
    residues = np.concatenate([a_bytes, b_bytes])
    pf_id = np.concatenate([r["id_f"], r["id_d"]])
    mass32 = np.concatenate([r["mass32_f"], r["mass32_d"]])
    seq_off = np.concatenate([a_off[seq_i], a_off[-1] + a_off[seq_i]])
    length = np.concatenate([lens[seq_i], lens[seq_i]]).astype(np.uint8)
    mask = np.concatenate([r["mask_f"], r["mask_d"]])
    n_var = np.concatenate([r["n_var"], r["n_var"]])
    # missed cleavages: a pair attribute (forward target's digest count), identical for both members
    fwd_of = np.array([x if not swap_labels else y for x, y in zip(a_seq, b_seq)], dtype=object)
    tgt_of_fwd = [fwd_fam[f][2] for f in fwd_of]
    missed_seq = np.array([t_missed[t] for t in tgt_of_fwd], np.uint8)
    missed = np.concatenate([missed_seq[seq_i], missed_seq[seq_i]])
    side_b = np.concatenate([np.zeros(n_pairs, bool), np.ones(n_pairs, bool)])
    pair_row = np.concatenate([np.arange(n_pairs), np.arange(n_pairs)])
    assert len(np.unique(pf_id)) == len(pf_id), "peptidoform id collision"
    srt = np.lexsort((pf_id, mass32))
    inv = np.empty_like(srt)
    inv[srt] = np.arange(len(srt))
    partner_unsorted = np.concatenate([np.arange(n_pairs) + n_pairs, np.arange(n_pairs)])
    partner = inv[partner_unsorted][srt]

    # canonical residue table: rebuild residues in sorted order so the index side is order-free
    lens_s = length[srt].astype(np.int64)
    new_off = np.zeros(len(srt) + 1, np.int64)
    np.cumsum(lens_s, out=new_off[1:])
    gather = (np.repeat(seq_off[srt] - new_off[:-1], lens_s) + np.arange(new_off[-1]))
    residues_s = residues[gather]
    idx_tbl = pa.table({
        "pf_id": pa.array(pf_id[srt], pa.uint64()),
        "mass32": pa.array(mass32[srt], pa.float32()),
        "seq_off": pa.array(new_off[:-1], pa.int64()),
        "length": pa.array(length[srt], pa.uint8()),
        "mask": pa.array(mask[srt], pa.uint64()),
        "n_var": pa.array(n_var[srt], pa.uint8()),
        "missed": pa.array(missed[srt], pa.uint8()),
        "partner": pa.array(partner, pa.int64()),
    })
    timing["arrays_s"] = time.perf_counter() - t3

    # ---- labels (L4 only)
    fam_fwd = np.array([FAMILY[fwd_fam[f][0]] for f in fw], np.uint8)[order]
    seed_fwd = np.array([fwd_fam[f][1] for f in fw], np.int32)[order]
    tgt_fwd = np.array([fwd_fam[f][2] for f in fw], dtype=object)[order]
    fwd_seq = np.array(fw, dtype=object)[order]
    # forward side = the true forward peptide regardless of swap_labels
    is_dec = side_b if not swap_labels else ~side_b
    fam_side = np.concatenate([fam_fwd[seq_i], fam_fwd[seq_i]])
    fam = np.where(is_dec, fam_side + 2, fam_side)
    fwd_member_id = np.where(is_dec, pf_id[partner_unsorted], pf_id)
    tseq = np.concatenate([tgt_fwd[seq_i], tgt_fwd[seq_i]])
    entrap_pair = np.array([int.from_bytes(hashlib.blake2b(t.encode(), digest_size=8).digest(), "little")
                            for t in tgt_fwd], np.uint64)
    entrap_pair = np.concatenate([entrap_pair[seq_i], entrap_pair[seq_i]])
    lab = pa.table({
        "pf_id": pa.array(pf_id[srt], pa.uint64()),
        "family": pa.array(fam[srt].astype(np.uint8)),
        "entrap_seed": pa.array(np.concatenate([seed_fwd[seq_i], seed_fwd[seq_i]])[srt], pa.int32()),
        "is_decoy": pa.array(is_dec[srt]),
        "is_entrapment": pa.array((fam[srt] % 2) == 1),
        "tdc_pair": pa.array(fwd_member_id[srt], pa.uint64()),
        "entrap_pair": pa.array(entrap_pair[srt], pa.uint64()),
        "fwd_seq": pa.array(np.concatenate([fwd_seq[seq_i], fwd_seq[seq_i]])[srt].tolist(), pa.string()),
        "target_seq": pa.array(tseq[srt].tolist(), pa.string()),
    })
    del pair_row
    prot = pa.table({
        "peptide": pa.array(targets, pa.string()),
        "prot_off": pa.array(dg["prot_off"][:-1], pa.int64()),
        "prot_len": pa.array(np.diff(dg["prot_off"]), pa.int64()),
    })
    os.makedirs(out_dir, exist_ok=True)
    _write(idx_tbl, os.path.join(out_dir, "peptidoforms.arrow"))
    np.save(os.path.join(out_dir, "residues.npy"), residues_s)
    _write(lab, os.path.join(out_dir, "labels.arrow"))
    _write(prot, os.path.join(out_dir, "proteins.arrow"))
    np.save(os.path.join(out_dir, "protein_idx.npy"), dg["prot_idx"])
    with open(os.path.join(out_dir, "accessions.txt"), "w") as fh:
        fh.write("\n".join(dg["accessions"]) + "\n")
    timing["total_s"] = time.perf_counter() - t0
    fams = np.bincount(fam.astype(np.int64), minlength=4)
    meta = {
        "index_version": INDEX_VERSION, "sage_commit": core.SAGE_COMMIT, "params": p,
        "fasta": [os.path.abspath(f) for f in fasta_paths],
        "content_hash": content_hash(fasta_paths, {**p, "foreign_sha256": ent_stats["foreign"]["sha256"]}),
        "n_targets": len(targets), "n_nonstandard_dropped": int(dg["n_nonstandard"]),
        "n_forward_peptides": len(fw), "n_decoy_shuffled": n_dec_shuffled,
        "n_forward_without_decoy": len(fwd_fam) - len(fw),
        "entrapment": ent_stats,
        "n_peptidoforms": int(len(pf_id)),
        "n_peptidoforms_by_family": {k: int(fams[v]) for k, v in FAMILY.items()},
        "timing_s": timing,
        "g0_test_knobs": {"shuffle_rows_seed": shuffle_rows_seed, "swap_labels": swap_labels},
    }
    with open(os.path.join(out_dir, "meta.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    return meta


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", nargs="+", required=True)
    ap.add_argument("--foreign", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    print(json.dumps(build_foreign(a.fasta, a.out, a.foreign), indent=1)[:2000])


if __name__ == "__main__":
    main()
