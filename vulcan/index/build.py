"""L1 proteome index builder (DESIGN.md §4 L1).

    python -m vulcan.index.build --fasta a.fasta [b.fasta ...] --out DIR [--entrap-seeds 2001 ...]

Pipeline (no per-spectrum work; per-peptide Python only for decoy/entrapment recipes):
  1. sage-core digestion (vulcan_core.digest_fasta): unique I/L-collapsed peptides + protein CSR.
  2. families: T (targets), E_s (paired shuffled entrapment per seed s), dec(T), dec(E_s) (vulcan.index.decoys).
  3. vulcan_core.expand_pairs: Sage's variable-mod enumeration on the forward side, carried through the decoy
     permutation -> one decoy peptidoform per forward peptidoform (same composition, same mass);
     hash-assigned u64 ids (xxhash64 of the canonical peptidoform string), Sage-f32 masses.
  4. sort by (Sage f32 mass, id): storage order is content-defined, never label- or input-order-defined.

Output (directory, Arrow IPC files):
  peptidoforms.arrow  INDEX side, label-free (what L1-L3 see): pf_id, mass32 (Sage f32), seq_off, length, mask,
                      n_var, missed, partner (row of the pair partner; symmetric: partner[partner[i]] = i)
  residues.npy        u8 ASCII residues of all sequences (rows point into it via seq_off/length)
  labels.arrow        L4 ONLY: pf_id, family (0 T, 1 E, 2 dec(T), 3 dec(E)), entrap_seed, is_decoy,
                      is_entrapment, tdc_pair (pf_id of the forward member), entrap_pair (u64 key shared by a
                      target and its entrapment, same for their decoys), fwd_seq (str), target_seq (str)
  proteins.arrow      L4 ONLY: target peptide (I/L collapsed) -> protein accessions
  meta.json           params, counts, content hash, timings
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

INDEX_VERSION = "vulcan-l1-0.1"
FAMILY = {"T": 0, "E": 1, "DT": 2, "DE": 3}


def _write(table: pa.Table, path: str) -> None:
    with pa.OSFile(path + ".partial", "wb") as sink, pa.ipc.new_file(sink, table.schema) as w:
        w.write_table(table)
    os.replace(path + ".partial", path)


def _read(path: str) -> pa.Table:
    with pa.OSFile(path, "rb") as src:
        return pa.ipc.open_file(src).read_all()


def default_params() -> dict:
    """M1e baseline search space (scripts/baselines/run_searches.py): trypsin KR/P, 2 missed, len 7-50,
    500-5000 Da, CAM static, Met-ox variable (<= 2)."""
    return {
        "enzyme": {"cleave_at": "KR", "restrict": "P", "c_terminal": True, "missed_cleavages": 2,
                   "min_len": 7, "max_len": 50},
        "static_mods": {"C": 57.021464},
        "variable_mods": {"M": 15.9949},
        "max_variable_mods": 2,
        "peptide_min_mass": 500.0, "peptide_max_mass": 5000.0,
        "decoy_seed": D.DECOY_SEED, "entrap_seeds": [], "id_seed": 0,
    }


def content_hash(fasta_paths, params) -> str:
    h = hashlib.sha256()
    for p in fasta_paths:
        with open(p, "rb") as fh:
            for b in iter(lambda: fh.read(1 << 22), b""):
                h.update(b)
    h.update(json.dumps(params, sort_keys=True).encode())
    h.update(INDEX_VERSION.encode())
    h.update(core.SAGE_COMMIT.encode())
    return h.hexdigest()[:16]


def build(fasta_paths, out_dir: str, params: dict | None = None, workers: int | None = None,
          shuffle_rows_seed: int | None = None, swap_labels: bool = False) -> dict:
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
        "content_hash": content_hash(fasta_paths, p),
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


def _missed(pep: str, e: dict) -> int:
    """Missed cleavages of a peptide not from the digest (arm-C foreign peptides): internal cleavage sites."""
    cut, rest = e["cleave_at"], e.get("restrict") or ""
    return sum(1 for i in range(len(pep) - 1) if pep[i] in cut and pep[i + 1] not in rest)


def build_multi_decoy(fasta_paths, out_dir: str, params: dict | None = None, workers: int | None = None,
                      foreign: list[str] | None = None) -> dict:
    """Multi-replicate index for family top-2 (DESIGN §4 L1, §6.3): DB = T ∪ E_s ∪ dec_j(T) ∪ dec_j(E_s) for every
    entrapment seed s in params['entrap_seeds'] and every decoy seed j in params['decoy_seeds'] (C-term-fixed
    hash(peptide, seed_j, k) shuffles, vulcan.index.decoys.make_shuffled_decoys; no pseudo-reverse set).

    Same output layout as build(), plus:
      peptidoforms.arrow  partner   = pair partner in decoy set 0 (self for decoy rows of other sets)
                          partner_j = pair partner in decoy set j >= 1 (same convention)
      labels.arrow        decoy_set (-1 forward, j for dec_j), decoy_seed (0 forward)
    One search over this index + vulcan.gpu.family reconstructs every (decoy seed x entrapment seed) replicate.

    foreign: optional arm-C natural foreign peptides (I/L collapsed, neighbour-filtered upstream, e.g. by
    scripts/entrapment/build_dbs.py --arms C). They become forward family E with entrap_seed -1 (target_seq = itself,
    no pairing) and get shuffled decoys of every decoy seed like any forward peptide. Foreign peptides equal to a
    target are dropped; paired entrapments are drawn after them, so no entrapment equals a foreign peptide."""
    p = default_params()
    p.update(params or {})
    dseeds = list(p["decoy_seeds"])
    J = len(dseeds)
    assert J >= 1
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
    t1 = time.perf_counter()
    taken = set(targets)
    fwd_fam = {t: ("T", 0, t) for t in targets}
    ent_stats = {}
    if foreign:
        fset = sorted(set(foreign) - taken)
        ent_stats["foreign"] = {"n_input": len(set(foreign)), "n": len(fset)}
        for f in fset:
            fwd_fam[f] = ("E", -1, f)
        taken.update(fset)
    for s_ in p["entrap_seeds"]:
        m, nr, nd = D.make_entrapment(targets, s_, taken, workers=workers)
        ent_stats[str(s_)] = {"n": len(m), "redrawn": nr, "dropped_unpairable": nd}
        for t, (ent, _) in m.items():
            fwd_fam[ent] = ("E", s_, t)
    forward = sorted(fwd_fam)
    used: set[str] = set()
    decs, dec_stats = [], {}
    for ds in dseeds:
        dj, nr, nd = D.make_shuffled_decoys(forward, ds, used, workers=workers)
        decs.append(dj)
        dec_stats[str(ds)] = {"n": len(dj), "redrawn": nr, "dropped": nd}
    fw = [f for f in forward if all(f in dj for dj in decs)]
    timing["families_s"] = time.perf_counter() - t1

    t2 = time.perf_counter()
    smods = sorted(p["static_mods"].items())
    vmods = sorted(p["variable_mods"].items())
    res = []
    for dj in decs:
        b_seq = [dj[f][0] for f in fw]
        pr = [dj[f][1] for f in fw]
        poff = np.zeros(len(pr) + 1, np.int64)
        np.cumsum([len(q) for q in pr], out=poff[1:])
        r = core.expand_pairs(fw, b_seq, np.concatenate(pr).astype(np.uint8), poff, smods, vmods,
                              p["max_variable_mods"], p["peptide_min_mass"], p["peptide_max_mass"], p["id_seed"])
        res.append((b_seq, r))
    r0 = res[0][1]
    for _, r in res[1:]:
        assert np.array_equal(r["seq"], r0["seq"]) and np.array_equal(r["id_f"], r0["id_f"]), "forward side differs"
    timing["expand_s"] = time.perf_counter() - t2

    t3 = time.perf_counter()
    nf = len(r0["seq"])
    seq_i = r0["seq"].astype(np.int64)
    lens = np.array([len(x) for x in fw], np.int64)
    a_off = np.zeros(len(fw) + 1, np.int64)
    np.cumsum(lens, out=a_off[1:])
    L = int(a_off[-1])
    residues = np.concatenate([np.frombuffer("".join(fw).encode(), np.uint8)] +
                              [np.frombuffer("".join(b).encode(), np.uint8) for b, _ in res])
    pf_id = np.concatenate([r0["id_f"]] + [r["id_d"] for _, r in res])
    mass32 = np.concatenate([r0["mass32_f"]] + [r["mass32_d"] for _, r in res])
    seq_off = np.concatenate([a_off[seq_i] + j * L for j in range(J + 1)])
    length = np.tile(lens[seq_i], J + 1).astype(np.uint8)
    mask = np.concatenate([r0["mask_f"]] + [r["mask_d"] for _, r in res])
    n_var = np.tile(r0["n_var"], J + 1)
    tgt_of_fwd = [fwd_fam[f][2] for f in fw]
    missed_seq = np.array([t_missed[t] if t in t_missed else _missed(t, e) for t in tgt_of_fwd], np.uint8)
    missed = np.tile(missed_seq[seq_i], J + 1)
    assert len(np.unique(pf_id)) == len(pf_id), "peptidoform id collision"
    N = len(pf_id)
    srt = np.lexsort((pf_id, mass32))
    inv = np.empty_like(srt)
    inv[srt] = np.arange(N)
    base = np.arange(nf)
    partners = []
    for j in range(J):
        pu = np.arange(N)
        pu[:nf] = (j + 1) * nf + base
        pu[(j + 1) * nf:(j + 2) * nf] = base
        partners.append(inv[pu][srt])
    lens_s = length[srt].astype(np.int64)
    new_off = np.zeros(N + 1, np.int64)
    np.cumsum(lens_s, out=new_off[1:])
    gather = np.repeat(seq_off[srt] - new_off[:-1], lens_s) + np.arange(new_off[-1])
    residues_s = residues[gather]
    cols = {
        "pf_id": pa.array(pf_id[srt], pa.uint64()), "mass32": pa.array(mass32[srt], pa.float32()),
        "seq_off": pa.array(new_off[:-1], pa.int64()), "length": pa.array(length[srt], pa.uint8()),
        "mask": pa.array(mask[srt], pa.uint64()), "n_var": pa.array(n_var[srt], pa.uint8()),
        "missed": pa.array(missed[srt], pa.uint8()), "partner": pa.array(partners[0], pa.int64()),
    }
    for j in range(1, J):
        cols[f"partner_{j}"] = pa.array(partners[j], pa.int64())
    idx_tbl = pa.table(cols)
    timing["arrays_s"] = time.perf_counter() - t3

    fam_fwd = np.array([FAMILY[fwd_fam[f][0]] for f in fw], np.uint8)
    seed_fwd = np.array([fwd_fam[f][1] for f in fw], np.int32)
    side = np.repeat(np.arange(-1, J), nf).astype(np.int8)       # -1 forward, j decoy set
    fam = np.tile(fam_fwd[seq_i], J + 1) + np.where(side >= 0, 2, 0).astype(np.uint8)
    tseq = np.array(tgt_of_fwd, dtype=object)[seq_i]
    ep = np.array([int.from_bytes(hashlib.blake2b(t.encode(), digest_size=8).digest(), "little")
                   for t in tgt_of_fwd], np.uint64)[seq_i]
    lab = pa.table({
        "pf_id": pa.array(pf_id[srt], pa.uint64()),
        "family": pa.array(fam[srt].astype(np.uint8)),
        "entrap_seed": pa.array(np.tile(seed_fwd[seq_i], J + 1)[srt], pa.int32()),
        "decoy_set": pa.array(side[srt], pa.int8()),
        "decoy_seed": pa.array(np.array([0] + dseeds, np.int32)[side[srt] + 1], pa.int32()),
        "is_decoy": pa.array(side[srt] >= 0),
        "is_entrapment": pa.array((fam[srt] % 2) == 1),
        "tdc_pair": pa.array(np.tile(r0["id_f"], J + 1)[srt], pa.uint64()),
        "entrap_pair": pa.array(np.tile(ep, J + 1)[srt], pa.uint64()),
        "fwd_seq": pa.array(np.tile(np.array(fw, dtype=object)[seq_i], J + 1)[srt].tolist(), pa.string()),
        "target_seq": pa.array(np.tile(tseq, J + 1)[srt].tolist(), pa.string()),
    })
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
        "index_version": INDEX_VERSION + "+multidecoy", "sage_commit": core.SAGE_COMMIT, "params": p,
        "fasta": [os.path.abspath(f) for f in fasta_paths], "content_hash": content_hash(fasta_paths, p),
        "n_targets": len(targets), "n_forward_peptides": len(fw), "decoy_sets": dec_stats,
        "n_forward_without_decoy": len(forward) - len(fw), "entrapment": ent_stats,
        "n_peptidoforms": int(N), "n_peptidoforms_by_family": {k: int(fams[v]) for k, v in FAMILY.items()},
        "n_decoy_sets": J, "timing_s": timing,
    }
    with open(os.path.join(out_dir, "meta.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    return meta


class L1Index:
    """Index side of a built L1 index (label-free)."""

    def __init__(self, path: str):
        self.path = path
        t = _read(os.path.join(path, "peptidoforms.arrow"))
        self.pf_id = t.column("pf_id").to_numpy()
        self.mass32 = t.column("mass32").to_numpy()
        self.seq_off = t.column("seq_off").to_numpy()
        self.length = t.column("length").to_numpy()
        self.mask = t.column("mask").to_numpy()
        self.n_var = t.column("n_var").to_numpy()
        self.missed = t.column("missed").to_numpy()
        self.partner = t.column("partner").to_numpy()
        # multi-decoy indexes (build_multi_decoy): partner per decoy set; [partner] for single-decoy indexes
        self.partners = [self.partner] + [t.column(f"partner_{j}").to_numpy()
                                          for j in range(1, 64) if f"partner_{j}" in t.column_names]
        self.residues = np.load(os.path.join(path, "residues.npy"))
        with open(os.path.join(path, "meta.json")) as fh:
            self.meta = json.load(fh)

    def __len__(self):
        return len(self.pf_id)

    def sequence(self, i: int) -> str:
        o = int(self.seq_off[i])
        return self.residues[o:o + int(self.length[i])].tobytes().decode()


def read_labels(path: str):
    """L4 only."""
    return _read(os.path.join(path, "labels.arrow"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--entrap-seeds", nargs="*", type=int, default=[])
    ap.add_argument("--decoy-seeds", nargs="*", type=int, default=[],
                    help="shuffled decoy sets (multi-replicate index for family top-2); default: pseudo-reverse")
    ap.add_argument("--workers", type=int, default=None)
    a = ap.parse_args(argv)
    if a.decoy_seeds:
        meta = build_multi_decoy(a.fasta, a.out, {"entrap_seeds": a.entrap_seeds, "decoy_seeds": a.decoy_seeds},
                                 workers=a.workers)
    else:
        meta = build(a.fasta, a.out, {"entrap_seeds": a.entrap_seeds}, workers=a.workers)
    print(json.dumps({k: v for k, v in meta.items() if k != "params"}, indent=1))


if __name__ == "__main__":
    main()
