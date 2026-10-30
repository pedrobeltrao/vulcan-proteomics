"""T1 driver, stage `db`: the identical external search databases for every method (DESIGN §1, §6.1-6.3, A1.9).

One multi-replicate L1 index (vulcan.index.build.build_multi_decoy) holds every forward family and every decoy set:
    T                targets (I/L collapsed tryptic peptides of the dataset FASTA)
    F                arm-C natural foreign peptides (neighbour-filtered by scripts/entrapment/build_dbs.py --arms C)
    E_s              paired shuffled entrapment, s in the entrapment seeds (dev 2001-2003; sealed 2004-2006)
    dec_j(.)         C-term-fixed hash(peptide, j, k) shuffled decoys of every forward peptide, j in DECOY_SEEDS 1001-1005
Replicate databases are row subsets of it (materialised as ordinary index directories, so every tool in the repo can
read them) and are written as peptide-level FASTAs for the comparators (one entry per sequence; decoys prefixed rev_):
    prod_j      T + dec_j(T)                              production search (sensitivity)
    P_j_s       T + E_s + dec_j(T) + dec_j(E_s)           arm P (paired), 5 x 3 replicates
    C_j         T + F + dec_j(T) + dec_j(F)               arm C (combined), 5 replicates
Every method searches exactly these peptide sets with the same fixed/variable modifications, so the external decoys
and entrapments are identical by construction. verify() checks §6.1 uniqueness on the final sequence sets.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import time
import zipfile

import numpy as np
import pyarrow as pa

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def seq_bytes(residues, off, ln, width: int = 64) -> np.ndarray:
    """Fixed-width byte strings (dtype S<width>) of index rows (each row owns its residues slice)."""
    off = np.asarray(off, np.int64)
    ln = np.asarray(ln, np.int64)
    assert ln.max(initial=0) <= width
    out = np.zeros((len(off), width), np.uint8)
    res = np.asarray(residues)
    for k in range(int(ln.max(initial=0))):
        m = ln > k
        out[m, k] = res[off[m] + k]
    return out.view(f"S{width}").ravel()


def replicate_names(decoy_seeds, entrap_seeds):
    out = []
    for j in decoy_seeds:
        out.append(f"prod_{j}")
        out += [f"P_{j}_{s}" for s in entrap_seeds]
        out.append(f"C_{j}")
    return out


def parse_rep(name: str):
    """-> (kind, decoy seed, entrapment seed or None)."""
    p = name.split("_")
    return p[0], int(p[1]), (int(p[2]) if len(p) > 2 else None)


def target_fasta(cfg: dict, work: str) -> str:
    out = os.path.join(work, "target.fasta")
    if os.path.exists(out):
        return out
    os.makedirs(work, exist_ok=True)
    with open(out + ".partial", "w") as oh:
        for f in cfg["fasta"]:
            p = os.path.join(ROOT, f)
            if p.endswith(".zip"):
                with zipfile.ZipFile(p) as z:
                    for n in sorted(z.namelist()):
                        if n.lower().endswith((".fasta", ".fa")):
                            txt = z.read(n).decode()
                            oh.write(txt if txt.endswith("\n") else txt + "\n")
            else:
                txt = open(p).read()
                oh.write(txt if txt.endswith("\n") else txt + "\n")
    os.replace(out + ".partial", out)
    return out


def foreign_peptides(cfg: dict, work: str, tfa: str, log) -> tuple[list[str], dict]:
    """Arm-C foreign peptides via the M1e builder (I/L, single-substitution, isobaric filters; §6.3)."""
    d = os.path.join(work, "armC")
    tsv = os.path.join(d, "ds.C.peptides.tsv")
    if not os.path.exists(tsv):
        os.makedirs(d, exist_ok=True)
        cmd = [sys.executable, os.path.join(ROOT, "scripts", "entrapment", "build_dbs.py"), "--name", "ds",
               "--target", tfa, "--foreign", *[os.path.join(ROOT, f) for f in cfg["foreign"]], "--arms", "C",
               "--out", d]
        log("armC: " + " ".join(cmd))
        with open(os.path.join(d, "build.log"), "w") as fh:
            subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.STDOUT, cwd=ROOT)
    out = []
    with open(tsv) as fh:
        next(fh)
        for line in fh:
            seq, kind = line.split("\t", 2)[:2]
            if kind == "E":
                out.append(seq)
    meta = json.load(open(os.path.join(d, "ds.C.json")))
    return out, meta


def build_full(cfg: dict, work: str, decoy_seeds, entrap_seeds, log) -> str:
    from vulcan.index.build import build_multi_decoy
    d = os.path.join(work, "index", "full")
    if os.path.exists(os.path.join(d, "meta.json")):
        return d
    tfa = target_fasta(cfg, work)
    foreign, cmeta = foreign_peptides(cfg, work, tfa, log)
    fpath = os.path.join(work, "armC", "foreign.txt")
    with open(fpath, "w") as fh:
        fh.write("\n".join(foreign) + "\n")
    t = time.perf_counter()
    meta = build_multi_decoy([tfa], d, {"entrap_seeds": list(entrap_seeds), "decoy_seeds": list(decoy_seeds),
                                        "foreign_sha256": sha256(fpath)}, foreign=foreign)
    meta["armC_builder"] = {k: v for k, v in cmeta.items() if k not in ("fasta", "fasta_sha256")}
    meta["build_wall_s"] = time.perf_counter() - t
    with open(os.path.join(d, "meta.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    log(f"full index built in {meta['build_wall_s']:.0f}s: {meta['n_peptidoforms']} peptidoforms")
    return d


def _read(path):
    with pa.OSFile(path, "rb") as src:
        return pa.ipc.open_file(src).read_all()


def _write(table, path):
    with pa.OSFile(path + ".partial", "wb") as sink, pa.ipc.new_file(sink, table.schema) as w:
        w.write_table(table)
    os.replace(path + ".partial", path)


def materialize(full: str, work: str, decoy_seeds, entrap_seeds, log) -> dict:
    """Write every replicate as an index directory (row subset, partner of its decoy set, shared residues) and as a
    peptide-level FASTA. Returns {rep: {"index": dir, "fasta": path, ...}}."""
    reps = replicate_names(decoy_seeds, entrap_seeds)
    out = {r: {"index": os.path.join(work, "index", r), "fasta": os.path.join(work, "fasta", r + ".fasta")}
           for r in reps}
    if all(os.path.exists(os.path.join(v["index"], "meta.json")) and os.path.exists(v["fasta"]) for v in out.values()):
        return out
    pf = _read(os.path.join(full, "peptidoforms.arrow"))
    lab = _read(os.path.join(full, "labels.arrow"))
    lab_schema = lab.schema
    # string columns of a large index exceed 2 GB: take() through large_string, cast back per replicate
    lab = lab.cast(pa.schema([pa.field(f.name, pa.large_string() if f.type == pa.string() else f.type)
                              for f in lab_schema]))
    meta = json.load(open(os.path.join(full, "meta.json")))
    residues = np.load(os.path.join(full, "residues.npy"), mmap_mode="r")
    fam = lab.column("family").to_numpy()
    eseed = lab.column("entrap_seed").to_numpy()
    dset = lab.column("decoy_set").to_numpy()
    seq_off = pf.column("seq_off").to_numpy()
    length = pf.column("length").to_numpy()
    fwdfam = np.where(fam % 2 == 0, 0, np.where(eseed == -1, -1, eseed))      # 0 T, -1 foreign, s entrapment
    partners = [pf.column("partner").to_numpy()] + [pf.column(f"partner_{j}").to_numpy()
                                                     for j in range(1, len(decoy_seeds))]
    dseeds = list(meta["params"]["decoy_seeds"])
    assert dseeds == list(decoy_seeds), (dseeds, decoy_seeds)
    os.makedirs(os.path.join(work, "fasta"), exist_ok=True)
    for r in reps:
        kind, j_seed, s = parse_rep(r)
        j = dseeds.index(j_seed)
        fams = {"prod": [0], "P": [0, s], "C": [0, -1]}[kind]
        keep = np.isin(fwdfam, fams) & np.isin(dset, [-1, j])
        rows = np.flatnonzero(keep)
        remap = np.full(len(fam), -1, np.int64)
        remap[rows] = np.arange(len(rows))
        pj = remap[partners[j][rows]]
        assert (pj >= 0).all(), f"{r}: partner outside the replicate"
        d = out[r]["index"]
        os.makedirs(d, exist_ok=True)
        cols = {n: pf.column(n).take(pa.array(rows)) for n in ("pf_id", "mass32", "seq_off", "length", "mask",
                                                              "n_var", "missed")}
        cols["partner"] = pa.array(pj, pa.int64())
        _write(pa.table(cols), os.path.join(d, "peptidoforms.arrow"))
        _write(lab.take(pa.array(rows)).cast(lab_schema), os.path.join(d, "labels.arrow"))
        for f in ("residues.npy", "proteins.arrow", "protein_idx.npy", "accessions.txt"):
            dst = os.path.join(d, f)
            if not os.path.lexists(dst):
                os.symlink(os.path.join(os.path.abspath(full), f), dst)
        m = {k: v for k, v in meta.items() if k not in ("timing_s",)}
        m.update({"replicate": r, "parent_index": os.path.abspath(full), "parent_content_hash": meta["content_hash"],
                  "n_peptidoforms": int(len(rows)), "decoy_seed": j_seed, "entrap_seed": s,
                  "forward_families": fams})
        with open(os.path.join(d, "meta.json") + ".partial", "w") as fh:
            json.dump(m, fh, indent=1)
        os.replace(os.path.join(d, "meta.json") + ".partial", os.path.join(d, "meta.json"))
        # FASTA: one entry per distinct sequence
        sb = seq_bytes(residues, seq_off[rows], length[rows])
        _, first = np.unique(sb, return_index=True)
        first = np.sort(first)
        rr = rows[first]
        with open(out[r]["fasta"] + ".partial", "w") as fh:
            for i, (x, b) in enumerate(zip(rr, sb[first])):
                ff = int(fwdfam[x])
                pre = ("rev_" if dset[x] >= 0 else "") + ("T" if ff == 0 else "F" if ff == -1 else "E")
                fh.write(f">{pre}_{i:08d}\n{b.decode()}\n")
        os.replace(out[r]["fasta"] + ".partial", out[r]["fasta"])
        out[r].update({"n_peptidoforms": int(len(rows)), "n_sequences": int(len(rr))})
        log(f"materialized {r}: {len(rows)} peptidoforms, {len(rr)} sequences")
    return out


def verify(full: str) -> dict:
    """§6.1 on the final sequence sets: no sequence belongs to two forward families, to a forward family and a decoy
    set, or to two decoy sets; I/L collapsed; every pair partner is a same-composition, C-term-fixed, different
    sequence."""
    pf = _read(os.path.join(full, "peptidoforms.arrow"))
    lab = _read(os.path.join(full, "labels.arrow"))
    residues = np.load(os.path.join(full, "residues.npy"), mmap_mode="r")
    seq_off = pf.column("seq_off").to_numpy()
    length = pf.column("length").to_numpy()
    fam = lab.column("family").to_numpy()
    eseed = lab.column("entrap_seed").to_numpy()
    dset = lab.column("decoy_set").to_numpy()
    sb = seq_bytes(residues, seq_off, length)
    dset = dset.astype(np.int64)
    kcode = np.where(dset >= 0, 100000 + dset, np.where(fam % 2 == 0, 0, eseed.astype(np.int64)))
    order = np.argsort(sb, kind="stable")
    s_sorted = sb[order]
    grp = np.r_[True, s_sorted[1:] != s_sorted[:-1]]
    gid = np.cumsum(grp) - 1
    ng = int(gid[-1]) + 1
    kmin = np.full(ng, np.iinfo(np.int64).max)
    kmax = np.full(ng, np.iinfo(np.int64).min)
    np.minimum.at(kmin, gid, kcode[order])
    np.maximum.at(kmax, gid, kcode[order])
    kg = kmin
    res = {"n_sequences": {}, "violations": {}}
    for k in np.unique(kg):
        name = (f"dec_set{int(k) - 100000}" if k >= 100000 else "T" if k == 0 else "F" if k == -1 else f"E{int(k)}")
        res["n_sequences"][name] = int((kg == k).sum())
    res["violations"]["sequence_in_two_families_or_decoy_sets"] = int((kmin != kmax).sum())
    res["violations"]["contains_I"] = int(sum(1 for x in s_sorted[grp] if b"I" in x))
    partners = [pf.column("partner").to_numpy()] + [pf.column(f"partner_{j}").to_numpy()
                                                     for j in range(1, 64) if f"partner_{j}" in pf.column_names]
    rng = np.random.default_rng(11)
    fw = np.flatnonzero(dset < 0)
    smp = rng.choice(fw, size=min(20000, len(fw)), replace=False)
    bad = 0
    for p in partners:
        for x in smp:
            a, b = sb[x].decode(), sb[p[x]].decode()
            bad += (sorted(a) != sorted(b)) or a[-1] != b[-1] or a == b
    res["violations"]["partner_not_cterm_fixed_shuffle"] = int(bad)
    res["partner_rows_checked"] = int(len(smp) * len(partners))
    res["ok"] = all(v == 0 for v in res["violations"].values())
    return res
