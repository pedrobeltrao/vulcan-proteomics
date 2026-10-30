"""Run Sage 0.14.7 on exactly our spectra and our index (E2 comparator and scoring oracle).

* spectra: the store's MS2 scans written as MGF (TITLE = row index in the run's MS2 table, f32-exact m/z text,
  CHARGE only when the instrument reported one -> Sage then tries precursor_charge range like we do).
* index: one FASTA entry per peptide sequence of our L1 index (targets "T_", entrapments "E_", decoys "rev_"),
  searched without cleavage (cleave_at "$"), generate_decoys false, same mods/limits -> the same peptidoforms.
  (The FASTA is built from labels.arrow: this is L4/benchmark plumbing, not part of the search path.)
"""
from __future__ import annotations

import json
import os
import re

import numpy as np
import pyarrow as pa
import pyarrow.csv as pacsv

from vulcan.index.build import _read

SAGE_BIN = os.path.join(os.path.dirname(__file__), "..", "..", "tools", "bin", "sage")


def write_mgf(sp, path: str) -> int:
    """All MS2 spectra of `sp` (vulcan.gpu.engine.Spectra). Returns bytes written."""
    lines = []
    off = sp.offsets
    mzs = np.char.mod("%.9g", sp.mz.astype(np.float32))
    its = np.char.mod("%.9g", sp.intensity.astype(np.float32))
    pk = np.char.add(np.char.add(mzs, " "), its)
    pmz = np.char.mod("%.9g", sp.precursor_mz.astype(np.float32))
    for i in range(len(sp)):
        lines.append("BEGIN IONS")
        lines.append(f"TITLE={i}")
        lines.append(f"RTINSECONDS={sp.rt[i] * 60.0:.6f}")
        lines.append(f"PEPMASS={pmz[i]}")
        if sp.charge[i] > 0:
            lines.append(f"CHARGE={int(sp.charge[i])}+")
        lines.extend(pk[off[i]:off[i + 1]].tolist())
        lines.append("END IONS")
    txt = "\n".join(lines) + "\n"
    with open(path + ".partial", "w") as fh:
        fh.write(txt)
    os.replace(path + ".partial", path)
    return len(txt)


def index_sequences(index_dir: str):
    """(sequence, is_decoy, is_entrapment) for every distinct sequence of the index."""
    from vulcan.index.build import L1Index
    l1 = L1Index(index_dir)
    lab = _read(os.path.join(index_dir, "labels.arrow"))
    seqs = np.array([l1.sequence(i) for i in range(len(l1))], dtype=object)
    dec = lab.column("is_decoy").to_numpy(zero_copy_only=False)
    ent = lab.column("is_entrapment").to_numpy(zero_copy_only=False)
    _, first = np.unique(seqs, return_index=True)
    return seqs[first], dec[first], ent[first]


def write_peptide_fasta(index_dir: str, path: str) -> dict:
    seqs, dec, ent = index_sequences(index_dir)
    with open(path + ".partial", "w") as fh:
        for i, (s, d, e) in enumerate(zip(seqs, dec, ent)):
            pre = ("rev_" if d else "") + ("E_" if e else "T_")
            fh.write(f">{pre}{i:08d}\n{s}\n")
    os.replace(path + ".partial", path)
    return {"n_sequences": int(len(seqs)), "n_decoy": int(dec.sum()), "n_entrapment_fwd": int((ent & ~dec).sum())}


def sage_config(fasta: str, mgfs: list[str], outdir: str, p, static_mods: dict, var_mods: dict, max_var: int,
                min_mass=500.0, max_mass=5000.0) -> dict:
    """Sage JSON mirroring vulcan.gpu.engine.SearchParams `p` (peptide-level DB, no cleavage)."""
    return {
        "database": {
            "bucket_size": 8192,
            "enzyme": {"missed_cleavages": 0, "min_len": 7, "max_len": 50, "cleave_at": "$", "restrict": None,
                       "c_terminal": True},
            "peptide_min_mass": min_mass, "peptide_max_mass": max_mass,
            "fragment_min_mz": p.fragment_min_mz, "fragment_max_mz": p.fragment_max_mz,
            "ion_kinds": ["b", "y"], "min_ion_index": p.min_ion_index,
            "static_mods": static_mods, "variable_mods": {k: [v] for k, v in var_mods.items()},
            "max_variable_mods": max_var, "decoy_tag": "rev_", "generate_decoys": False, "fasta": fasta,
        },
        "precursor_tol": {"ppm": list(p.precursor_ppm)},
        "fragment_tol": {"ppm": list(p.fragment_ppm)},
        "precursor_charge": list(p.precursor_charge),
        "isotope_errors": list(p.isotope_errors),
        "deisotope": False, "chimera": False, "wide_window": False, "predict_rt": False,
        "min_peaks": p.min_peaks, "max_peaks": p.max_peaks, "min_matched_peaks": p.min_matched_peaks,
        "max_fragment_charge": p.max_fragment_charge, "report_psms": p.report_psms,
        "output_directory": outdir, "mzml_paths": mgfs,
    }


_MOD = re.compile(r"\[[^\]]*\]")


def parse_peptide(s: str, var_mods: dict) -> tuple[str, int]:
    """Sage peptide string -> (plain sequence, variable-mod position mask)."""
    seq, mask, k = [], 0, 0
    i = 0
    while i < len(s):
        c = s[i]
        i += 1
        if i < len(s) and s[i] == "[":
            j = s.index("]", i)
            m = float(s[i + 1:j])
            i = j + 1
            if c in var_mods and abs(m - var_mods[c]) < 1e-3:
                mask |= 1 << k
        seq.append(c)
        k += 1
    return "".join(seq), mask


def read_results(outdir: str) -> pa.Table:
    t = pacsv.read_csv(os.path.join(outdir, "results.sage.tsv"), parse_options=pacsv.ParseOptions(delimiter="\t"),
                       convert_options=pacsv.ConvertOptions(column_types={"scannr": pa.string()}))
    return t


def write_config(cfg: dict, path: str):
    with open(path, "w") as fh:
        json.dump(cfg, fh, indent=1)
