#!/usr/bin/env python3
"""Sage 0.14.7 semi-specific comparator on D7 (same spectra: MGFs written from our stores, as in E2/M4).

Sage's own semi-specific mode: FASTA = yeast + contaminants + P37840, trypsin KR/P, 2 missed, length 7-50,
semi_enzymatic true, generate_decoys true (Sage reverses proteins, i.e. decoy *proteins* digested semi-specifically),
CAM static, Met-ox variable (<= 2), our tolerances (precursor +-10 ppm, fragment +-20 ppm, isotope errors 0..2),
predict_rt on (strongest Sage configuration, M3b review), LDA rescoring, Sage's own q-values.
One Sage run per digest group (LiP 6 files, TrP 6 files), timed with scripts/baselines/timed.py --lock.
    python scripts/m5lc/sage_lip.py mgf            # write MGFs
    python scripts/m5lc/sage_lip.py run LiP|TrP [threads]
"""
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

from vulcan.gpu.host import SearchParams  # noqa: E402


def mgf_path(run):
    return os.path.join(C.WORK, "mgf", run + ".mgf")


def write_mgfs():
    from vulcan.gpu import sage_bridge as SB
    from vulcan.gpu.engine import load_store
    os.makedirs(os.path.join(C.WORK, "mgf"), exist_ok=True)
    for run in sorted(C.D7_DESIGN):
        if not os.path.exists(mgf_path(run)):
            SB.write_mgf(load_store(C.store(run)), mgf_path(run))
            print("mgf", run, flush=True)


def config(group, outdir):
    p = SearchParams()
    runs = [r for r in sorted(C.D7_DESIGN) if C.D7_DESIGN[r][0] == group]
    return {
        "database": {
            "bucket_size": 8192,
            "enzyme": {"missed_cleavages": 2, "min_len": 7, "max_len": 50, "cleave_at": "KR", "restrict": "P",
                       "c_terminal": True, "semi_enzymatic": True},
            "peptide_min_mass": 500.0, "peptide_max_mass": 5000.0,
            "fragment_min_mz": p.fragment_min_mz, "fragment_max_mz": p.fragment_max_mz,
            "ion_kinds": ["b", "y"], "min_ion_index": p.min_ion_index,
            "static_mods": C.STATIC, "variable_mods": {k: [v] for k, v in C.VAR.items()}, "max_variable_mods": 2,
            "decoy_tag": "rev_", "generate_decoys": True, "fasta": C.LIP_FASTA,
        },
        "precursor_tol": {"ppm": list(p.precursor_ppm)}, "fragment_tol": {"ppm": list(p.fragment_ppm)},
        "precursor_charge": list(p.precursor_charge), "isotope_errors": list(p.isotope_errors),
        "deisotope": False, "chimera": False, "wide_window": False, "predict_rt": True,
        "min_peaks": p.min_peaks, "max_peaks": p.max_peaks, "min_matched_peaks": p.min_matched_peaks,
        "max_fragment_charge": p.max_fragment_charge, "report_psms": 1,
        "output_directory": outdir, "mzml_paths": [mgf_path(r) for r in runs],
    }


def run(group, threads):
    out = os.path.join(C.WORK, "sage", f"semi_{group}_t{threads}")
    os.makedirs(out, exist_ok=True)
    json.dump(config(group, out), open(os.path.join(out, "config.json"), "w"), indent=1)
    os.makedirs(C.RES, exist_ok=True)
    subprocess.run([sys.executable, os.path.join(C.ROOT, "scripts/baselines/timed.py"), "--label",
                    f"M5_sage_semi_D7_{group}_t{threads}", "--log", os.path.join(C.RES, "timings.jsonl"), "--lock",
                    "--threads", str(threads), "--stdout", os.path.join(out, "sage.log"),
                    "--meta", json.dumps({"dataset": "D7", "group": group, "tool": "sage-0.14.7", "mode": "semi",
                                          "threads": threads, "outdir": out}), "--",
                    os.path.join(C.ROOT, "tools/bin/sage"), os.path.join(out, "config.json"),
                    "--disable-telemetry-i-dont-want-to-improve-sage"], check=True, cwd=C.ROOT)


if __name__ == "__main__":
    if sys.argv[1] == "mgf":
        write_mgfs()
    else:
        run(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 32)
