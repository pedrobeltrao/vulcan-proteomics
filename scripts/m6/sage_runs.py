#!/usr/bin/env python3
"""Sage comparator runs for M6 (configs written to results/.../configs, timed under the bench lock).

  D9: Sage LFQ (built in, defaults) on the 6 ProteoBench DDA-Astral HYE mzMLs, HYE FASTA, same search space as ours
      (trypsin KR/P, 2 missed, 7-50, CAM static, Met-ox var <=2, 10 ppm precursor, 20 ppm fragment, z 2-4, iso 0-2).
  D4: Sage TMT11 SPS-MS3 (level 3) on the 12 PXD016766 mzMLs of M1e, human+contaminant FASTA, TMT K/N-term static,
      0.5 Da fragment (IT-CID).
    python scripts/m6/sage_runs.py D9|D4 [--threads N] [--dry]
"""
import argparse, json, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

SAGE = os.path.join(C.ROOT, "tools", "bin", "sage")


def cfg(ds, out):
    base = {
        "database": {"bucket_size": 8192,
                     "enzyme": {"missed_cleavages": 2, "min_len": 7, "max_len": 50, "cleave_at": "KR",
                                "restrict": "P", "c_terminal": True},
                     "peptide_min_mass": 500.0, "peptide_max_mass": 5000.0, "ion_kinds": ["b", "y"],
                     "min_ion_index": 2, "static_mods": {"C": 57.021464}, "variable_mods": {"M": [15.9949]},
                     "max_variable_mods": 2, "decoy_tag": "rev_", "generate_decoys": True},
        "precursor_tol": {"ppm": [-10.0, 10.0]},
        "precursor_charge": [2, 4], "isotope_errors": [0, 2], "deisotope": False, "chimera": False,
        "wide_window": False, "predict_rt": True, "min_peaks": 15, "max_peaks": 150, "min_matched_peaks": 4,
        "max_fragment_charge": 1, "report_psms": 1, "output_directory": out,
    }
    if ds == "D9":
        base["database"]["fasta"] = C.D9_FASTA
        base["fragment_tol"] = {"ppm": [-20.0, 20.0]}
        base["quant"] = {"lfq": True}
        a, _ = C.d9_raws()
        base["mzml_paths"] = [r[:-4] + ".mzML" for r in a]
    else:
        base["database"]["fasta"] = C.D4_FASTA
        base["database"]["bucket_size"] = 32768
        base["database"]["static_mods"].update({"K": C.TMT, "^": C.TMT})
        base["fragment_tol"] = {"da": [-0.5, 0.5]}
        base["quant"] = {"tmt": "Tmt11", "tmt_settings": {"level": 3, "sn": False}}
        base["mzml_paths"] = sorted(os.path.join("/tmp/vulcan_m1e/mzml", os.path.basename(r)[:-4] + ".mzML")
                                    for r in C.d4_raws())
    return base


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ds")
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    out = os.path.join(C.WORK, "sage", a.ds)
    os.makedirs(out, exist_ok=True)
    c = cfg(a.ds, out)
    cd = os.path.join(C.RES, "configs")
    os.makedirs(cd, exist_ok=True)
    p = os.path.join(cd, f"sage_{a.ds}_t{a.threads}.json")
    json.dump(c, open(p, "w"), indent=1)
    if a.dry:
        print(p); return
    cmd = [sys.executable, os.path.join(C.ROOT, "scripts", "baselines", "timed.py"), "--label", f"sage_{a.ds}_t{a.threads}",
           "--log", os.path.join(C.RES, "timings.jsonl"), "--lock", "--threads", str(a.threads),
           "--stdout", os.path.join(C.WORK, "logs", f"sage_{a.ds}.log"), "--",
           SAGE, p, "--batch-size", str(min(len(c["mzml_paths"]), a.threads)), "--disable-telemetry-i-dont-want-to-improve-sage"]
    env = dict(os.environ, RAYON_NUM_THREADS=str(a.threads))
    sys.exit(subprocess.call(cmd, env=env))


if __name__ == "__main__":
    main()
