#!/usr/bin/env python3
"""M1e baseline searches: Sage (+LDA, +mokapot) and Comet (+mokapot) on D0/D1(/D4) with and without entrapment.

    python scripts/baselines/run_searches.py [--datasets D0 D1] [--tools sage comet] [--dbs prod target P2001 ...]
                                             [--lock-dbs prod P2001] [--force]

Search databases (built by scripts/entrapment/run_all.sh):
  prod     protein FASTA (target + contaminants), trypsin, engine-native decoys
           (Sage: generate_decoys, reversed peptides; Comet: decoy_search = 1). Sensitivity numbers come from here.
  target   peptide-level DB, T + dec(T) only (control for the peptide-level set-up)
  P<seed>  peptide-level DB, T + E_seed + dec(T) + dec(E_seed)  (arm P, FDRBench, r = 1)
  C        peptide-level DB, T + foreign + decoys of both     (arm C)
Peptide-level DBs are searched without further cleavage (Sage cleave_at "$"; Comet enzyme 11 No_cut), with the
same length (7-50) and mass (500-5000 Da) limits; external decoys carry the "rev_" prefix.

Every command runs through scripts/baselines/timed.py (GNU time -v; threads = vulcan.hw.n_threads(), i.e. affinity or $VULCAN_THREADS). Runs on DBs listed in --lock-dbs
hold the node-wide benchmark lock; the others are correctness runs (timed, but flagged unlocked).
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from vulcan import hw  # noqa: E402

THREADS = hw.n_threads()  # affinity mask or $VULCAN_THREADS (=4 for M1-continuity runs)
WORK = os.environ.get("VULCAN_M1E_WORK", "/tmp/vulcan_m1e")
N_AFF = len(os.sched_getaffinity(0))
# search outputs on local disk; the full-budget pass (threads == affinity) uses cache/, other thread counts cache_t<N>/
CACHE = os.path.join(WORK, "cache" if THREADS == N_AFF else f"cache_t{THREADS}")
TTAG = f"@t{THREADS}"                         # label suffix: timings of different thread budgets coexist in timings.jsonl
ENT = os.path.join(ROOT, "data", "cache", "M1e", "entrapment")   # db/ and inputs/ are symlinks into $WORK
RES = os.path.join(ROOT, "results", "2026-10-27_M1e_baselines")
TIMINGS = os.path.join(RES, "timings.jsonl")
TIMED = [sys.executable, os.path.join(ROOT, "scripts", "baselines", "timed.py")]
SAGE = os.path.join(ROOT, "tools", "bin", "sage")
COMET = os.path.join(ROOT, "conda_envs", "comet", "bin", "comet")
MOKAPOT = os.path.join(os.environ.get("VULCAN_MOKAPOT_ENV", os.path.join(WORK, "envs", "mokapot")), "bin", "mokapot")  # see setup_mokapot_env.sh

# Instrument-specific tolerances (see summary.md "Search settings"; checked on a wide-tolerance dev search).
DATASETS = {
    "D0": {"accession": "PXD001468", "family": "human_contam",
           "files": ["b1906_293T_proteinID_01A_QE3_122212", "b1937_293T_proteinID_01B_QE3_122212"],
           "instrument": "Q Exactive, HCD, Orbitrap MS2",
           "prec_ppm": 10.0, "frag_ppm": 20.0, "comet_frag_bin_tol": 0.02, "comet_frag_bin_offset": 0.0,
           "comet_theo_ions": 0},
    "D1": {"accession": "PXD001077", "family": "pfu_contam", "files": ["Velos005137"],
           "instrument": "LTQ Orbitrap Velos, HCD, Orbitrap MS2 (R=7500)",
           "prec_ppm": 10.0, "frag_ppm": 20.0, "comet_frag_bin_tol": 0.02, "comet_frag_bin_offset": 0.0,
           "comet_theo_ions": 0},
    # D4: TMT11 SPS-MS3; MS2 = CID in the ion trap -> low-resolution fragment settings; TMT static on K and N-term.
    "D4": {"accession": "PXD016766", "family": "human_contam",
           # 2 of the 12 files (first of plate rows A and B): Comet on the peptide-level entrapment DBs costs 5-10 min per
           # file per DB at 32 threads (8x that at 4 threads), 12 files x 5 DBs x 2 thread budgets is days. Documented deviation.
           "files": ["dq_00095_11cell_90min_sps_A1", "dq_00101_11cell_90min_sps_B1"],
           "instrument": "Orbitrap Eclipse, TMT11 SPS-MS3, MS2 CID ion trap",
           "prec_ppm": 10.0, "frag_da": 0.5, "comet_frag_bin_tol": 1.0005, "comet_frag_bin_offset": 0.4,
           "comet_theo_ions": 1, "tmt": 229.162932, "bucket_size": 32768},
}
SEEDS = [2001, 2002, 2003]
ALL_DBS = ["prod", "target"] + [f"P{s}" for s in SEEDS] + ["C"]


def db_fasta(ds, db):
    fam = DATASETS[ds]["family"]
    if db == "prod":
        return os.path.join(ENT, "inputs", f"{fam}.fasta")
    return os.path.join(ENT, "db", f"{fam}.{db}.fasta")


def mzmls(ds):
    return [os.path.join(WORK, "mzml", f"{f}.mzML") for f in DATASETS[ds]["files"]]


def sage_config(ds, db, outdir):
    d = DATASETS[ds]
    peptide_level = db != "prod"
    enzyme = ({"missed_cleavages": 0, "min_len": 7, "max_len": 50, "cleave_at": "$", "restrict": None,
               "c_terminal": True} if peptide_level else
              {"missed_cleavages": 2, "min_len": 7, "max_len": 50, "cleave_at": "KR", "restrict": "P",
               "c_terminal": True})
    return {
        "database": {
            "bucket_size": d.get("bucket_size", 8192),
            "enzyme": enzyme,
            "peptide_min_mass": 500.0, "peptide_max_mass": 5000.0,
            "ion_kinds": ["b", "y"], "min_ion_index": 2,
            "static_mods": {"C": 57.021464, **({"K": d["tmt"], "^": d["tmt"]} if "tmt" in d else {})},
            "variable_mods": {"M": [15.9949]},
            "max_variable_mods": 2,
            "decoy_tag": "rev_",
            "generate_decoys": not peptide_level,
            "fasta": db_fasta(ds, db),
        },
        "precursor_tol": {"ppm": [-d["prec_ppm"], d["prec_ppm"]]},
        "fragment_tol": ({"da": [-d["frag_da"], d["frag_da"]]} if "frag_da" in d else
                         {"ppm": [-d["frag_ppm"], d["frag_ppm"]]}),
        "precursor_charge": [2, 4],
        "isotope_errors": [0, 2],
        "deisotope": False, "chimera": False, "wide_window": False, "predict_rt": True,
        "min_peaks": 15, "max_peaks": 150, "min_matched_peaks": 4, "max_fragment_charge": 1,
        "report_psms": 1,
        "output_directory": outdir,
        "mzml_paths": mzmls(ds),
    }


COMET_TEMPLATE = """# comet_version 2026.01 rev. 1 (e4f767c)
# written by scripts/baselines/run_searches.py (M1e)
database_name = {db}
decoy_search = {decoy_search}
decoy_prefix = {decoy_prefix}
num_threads = {threads}
peptide_mass_tolerance_upper = {prec}
peptide_mass_tolerance_lower = -{prec}
peptide_mass_units = 2
precursor_tolerance_type = 1
isotope_error = 2
search_enzyme_number = {enzyme}
search_enzyme2_number = 0
sample_enzyme_number = 1
num_enzyme_termini = 2
allowed_missed_cleavage = {missed}
variable_mod01 = 15.9949 M 0 2 -1 0 0 0.0
variable_mod02 = 0.0 X 0 3 -1 0 0 0.0
variable_mod03 = 0.0 X 0 3 -1 0 0 0.0
variable_mod04 = 0.0 X 0 3 -1 0 0 0.0
variable_mod05 = 0.0 X 0 3 -1 0 0 0.0
max_variable_mods_in_peptide = 2
require_variable_mod = 0
fragment_bin_tol = {fbt}
fragment_bin_offset = {fbo}
theoretical_fragment_ions = {theo}
use_A_ions = 0
use_B_ions = 1
use_C_ions = 0
use_X_ions = 0
use_Y_ions = 1
use_Z_ions = 0
use_Z1_ions = 0
use_NL_ions = 0
output_sqtfile = 0
output_txtfile = 0
output_pepxmlfile = 0
output_mzidentmlfile = 0
output_percolatorfile = 1
num_output_lines = 1
scan_range = 0 0
precursor_charge = 0 0
override_charge = 0
ms_level = 2
activation_method = ALL
digest_mass_range = 500.0 5000.0
peptide_length_range = 7 50
max_duplicate_proteins = 10
max_fragment_charge = 3
min_precursor_charge = 1
max_precursor_charge = 6
clip_nterm_methionine = 0
spectrum_batch_size = 15000
equal_I_and_L = 1
mass_offsets =
minimum_peaks = 10
minimum_intensity = 0
remove_precursor_peak = 0
remove_precursor_tolerance = 1.5
clear_mz_range = 0.0 0.0
percentage_base_peak = 0.0
add_Cterm_peptide = 0.0
add_Cterm_protein = 0.0
add_Nterm_protein = 0.0
add_C_cysteine = 57.021464
add_K_lysine = {tmt}
add_Nterm_peptide = {tmt}

[COMET_ENZYME_INFO]
0.  Cut_everywhere         0      -           -
1.  Trypsin                1      KR          P
2.  Trypsin/P              1      KR          -
3.  Lys_C                  1      K           P
4.  Lys_N                  0      K           -
5.  Arg_C                  1      R           P
6.  Asp_N                  0      DN          -
7.  CNBr                   1      M           -
8.  Asp-N_ambic            1      DE          -
9.  PepsinA                1      FL          -
10. Chymotrypsin           1      FWYL        P
11. No_cut                 1      @           @
"""


def comet_params(ds, db):
    d = DATASETS[ds]
    peptide_level = db != "prod"
    return COMET_TEMPLATE.format(
        db=db_fasta(ds, db), decoy_search=0 if peptide_level else 1,
        decoy_prefix="rev_" if peptide_level else "DECOY_", prec=d["prec_ppm"],
        enzyme=11 if peptide_level else 1, missed=0 if peptide_level else 2,
        fbt=d["comet_frag_bin_tol"], fbo=d["comet_frag_bin_offset"], theo=d["comet_theo_ions"],
        tmt=d.get("tmt", 0.0), threads=THREADS)


def timed(label, cmd, log, lock, meta, cwd=None):
    meta = {**meta, "threads": THREADS}
    full = TIMED + ["--label", label + TTAG, "--threads", str(THREADS), "--log", TIMINGS, "--stdout", log, "--meta", json.dumps(meta)] + \
        (["--lock"] if lock else []) + ["--"] + cmd
    print("RUN", label, "(lock)" if lock else "", flush=True)
    rc = subprocess.call(full, cwd=cwd)
    if rc != 0:
        raise SystemExit(f"{label} failed rc={rc}; see {log}")


def run_mokapot(tool, ds, db, pins, lock, force):
    out = os.path.join(CACHE, f"{tool}_mokapot", ds, db)
    if os.path.exists(os.path.join(out, "mokapot.peptides.txt")) and not force:
        return
    os.makedirs(out, exist_ok=True)
    decoy_prefix = "rev_" if (db != "prod" or tool == "sage") else "DECOY_"
    cmd = [MOKAPOT, "--dest_dir", out, "--seed", "11", "--max_workers", "1", "--keep_decoys",
           "--decoy_prefix", decoy_prefix] + (["--aggregate"] if len(pins) > 1 else []) + [*pins]
    timed(f"mokapot:{tool}:{ds}:{db}", cmd, os.path.join(out, "mokapot.log"), lock,
          {"stage": "rescoring", "tool": f"{tool}+mokapot", "dataset": ds, "db": db, "n_files": len(pins)})


def run_sage(ds, db, lock, force):
    out = os.path.join(CACHE, "sage", ds, db)
    if os.path.exists(os.path.join(out, "results.sage.tsv")) and not force:
        return out
    os.makedirs(out, exist_ok=True)
    cfg = sage_config(ds, db, out)
    cpath = os.path.join(RES, "configs", f"sage_{ds}_{db}_t{THREADS}.json")
    with open(cpath, "w") as fh:
        json.dump(cfg, fh, indent=1)
    env_cmd = ["env", "RUST_LOG=info", f"RAYON_NUM_THREADS={THREADS}", SAGE, cpath, "--write-pin",
               "--batch-size", str(max(1, min(len(cfg["mzml_paths"]), THREADS))), "--disable-telemetry-i-dont-want-to-improve-sage"]
    timed(f"sage:{ds}:{db}", env_cmd, os.path.join(out, "sage.log"), lock,
          {"stage": "search", "tool": "sage", "dataset": ds, "db": db, "n_files": len(cfg["mzml_paths"]),
           "entrapment": db not in ("prod", "target")})
    return out


def run_comet(ds, db, lock, force):
    out = os.path.join(CACHE, "comet", ds, db)
    os.makedirs(out, exist_ok=True)
    ppath = os.path.join(RES, "configs", f"comet_{ds}_{db}_t{THREADS}.params")
    with open(ppath, "w") as fh:
        fh.write(comet_params(ds, db))
    pins = []
    for mz in mzmls(ds):
        base = os.path.join(out, os.path.basename(mz)[:-5])
        pins.append(base + ".pin")
        if os.path.exists(base + ".pin") and not force:
            continue
        timed(f"comet:{ds}:{db}:{os.path.basename(base)}", [COMET, f"-P{ppath}", f"-N{base}", mz],
              base + ".comet.log", lock,
              {"stage": "search", "tool": "comet", "dataset": ds, "db": db, "file": os.path.basename(base),
               "entrapment": db not in ("prod", "target")}, cwd=out)
    return out, pins


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=["D0", "D1"])
    ap.add_argument("--tools", nargs="+", default=["sage", "comet"])
    ap.add_argument("--dbs", nargs="+", default=ALL_DBS)
    ap.add_argument("--lock-dbs", nargs="*", default=["prod", "P2001"])
    ap.add_argument("--no-mokapot", action="store_true")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    os.makedirs(os.path.join(RES, "configs"), exist_ok=True)
    for ds in a.datasets:
        for mz in mzmls(ds):
            assert os.path.exists(mz), f"missing {mz}; run scripts/baselines/01_convert.sh"
        for db in a.dbs:
            if not os.path.exists(db_fasta(ds, db)):
                print("skip (no DB)", ds, db)
                continue
            lock = db in a.lock_dbs
            if "sage" in a.tools:
                out = run_sage(ds, db, lock, a.force)
                if not a.no_mokapot:
                    pins = sorted(glob.glob(os.path.join(out, "*.pin")))
                    run_mokapot("sage", ds, db, pins, lock, a.force)
            if "comet" in a.tools:
                out, pins = run_comet(ds, db, lock, a.force)
                if not a.no_mokapot:
                    run_mokapot("comet", ds, db, pins, lock, a.force)


if __name__ == "__main__":
    main()
