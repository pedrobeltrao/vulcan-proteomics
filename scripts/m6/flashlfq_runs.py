#!/usr/bin/env python3
"""FlashLFQ 2.1.5 comparator runs on D9 (timed under the bench lock, 32 threads, MBR on = FlashLFQ default).

    python scripts/m6/flashlfq_runs.py ours_hye|sage_hye|ours_all9

  ours_hye : our accepted PSMs (ids over the 6 HYE runs, static_gbm_pred, peptidoform q <= 1%) on the 6 HYE runs
  sage_hye : Sage 0.14.6 PSMs (label 1, spectrum_q and peptide_q <= 1%) on the same 6 runs
  ours_all9: our PSMs from the 9-run ID pass (6 HYE + 3 human-only) -> FlashLFQ MBR false-transfer check
FlashLFQ's bundled Thermo reader cannot open .raw on this Linux node (mzLib "Error opening RAW file!", with and
without ICU), so it reads mzML: the ProteoBench-provided mzMLs for the 6 HYE runs (the same files Sage reads) and
ThermoRawFileParser 2.0 (-f 2, vendor centroiding) mzMLs for the 3 human-only runs.  Input: FlashLFQ generic TSV.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

import pandas as pd  # noqa: E402

FL = os.path.join(C.ROOT, "tools", "flashlfq", "CMD.dll")
DOTNET = os.path.join(C.ROOT, "tools", "dotnet10", "dotnet")
SPECDIR = os.path.join(C.WORK, "mzml_local")   # FlashLFQ's Linux RawFileReader path fails ("Error opening RAW file")


def ours(ids_dir):
    ids = pd.read_parquet(os.path.join(ids_dir, "ids.parquet"))
    pep = pd.read_parquet(os.path.join(ids_dir, "peptides.parquet")).set_index("pf_id")
    a = ids[ids["pass_0.01"] & ~ids["is_decoy"]]
    ms = pep["modseq"].reindex(a["pf_id"].to_numpy()).to_numpy()
    return pd.DataFrame({"File Name": a["run"].to_numpy(),
                         "Base Sequence": [re.sub(r"\[.*?\]", "", s) for s in ms],
                         "Full Sequence": ms,
                         "Peptide Monoisotopic Mass": a["calcmass"].astype(float).to_numpy(),
                         "Scan Retention Time": a["rt"].to_numpy(),
                         "Precursor Charge": a["charge"].astype(int).to_numpy(),
                         "Protein Accession": pep["proteins"].reindex(a["pf_id"].to_numpy()).to_numpy()})


def sage():
    t = pd.read_csv(os.path.join(C.WORK, "sage", "D9", "results.sage.tsv"), sep="\t")
    t = t[(t["label"] == 1) & (t["spectrum_q"] <= 0.01) & (t["peptide_q"] <= 0.01)]
    return pd.DataFrame({"File Name": t["filename"].str.replace(".mzML", "", regex=False),
                         "Base Sequence": t["peptide"].str.replace(r"\[.*?\]", "", regex=True),
                         "Full Sequence": t["peptide"], "Peptide Monoisotopic Mass": t["calcmass"],
                         "Scan Retention Time": t["rt"], "Precursor Charge": t["charge"],
                         "Protein Accession": t["proteins"]})


def main():
    kind = sys.argv[1]
    out = os.path.join(C.WORK, "flashlfq", kind)
    os.makedirs(out, exist_ok=True)
    if kind == "ours_hye":
        df = ours(os.path.join(C.WORK, "ids", "d9_hye"))
    elif kind == "ours_all9":
        df = ours(os.path.join(C.WORK, "ids", "d9_all"))
    else:
        df = sage()
    idt = os.path.join(out, "identifications.tsv")
    df.to_csv(idt, sep="\t", index=False)
    rep = os.path.join(out, "mzml")
    os.makedirs(rep, exist_ok=True)
    for f in sorted(set(df["File Name"])):
        dst = os.path.join(rep, f + ".mzML")
        if not os.path.exists(dst):
            os.symlink(os.path.join(SPECDIR, f + ".mzML"), dst)
    env = dict(os.environ, DOTNET_ROOT=os.path.dirname(DOTNET), DOTNET_SYSTEM_GLOBALIZATION_INVARIANT="1")
    cmd = [sys.executable, os.path.join(C.ROOT, "scripts", "baselines", "timed.py"), "--label", f"flashlfq_{kind}_t32",
           "--log", os.path.join(C.RES, "timings.jsonl"), "--lock", "--threads", "32",
           "--stdout", os.path.join(C.WORK, "logs", f"flashlfq_{kind}.log"), "--",
           DOTNET, FL, "--idt", idt, "--rep", rep, "--out", out, "--thr", "32", "--mbr", "true",
           "--sil", "true"]
    sys.exit(subprocess.call(cmd, env=env))


if __name__ == "__main__":
    main()
