"""T1 driver, stage `comparators` (DESIGN §1 "the bar", §6.0 A1.9: the comparators at their strongest).

All comparators search the SAME peptide-level replicate FASTAs as Vulcan (bench/t1/dbs.py; external shuffled decoys
prefixed rev_, no cleavage inside an entry, same length/mass limits, CAM fixed, Met-ox variable <= 2, same
precursor/fragment tolerances), from ThermoRawFileParser mzML (vendor centroids), like the M1e baselines.

  sage            Sage 0.14.7 with its built-in LDA, predict_rt = true (checked to take effect: rt_check())
  sage+mokapot    mokapot 0.10.0 on Sage's PIN (includes Sage's RT-model features)
  sage+pred+mokapot  Sage PSMs + the SAME AlphaPeptDeep-derived L3 features Vulcan uses (pred_sa, pred_sa_sqrt,
                  pred_pearson, pred_entropy, pred_frac_obs, pred_top6, pred_obs_tic, rt_abs_resid, rt_resid_z;
                  vulcan.rescore.features on Sage's top-1 PSMs, label-blind) appended to Sage's PIN -> mokapot.
                  An MS2Rescore/Oktoberfest-style predicted-feature rescorer on the comparator's own PSMs.
  comet+mokapot   Comet 2026.01 (spectrum_batch_size 0: output identical to the M1e batching, ~2x faster) -> mokapot
                  --aggregate over the per-file PINs.
mokapot: --seed 11 --max_workers 1 --keep_decoys --decoy_prefix rev_ (as M1e).
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SAGE = os.path.join(ROOT, "tools", "bin", "sage")
COMET = os.path.join(ROOT, "conda_envs", "comet", "bin", "comet")
TRFP = os.path.join(ROOT, "conda_envs", "trfp", "bin", "ThermoRawFileParser")
MOKAPOT_ENV = os.environ.get("VULCAN_MOKAPOT_ENV", "/tmp/vulcan_m1e/envs/mokapot")
METHODS = ["sage", "sage+mokapot", "sage+pred+mokapot", "comet+mokapot"]
PRED_PIN_FEATURES = ["pred_sa", "pred_sa_sqrt", "pred_pearson", "pred_entropy", "pred_frac_obs", "pred_top6",
                     "pred_obs_tic", "rt_abs_resid", "rt_resid_z"]


def timed_run(cmd, log_path, cwd=None, env=None, lock=False) -> dict:
    """Run a command under /usr/bin/time -v; optional flock on the node-wide benchmark lock."""
    tf = log_path + ".time"
    full = ["/usr/bin/time", "-v", "-o", tf] + list(cmd)
    if lock:
        full = ["flock", os.environ.get("VULCAN_BENCH_LOCK", "/tmp/vulcan_bench.lock")] + full
    t0 = time.time()
    with open(log_path, "w") as fh:
        rc = subprocess.call(full, stdout=fh, stderr=subprocess.STDOUT, cwd=cwd, env=env)
    rec = {"wall_s": time.time() - t0, "rc": rc, "loadavg_end": os.getloadavg()}
    try:
        txt = open(tf).read()
        u = float(re.search(r"User time \(seconds\): ([\d.]+)", txt).group(1))
        s = float(re.search(r"System time \(seconds\): ([\d.]+)", txt).group(1))
        rec["core_s"] = u + s
        rec["max_rss_mb"] = float(re.search(r"Maximum resident set size \(kbytes\): (\d+)", txt).group(1)) / 1024
    except Exception:  # noqa: BLE001
        pass
    if rc != 0:
        raise RuntimeError(f"{cmd[0]} failed rc={rc}; see {log_path}")
    return rec


# ------------------------------------------------------------------------------------------------ mzML
def ensure_mzml(raws, work, mzml_cache=None, threads=4, log=print) -> dict:
    """RAW -> indexed mzML (TRFP 2.0, vendor peak picking, -f 2), staged on local disk. mzml_cache: optional directory
    of previously converted files with the same flags (reused by symlink; recorded)."""
    d = os.path.join(work, "mzml")
    os.makedirs(d, exist_ok=True)
    out, info = {}, {}

    def one(raw):
        run = os.path.basename(raw).rsplit(".", 1)[0]
        mz = os.path.join(d, run + ".mzML")
        if os.path.exists(mz):
            return run, mz, {"cached": True}
        if mzml_cache and os.path.exists(os.path.join(mzml_cache, run + ".mzML")):
            os.symlink(os.path.join(mzml_cache, run + ".mzML"), mz)
            return run, mz, {"reused_from": mzml_cache}
        loc = os.path.join(d, run + ".raw")
        shutil.copyfile(raw, loc)
        try:
            rec = timed_run([TRFP, "-i", loc, "-o", d, "-f", "2", "-l", "3"], os.path.join(d, run + ".trfp.log"))
        finally:
            os.remove(loc)
        return run, mz, rec

    with ThreadPoolExecutor(threads) as ex:
        for run, mz, rec in ex.map(one, raws):
            out[run] = mz
            info[run] = rec
            log(f"mzML {run}: {rec}")
    return {"files": out, "info": info}


# ------------------------------------------------------------------------------------------------ Sage
def sage_config(fasta, mzmls, outdir, ds_cfg, predict_rt=True, chimera=False, report_psms=1) -> dict:
    """chimera/report_psms: E7b fair comparator (results/2026-10-28_E7b_fair_speed/config.yaml); defaults = the T1
    comparator config."""
    return {
        "database": {
            "bucket_size": 8192,
            "enzyme": {"missed_cleavages": 0, "min_len": 7, "max_len": 50, "cleave_at": "$", "restrict": None,
                       "c_terminal": True},
            "peptide_min_mass": 500.0, "peptide_max_mass": 5000.0,
            "ion_kinds": ["b", "y"], "min_ion_index": 2,
            "static_mods": {"C": 57.021464}, "variable_mods": {"M": [15.9949]}, "max_variable_mods": 2,
            "decoy_tag": "rev_", "generate_decoys": False, "fasta": fasta,
        },
        "precursor_tol": {"ppm": [-ds_cfg["prec_ppm"], ds_cfg["prec_ppm"]]},
        "fragment_tol": {"ppm": [-ds_cfg["frag_ppm"], ds_cfg["frag_ppm"]]},
        "precursor_charge": [2, 4], "isotope_errors": [0, 2],
        "deisotope": False, "chimera": bool(chimera), "wide_window": False, "predict_rt": bool(predict_rt),
        "min_peaks": 15, "max_peaks": 150, "min_matched_peaks": 4, "max_fragment_charge": 1,
        "report_psms": int(report_psms), "output_directory": outdir, "mzml_paths": list(mzmls),
    }


def run_sage(rep, fasta, mzmls, out_root, ds_cfg, threads, predict_rt=True, tag="sage", lock=False, chimera=False,
             report_psms=1) -> dict:
    out = os.path.join(out_root, tag, rep)
    done = os.path.join(out, "done.json")
    if os.path.exists(done):
        return json.load(open(done))
    os.makedirs(out, exist_ok=True)
    cfg = sage_config(fasta, mzmls, out, ds_cfg, predict_rt, chimera, report_psms)
    cp = os.path.join(out, "sage_config.json")
    json.dump(cfg, open(cp, "w"), indent=1)
    env = {**os.environ, "RUST_LOG": "info", "RAYON_NUM_THREADS": str(threads)}
    rec = timed_run([SAGE, cp, "--write-pin", "--batch-size", str(max(1, min(len(mzmls), threads))),
                     "--disable-telemetry-i-dont-want-to-improve-sage"], os.path.join(out, "sage.log"), env=env,
                    lock=lock)
    rec["threads"] = threads
    json.dump(rec, open(done, "w"), indent=1)
    return rec


def rt_check(sage_dir) -> dict:
    """Does predict_rt take effect?  Sage fits its RT model only when predict_rt is on; then predicted_rt and
    delta_rt_model vary and correlate with the observed (aligned) RT for confident targets."""
    import pandas as pd
    t = pd.read_csv(os.path.join(sage_dir, "results.sage.tsv"), sep="\t",
                    usecols=["label", "peptide_q", "aligned_rt", "predicted_rt", "delta_rt_model"])
    good = t[(t.label == 1) & (t.peptide_q <= 0.01)]
    from scipy.stats import spearmanr
    return {"n_confident": int(len(good)),
            "frac_predicted_rt_nonzero": float((t.predicted_rt != 0).mean()),
            "sd_predicted_rt": float(t.predicted_rt.std()),
            "median_abs_delta_rt_model_confident": float(good.delta_rt_model.abs().median()) if len(good) else None,
            "spearman_aligned_vs_predicted_confident": float(spearmanr(good.aligned_rt, good.predicted_rt)[0])
            if len(good) > 10 else None,
            "pin_has_rt_features": "sqrt(delta_rt_model)" in open(os.path.join(sage_dir, "results.sage.pin")).readline()}


# ------------------------------------------------------------------------------------------------ Comet
COMET_TEMPLATE = """# comet_version 2026.01 rev. 1 (e4f767c); written by bench/t1/comparators.py (peptide-level DB, external decoys)
database_name = {db}
decoy_search = 0
decoy_prefix = rev_
num_threads = {threads}
peptide_mass_tolerance_upper = {prec}
peptide_mass_tolerance_lower = -{prec}
peptide_mass_units = 2
precursor_tolerance_type = 1
isotope_error = 2
search_enzyme_number = 11
search_enzyme2_number = 0
sample_enzyme_number = 1
num_enzyme_termini = 2
allowed_missed_cleavage = 0
variable_mod01 = 15.9949 M 0 2 -1 0 0 0.0
variable_mod02 = 0.0 X 0 3 -1 0 0 0.0
variable_mod03 = 0.0 X 0 3 -1 0 0 0.0
variable_mod04 = 0.0 X 0 3 -1 0 0 0.0
variable_mod05 = 0.0 X 0 3 -1 0 0 0.0
max_variable_mods_in_peptide = 2
require_variable_mod = 0
fragment_bin_tol = 0.02
fragment_bin_offset = 0.0
theoretical_fragment_ions = 0
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
spectrum_batch_size = 0
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
add_K_lysine = 0.0
add_Nterm_peptide = 0.0

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


def run_comet(rep, fasta, mzml, out_root, ds_cfg, threads) -> dict:
    out = os.path.join(out_root, "comet", rep)
    run = os.path.basename(mzml)[:-5]
    done = os.path.join(out, run + ".done.json")
    if os.path.exists(done):
        return json.load(open(done))
    os.makedirs(out, exist_ok=True)
    pp = os.path.join(out, "comet.params")
    if not os.path.exists(pp):
        open(pp + f".{run}", "w").write(COMET_TEMPLATE.format(db=fasta, threads=threads, prec=ds_cfg["prec_ppm"]))
        os.replace(pp + f".{run}", pp)
    rec = timed_run([COMET, f"-P{pp}", f"-N{os.path.join(out, run)}", mzml], os.path.join(out, run + ".comet.log"),
                    cwd=out)
    rec["threads"] = threads
    json.dump(rec, open(done, "w"), indent=1)
    return rec


# ------------------------------------------------------------------------------------------------ mokapot
def ensure_mokapot(log=print):
    if not os.path.exists(os.path.join(MOKAPOT_ENV, "bin", "mokapot")):
        log("building mokapot env (scripts/baselines/setup_mokapot_env.sh)")
        subprocess.run(["bash", os.path.join(ROOT, "scripts", "baselines", "setup_mokapot_env.sh")], check=True,
                       env={**os.environ, "VULCAN_MOKAPOT_ENV": MOKAPOT_ENV})
    return os.path.join(MOKAPOT_ENV, "bin", "mokapot")


def sorted_pin(src: str, dst: str) -> str:
    """Copy of a PIN with data rows in a canonical order (the row content after SpecId, then SpecId). Sage writes PIN
    rows, and numbers its SpecIds, in thread-scheduling order; mokapot's seeded fold split depends on row order, so
    without this mokapot results differ between identical searches."""
    with open(src) as fh:
        hdr = fh.readline()
        rows = fh.readlines()
    rows.sort(key=lambda l: (l.split("\t", 1)[-1], l))
    with open(dst + ".partial", "w") as fh:
        fh.write(hdr)
        fh.writelines(r if r.endswith("\n") else r + "\n" for r in rows)
    os.replace(dst + ".partial", dst)
    return dst


def run_mokapot(pins, out) -> dict:
    done = os.path.join(out, "done.json")
    if os.path.exists(done):
        return json.load(open(done))
    os.makedirs(out, exist_ok=True)
    pins = [sorted_pin(p, os.path.join(out, f"in{i}_" + os.path.basename(p))) for i, p in enumerate(pins)]
    cmd = [os.path.join(MOKAPOT_ENV, "bin", "python"), os.path.join(os.path.dirname(__file__), "mokapot_seeded.py"),
           "--dest_dir", out, "--seed", "11", "--max_workers", "1",
           "--keep_decoys", "--decoy_prefix", "rev_"] + (["--aggregate"] if len(pins) > 1 else []) + list(pins)
    env = {**os.environ, "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    rec = timed_run(cmd, os.path.join(out, "mokapot.log"), env=env)
    json.dump(rec, open(done, "w"), indent=1)
    return rec


# ------------------------------------------------------------------------------------------------ parsing
_MODRE = re.compile(r"\[([^\]]*)\]")


def parse_peptide(p: str):
    """Tool peptide string -> (I/L-collapsed plain sequence, Met-ox position mask). Accepts Sage ('PEPM[+15.9949]K'),
    Comet/mokapot ('-.PEPM[15.9949]K.-') notation; fixed CAM brackets are ignored."""
    if p.count(".") >= 2 and re.match(r"^[A-Za-z\-n]?\.", p):
        p = p[p.index(".") + 1: p.rindex(".")]
    seq, mask, i = [], 0, 0
    while i < len(p):
        c = p[i]
        i += 1
        if c == "[":                      # n-terminal bracket before the first residue (none expected)
            j = p.index("]", i)
            i = j + 1
            continue
        if c == "-":
            continue
        if i < len(p) and p[i] == "[":
            j = p.index("]", i)
            m = float(p[i + 1:j])
            i = j + 1
            if c == "M" and abs(m - 15.9949) < 1e-2:
                mask |= 1 << len(seq)
        seq.append(c)
    return "".join(seq).replace("I", "L"), mask


def read_ragged(path, last="Proteins"):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        n = len(hdr)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d = dict(zip(hdr[:n - 1], f[:n - 1]))
            d[hdr[n - 1]] = "\t".join(f[n - 1:]) if hdr[n - 1] == last else f[n - 1]
            yield d


# ------------------------------------------------------------------------------------------------ sage+pred
class RowMap:
    """(sequence, mask) -> row of a replicate index."""

    def __init__(self, l1):
        from .dbs import seq_bytes
        sb = seq_bytes(l1.residues, l1.seq_off, l1.length)
        self.order = np.lexsort((l1.mask, sb))
        self.k_seq = sb[self.order]
        self.k_mask = l1.mask[self.order]
        self.width = sb.dtype.itemsize

    def rows(self, seqs, masks):
        q = np.array([s.encode() for s in seqs], dtype=f"S{self.width}")
        lo = np.searchsorted(self.k_seq, q, side="left")
        hi = np.searchsorted(self.k_seq, q, side="right")
        out = np.full(len(seqs), -1, np.int64)
        for i in range(len(seqs)):
            if hi[i] > lo[i]:
                k = np.flatnonzero(self.k_mask[lo[i]:hi[i]] == np.uint64(masks[i]))
                if len(k):
                    out[i] = self.order[lo[i] + k[0]]
        return out


def sage_pred_pin(scorer, sage_dir, rep_index, stores_by_run, out_dir, log=print, chimeric=False) -> dict:
    """Sage PIN + Vulcan's predicted L3 features on Sage's own top-1 PSMs -> <out_dir>/sage_pred.pin.

    chimeric (E7b, Sage chimera mode, several PSMs per spectrum): (i) Sage's `rank` is passed to the feature builder,
    so the label-blind per-run calibrations use Sage's rank-1 PSMs only (identical to the default when every rank is
    1); (ii) ScanNr is rewritten to ScanNr * 16 + rank in the PIN, because mokapot keeps one PSM per
    (FileName, ScanNr, ExpMass) and would otherwise drop every rank >= 2 PSM."""
    import pandas as pd

    from vulcan.index.build import L1Index
    pin_out = os.path.join(out_dir, "sage_pred.pin")
    info_p = os.path.join(out_dir, "pin_info.json")
    if os.path.exists(info_p):
        return json.load(open(info_p))
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.perf_counter()
    l1 = L1Index(rep_index)
    rm = RowMap(l1)
    rows = list(read_ragged(os.path.join(sage_dir, "results.sage.pin")))
    pin = pd.DataFrame(rows)
    hdr = open(os.path.join(sage_dir, "results.sage.pin")).readline().rstrip("\n").split("\t")
    tsv = pd.read_csv(os.path.join(sage_dir, "results.sage.tsv"), sep="\t",
                      usecols=["psm_id", "filename", "scannr", "charge", "hyperscore", "rt", "peptide"])
    tsv["SpecId"] = tsv.psm_id.astype(str)
    pin = pin.merge(tsv[["SpecId", "scannr", "charge", "hyperscore", "rt"]], on="SpecId", how="left")
    assert pin.scannr.notna().all(), "PIN/TSV join failed"
    parsed = [parse_peptide(p) for p in pin.Peptide]
    pf_row = rm.rows([s for s, _ in parsed], [m for _, m in parsed])
    n_unmapped = int((pf_row < 0).sum())
    feats = {k: np.zeros(len(pin), np.float32) for k in PRED_PIN_FEATURES}
    pin["run"] = pin.FileName.str.replace(r"\.mzML$", "", regex=True)
    scan = pin.scannr.str.extract(r"scan=(\d+)")[0].astype(np.int64).values
    info = {"n_psms": int(len(pin)), "n_unmapped_peptides": n_unmapped, "runs": {}}
    for run, sp in stores_by_run.items():
        m = np.flatnonzero((pin.run.values == run) & (pf_row >= 0))
        pos = {int(s): i for i, s in enumerate(sp.scan)}
        si = np.array([pos.get(int(s), -1) for s in scan[m]], np.int64)
        ok = si >= 0
        m, si = m[ok], si[ok]
        psm = {"pf_row": pf_row[m], "charge": pin.charge.values[m].astype(np.int8), "spec_idx": si,
               "rank": (pin["rank"].values[m].astype(np.int64).astype(np.int8) if chimeric
                        else np.ones(len(m), np.int8)), "hyperscore": pin.hyperscore.values[m].astype(np.float64),
               "rt": sp.rt[si].astype(np.float64)}
        f, finfo = scorer.features(l1, psm, sp)
        for k in PRED_PIN_FEATURES:
            feats[k][m] = f[k]
        info["runs"][run] = {"n_psms": int(len(m)), "n_scan_unmatched": int((~ok).sum()), "nce": finfo["nce"],
                             "rt_mad_min": finfo["rt_mad_min"]}
    for k in PRED_PIN_FEATURES:
        pin[k] = feats[k]
    if chimeric:
        rk = pin["rank"].astype(np.int64)
        assert int(rk.max()) < 16, "rank >= 16 breaks the ScanNr * 16 + rank key"
        pin["ScanNr"] = pin["ScanNr"].astype(np.int64) * 16 + rk
        info["n_rank_ge2"] = int((rk >= 2).sum())
    cols = hdr[:-2] + PRED_PIN_FEATURES + hdr[-2:]          # features before Peptide, Proteins
    with open(pin_out + ".partial", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        sub = pin[cols]
        for row in sub.itertuples(index=False):
            fh.write("\t".join(str(x) for x in row) + "\n")
    os.replace(pin_out + ".partial", pin_out)
    info["wall_s"] = time.perf_counter() - t0
    json.dump(info, open(info_p, "w"), indent=1, default=float)
    log(f"sage+pred PIN {out_dir}: {info['n_psms']} PSMs, unmapped {n_unmapped}")
    return info
