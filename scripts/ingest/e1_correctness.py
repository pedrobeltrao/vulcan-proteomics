"""E1 correctness checks on the stores produced by e1_bench.py (not timed; no lock needed).

- direct vs TRFP-mzML-path stores of the same RAW: counts, every peak, scalar metadata
- pyteomics re-read of sampled spectra from the TRFP mzML vs both stores
- MS1/MS2 counts vs TRFP's own mzML spectrumList count and per-level counts
- Bruker: native vs timsrust_pyo3 stores, window/IM invariants
Writes <out>/correctness.json.
"""
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from vulcan.io import check  # noqa: E402
from vulcan.io import store as St  # noqa: E402

STORE = "/tmp/vulcan_store/e1"
STAGE = "/tmp/vulcan_stage/e1"


def mzml_counts(path):
    """Spectrum count and MS-level histogram straight from the mzML text (independent of both parsers)."""
    head = open(path, "rb").read(1 << 16)
    m = re.search(rb'<spectrumList count="(\d+)"', head)
    lv = {}
    with open(path, "rb") as fh:
        for line in fh:
            if b'accession="MS:1000511"' in line:
                v = re.search(rb'value="(\d+)"', line).group(1).decode()
                lv[v] = lv.get(v, 0) + 1
    return int(m.group(1)) if m else None, lv


def main():
    out = sys.argv[1]
    R = {}
    for k in ("01A", "01B"):
        a, b = f"{STORE}/{k}_direct.vspec.arrow", f"{STORE}/{k}_mzml.vspec.arrow"
        if os.path.exists(a) and os.path.exists(b):
            R[f"{k}_direct_vs_mzml_path"] = check.compare_stores(a, b)
            R[f"{k}_counts"] = {"direct": {x: St.run_meta(a)[x] for x in ("n_spectra", "n_ms1", "n_ms2", "n_peaks")},
                                "mzml_path": {x: St.run_meta(b)[x] for x in ("n_spectra", "n_ms1", "n_ms2", "n_peaks")}}
    mz = f"{STAGE}/b1906_293T_proteinID_01A_QE3_122212.mzML"
    if os.path.exists(mz):
        n, lv = mzml_counts(mz)
        R["01A_trfp_mzml_text_counts"] = {"spectrumList_count": n, "ms_level_hist": lv}
        for k, p in (("direct", f"{STORE}/01A_direct.vspec.arrow"), ("mzml_path", f"{STORE}/01A_mzml.vspec.arrow")):
            R[f"01A_pyteomics_check_{k}"] = check.check_against_pyteomics(p, mz, n=200, seed=11)
    # precursor semantic summary (direct): fraction with mono, charge distribution
    p = f"{STORE}/01A_direct.vspec.arrow"
    if os.path.exists(p):
        t = St.read_table(p, columns=["ms_level", "mono_mz", "precursor_charge", "master_scan", "activation", "nce",
                                      "injection_time_ms", "analyzer", "filter_string"])
        lvl = t.column("ms_level").to_numpy()
        ms2 = lvl == 2
        ch = St.column_numpy(t, "precursor_charge")[ms2]
        R["01A_direct_ms2_metadata"] = dict(
            frac_mono_present=float(np.mean(~np.isnan(St.column_numpy(t, "mono_mz")[ms2]))),
            charge_hist={str(int(c)) if not np.isnan(c) else "null": int(v) for c, v in zip(*np.unique(ch, return_counts=True))},
            frac_master_scan_present=float(np.mean(~np.isnan(St.column_numpy(t, "master_scan")[ms2]))),
            nce_values=sorted(set(np.round(St.column_numpy(t, "nce")[ms2], 2).tolist())),
            activation_codes=sorted(set(St.column_numpy(t, "activation")[ms2].tolist())),
            analyzer_codes=sorted(set(t.column("analyzer").to_numpy().tolist())),
            n_filter_strings_null=int(t.column("filter_string").null_count))
    a, b = f"{STORE}/A9_R1_frames.vspec.arrow", f"{STORE}/A9_R1_frames_pyo3.vspec.arrow"
    if os.path.exists(a) and os.path.exists(b):
        ta, tb = St.read_table(a), St.read_table(b)
        pa_, pb = St.peaks_numpy(ta, True), St.peaks_numpy(tb, True)
        same = dict(n_rows=ta.num_rows == tb.num_rows, offsets=bool(np.array_equal(pa_.offsets, pb.offsets)),
                    mz=bool(np.array_equal(pa_.mz, pb.mz)), intensity=bool(np.array_equal(pa_.intensity, pb.intensity)),
                    im_index=bool(np.array_equal(pa_.im_index, pb.im_index)))
        for c in ("scan", "ms_level", "isolation_target_mz", "window_index", "im_scan_begin", "im_scan_end", "nce", "rt_min"):
            same[c] = bool(ta.column(c).equals(tb.column(c)))
        ma, mb = St.run_meta(a), St.run_meta(b)
        lvl = ta.column("ms_level").to_numpy()
        sb, se = St.column_numpy(ta, "im_scan_begin"), St.column_numpy(ta, "im_scan_end")
        rows = np.nonzero(lvl == 2)[0]
        seg_ok = all(
            (lambda s: len(s) == 0 or (s.min() >= sb[i] and s.max() < se[i]))(pa_.im_index[pa_.offsets[i]:pa_.offsets[i + 1]])
            for i in rows[:: max(1, len(rows) // 2000)])
        R["A9_R1_bruker_native_vs_pyo3"] = dict(identical=same, all_identical=all(same.values()),
                                                 n_rows=ta.num_rows, n_ms1=ma["n_ms1"], n_ms2_rows=ma["n_ms2"],
                                                 n_events_stored=ma["n_peaks"],
                                                 events_outside_windows=[ma.get("events_outside_windows"), mb.get("events_outside_windows")],
                                                 im_index_within_window_rows=bool(seg_ok))
    with open(os.path.join(out, "correctness.json"), "w") as fh:
        json.dump(R, fh, indent=1, default=str)
    print(json.dumps({k: (v.get("ok") if isinstance(v, dict) and "ok" in v else v.get("n_bad") if isinstance(v, dict) and "n_bad" in v else "") for k, v in R.items()}))


if __name__ == "__main__":
    main()
