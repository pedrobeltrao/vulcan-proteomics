"""Reference parse speed: iterate the first N spectra of an mzML with pyteomics (E1 baseline).
usage: python scripts/ingest/pyteomics_iter.py <file.mzML> [N=10000]; prints JSON."""
import json
import sys
import time

from pyteomics import mzml

path = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 10000
t0, c0 = time.perf_counter(), time.process_time()
n = npk = 0
with mzml.MzML(path, decode_binary=True, use_index=False) as r:
    for s in r:
        npk += len(s["m/z array"])
        n += 1
        if n >= N:
            break
w, c = time.perf_counter() - t0, time.process_time() - c0
print(json.dumps(dict(n_spectra=n, n_peaks=npk, wall_s=w, cpu_s=c, spectra_per_s=n / w)))
