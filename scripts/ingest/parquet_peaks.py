"""TRFP Parquet (long format, one row per peak) -> per-spectrum offsets + numpy peaks (E1).
TRFP 2.0 mzparquet lacks activation/NCE/analyzer/filter/injection time/FAIMS, so this path can
not build a complete store; this measures only its parse cost. Prints JSON."""
import json
import sys
import time

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

pa.set_cpu_count(1)
pa.set_io_thread_count(1)
t0, c0 = time.perf_counter(), time.process_time()
t = pq.read_table(sys.argv[1], columns=["scan", "level", "rt", "mz", "intensity", "precursor_mz", "precursor_charge"],
                  use_threads=False)
scan = t.column("scan").to_numpy()
starts = np.flatnonzero(np.r_[True, scan[1:] != scan[:-1]])
offsets = np.r_[starts, len(scan)].astype(np.int64)
mz = t.column("mz").to_numpy()
it = t.column("intensity").to_numpy()
lvl = t.column("level").to_numpy()[starts]
w, c = time.perf_counter() - t0, time.process_time() - c0
print(json.dumps(dict(n_spectra=int(len(starts)), n_ms2=int((lvl == 2).sum()), n_peaks=int(len(mz)), wall_s=w, cpu_s=c)))
