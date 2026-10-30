"""Read-back throughput of a store into numpy (E1).

usage: python scripts/ingest/readback.py <store.vspec.arrow> [--cold] [--mmap auto|yes|no] [--reps N]
--cold evicts the file from the page cache first (posix_fadvise DONTNEED; works on local disk
and on the Lustre client). Prints JSON: wall/cpu for (a) Arrow table read and (b) flattening
peaks into contiguous numpy + scalar columns, plus GB/s and spectra/s.
"""
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from vulcan.io import store as St  # noqa: E402


def evict(path):
    fd = os.open(path, os.O_RDONLY)
    try:
        os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
    finally:
        os.close(fd)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("store")
    ap.add_argument("--cold", action="store_true")
    ap.add_argument("--mmap", default="auto", choices=["auto", "yes", "no"])
    ap.add_argument("--reps", type=int, default=1)
    a = ap.parse_args()
    St.configure_threads(1)
    mm = None if a.mmap == "auto" else a.mmap == "yes"
    size = os.path.getsize(a.store)
    out = []
    for _ in range(a.reps):
        if a.cold:
            evict(a.store)
        w0, c0 = time.perf_counter(), time.process_time()
        t = St.read_table(a.store, mmap=mm)
        w1, c1 = time.perf_counter(), time.process_time()
        P = St.peaks_numpy(t, with_im=t.column("im_index").null_count < t.num_rows)
        cols = {c: St.column_numpy(t, c) for c in ("ms_level", "rt_min", "precursor_mz", "precursor_charge", "scan")}
        # touch the data so mmap'd pages are actually read
        chk = float(P.mz.sum(dtype="f8")) + float(P.intensity[: len(P.intensity) // 2].sum(dtype="f8"))
        w2, c2 = time.perf_counter(), time.process_time()
        out.append(dict(read_table_wall=w1 - w0, read_table_cpu=c1 - c0, to_numpy_wall=w2 - w1, to_numpy_cpu=c2 - c1,
                        total_wall=w2 - w0, total_cpu=c2 - c0, n_spectra=t.num_rows, n_peaks=int(len(P.mz)),
                        bytes=size, gb_per_s=size / 1e9 / (w2 - w0), spectra_per_s=t.num_rows / (w2 - w0),
                        mmap=(mm if mm is not None else os.path.abspath(a.store).startswith("/tmp/")), cold=a.cold,
                        checksum=chk, n_cols=len(cols)))
    print(json.dumps(out[-1] if a.reps == 1 else out))


if __name__ == "__main__":
    main()
