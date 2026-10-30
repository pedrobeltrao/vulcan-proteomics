# E1: ingest throughput (DESIGN §10)

Written by the PI from metrics.json (`summary` block), because the agent could not write this file.

**Hardware:** 4 pinned cores. Timings were taken under the bench lock; core-s = user + sys.

| Path | Core-s per 10⁵ MS2 |
|---|---|
| **Thermo direct** (in-process RawFileReader, .NET 8 helper `vulcan/io/dotnet/VulcanRaw`) | **17.8** |
| TRFP → mzML → parser | 531 |
| TRFP → Parquet | 434 (incomplete metadata) |
| Bruker native Rust reader (`vulcan/io/rust/vulcan_tims`) | 26.9 (36.6 with unzip) |
| Bruker timsrust_pyo3 | 96.2 |

**Decision rule:** the best path is below 30 core-s per 10⁵ MS2, so the headline is **not ingest-bound** on Thermo QE data. Bruker sits near the threshold.

**Correctness:**
- Thermo: the direct and TRFP-mzML stores are bit-identical on all peaks and fields, and 200/200 spectra match pyteomics.
- Bruker: the native and timsrust readers are bit-identical.
- Known differences:
  - precursor_mz is the instrument monoisotopic value, not TRFP's selected ion. Both are stored.
  - 18.6% of diaPASEF events fall outside isolation windows and are not stored.

**Operational rules:**
- Sequentially prefetch RAW files; never let RawFileReader mmap cold Lustre (2–10× wall-time penalty).
- Keep Thermo stores uncompressed; use lz4 for Bruker.
- Run about 3 long-lived workers.

**Gaps:** Exploris, Astral, Lumos-IT, FAIMS, SPS, EThcD and ddaPASEF are implemented but not yet measured.
