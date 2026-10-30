# Benchmark conventions

- **Exclusive timing.** Any run whose timing will be reported must hold the node-wide benchmark lock:
  `flock ${VULCAN_BENCH_LOCK:-/tmp/vulcan_bench.lock} <command>` (or `scripts/baselines/timed.py --lock`).
  Development and correctness runs don't need the lock. Agents share the node, so unlocked timings are noisy and must not be reported.
- **Thread pinning.** Size every pool from `vulcan.hw.n_threads()` (affinity mask, or `$VULCAN_THREADS`) and apply it with `vulcan.hw.pin_thread_pools(n)` / `hw.thread_env(n)`: `OMP_NUM_THREADS`, `NUMBA_NUM_THREADS`, `MKL_NUM_THREADS`, `RAYON_NUM_THREADS`, Comet `num_threads`, DIA-NN `--threads N`, `torch.set_num_threads`. Never use `os.cpu_count()`. For M1-continuity comparisons run with `VULCAN_THREADS=4` (or `timed.py --threads 4`); report full-node and 4-thread runs separately.
- **Hardware fingerprint.** Every result record carries `vulcan.hw.fingerprint()` (CPU model, affinity count, AVX-512 flags, GPU, driver, RAM, hostname); `python bench/hwinfo.py out.json` writes it standalone. `timed.py` adds it as `hw` to each timings.jsonl line. Speed comparisons are only valid within one node and session.
- **Record keeping.** Record the wall time, CPU core-seconds (`/usr/bin/time -v` user+sys), GPU-seconds, peak RSS and peak GPU memory.
- **Output layout.** Results go in `results/<date>_<experiment>/` as `metrics.json`, `config.*`, `env.txt` and `git_sha`.
- **Sealed data.** Files in `data/sealed/` are mode 000 and may only be read through `python3 bench/sealed_run.py --dataset T1 --purpose "..." -- <command, {files} expands to the sealed paths>`.
  It refuses unless the git tree is clean and HEAD (ignoring its own ledger commits) has a `milestone-*` tag, appends a line to `results/SEALED_LEDGER.md` and commits it, unlocks the files, runs the command, and re-locks.
  Use `--dry-run` to check preconditions. Never open sealed raw files any other way.
