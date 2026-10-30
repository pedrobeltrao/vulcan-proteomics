# Vulcan: GPU/CPU co-designed engine for repository-scale MS proteomics

Vulcan reprocesses data-dependent acquisition (DDA) runs from Thermo RAW files to scored peptides and protein groups. The goal is to run fast enough for whole-repository reanalysis (PRIDE) while keeping error rates under statistical control.

> [!WARNING]
> **Disclaimer: this is the output of an AI-agent experiment.** An AI agent team (Claude Opus 5.5 / Sonnet 5.5) produced the code, analyses, reports and manuscript in this repository. It worked from the prompt in [PROMPT.md](PROMPT.md). **Human experts have not reviewed or verified any of it.** It may contain errors in code, statistics or interpretation.
>
> **There is no held-out test-set comparison.** Every result below comes from development datasets, several of which were used for tuning. The single pre-registered sealed evaluation (T5 + T1) was started. The compute node was lost while T1 was running, so the evaluation was closed as "failed (unterminated)" and its outputs were deleted unread. T1–T5 remain unseen. No false-discovery-rate (FDR) claim has passed the pre-registered null-spectrum gates. See [REPORT.md](REPORT.md) §0 and §4.

## Documents

| Document | Description | Files |
|---|---|---|
| Prompt | Original task given to the agent team (verbatim) | [PROMPT.md](PROMPT.md) |
| Final report | Approach, best methods, results, gate history, limitations | [REPORT.md](REPORT.md) |
| Research tree | Every branch with its status and decision notes D0–D40 | [RESEARCH_TREE.md](RESEARCH_TREE.md) |
| Briefing | Short summary for non-specialists | [PDF](reports/briefing/BRIEFING.pdf) · [Markdown](reports/briefing/BRIEFING.md) |
| Technical report | Detailed expert report, including an appendix of inconsistencies | [PDF](reports/technical/TECHNICAL_REPORT.pdf) · [Markdown](reports/technical/TECHNICAL_REPORT.md) |
| Manuscript draft | Journal-style write-up (dev data only) | [PDF](reports/manuscript/MANUSCRIPT.pdf) · [Markdown](reports/manuscript/MANUSCRIPT.md) |
| Design / literature / data | Pre-registered design, literature review, dataset provenance and licences | [DESIGN.md](DESIGN.md) · [LITERATURE.md](LITERATURE.md) · [DATA.md](DATA.md) |
| Process lessons | Rules the agents kept for running the project | [LESSONS.md](LESSONS.md) |

## Key results (development data only; no held-out test-set comparison)

The table is copied from [REPORT.md](REPORT.md) §0. Hardware: a shared node with 32 cores of a 2× AMD EPYC 9354 (Zen 4) and one NVIDIA A100 80 GB PCIe. "Core-s" means CPU core-seconds; it excludes GPU time, which Vulcan uses and the CPU-only comparators do not.

| Workflow | Accuracy vs comparator (dev) | Speed vs comparator (dev, 32-core EPYC + A100, shared node) |
|---|---|---|
| Standard DDA search (frozen v1) | 0.990× [0.989, 0.992] peptidoforms vs Sage + predicted features + mokapot at matched 1% FDP (D0); +2.8% vs plain Sage | search stage about 23–40× fewer CPU core-s than Sage at 4–32 threads (measured single-process 6.0 core-s) (CPU only; excludes GPU time); RAW → proteins about 5,700 runs/h at about 9 core-s/run vs Sage 70–120 core-s (4–32 threads, from mzML) + about 250 for RAW→mzML |
| Chimeric search (v2_sep, candidate) | 1.150× [1.069, 1.203] peptidoforms vs Sage chimera + predicted features + mokapot on D8h (10 dev runs held out from scorer selection; not a sealed test set); protein groups 0.999× | about 3,900 runs/h at about 15 core-s/run |
| Ion-trap CID / TMT | 1.034× (D4) and 1.093× (D10cf) vs Sage + Prosit features + mokapot | about 0.9 core-s per fraction vs Sage about 54 core-s per file (16 threads, unlocked, indicative) |
| Label-free protein abundance | parity with Sage on the HYE benchmark (matched-set error 0.118 vs 0.1145) | about 1,200–1,280 runs/h at about 15 core-s/run vs Sage search+LFQ 266 core-s/run (+ about 250 for mzML conversion) |
| Phospho + site localization | Ala-decoy FLR point estimates near true; stricter upper-bound rule not met on PXD009449 (n = 183); mixed vs pyAscore; not run-to-run reproducible on about 0.7% of PSMs | about 3,200 runs/h; Sage 3,087 core-s for 2 runs, without localization |
| Open search | −4% peptide stems vs Sage-open | search stage about 65× fewer core-s than Sage-open (32 threads) |
| Cross-linking | DSSO (D3b): about 73% of xiSEARCH heteromeric pairs, at lower true FDP | DSS (D3): about 460× less wall time than xiSEARCH for the whole pipeline, about 8,900× for the search stage |
| Repository-scale FDR (L6) | global peptide FDP̂ 1.16% (arm P) / 1.21% (arm C) and protein 1.00% at N = 1,000 (982 runs completed), vs 9.3% (peptide) and 19.1% (protein) for the naive union | streaming from PRIDE is download-bound, about 580 runs/h |

FDP = false-discovery proportion. FDP̂ is estimated by entrapment; it is not a verified error rate. Arms P and C are the shuffled and natural-foreign entrapment databases. The frozen v1 scorer **failed** its pre-registered null-spectrum gates (tail-G2, G1, G3). Later amendments were either not evaluable or not run ([REPORT.md](REPORT.md) §1 and §4).

## Installation

### Requirements

- Linux x86-64 and an NVIDIA GPU with compute capability 8.0 or 8.6 (for example A100, A40 or RTX 30xx). The CUDA kernels are built for sm_80 + sm_86 + PTX, and you need a driver that supports CUDA 12.x. Development used an A100 80 GB. Large indexes, such as the phospho index (34 GB), need a GPU with a lot of memory.
- conda or mamba (conda-forge).
- A Rust toolchain ([rustup](https://rustup.rs), stable) to build the PyO3 crates.
- For Thermo RAW input: the .NET 8 runtime/SDK and Thermo's RawFileReader assemblies, which you obtain separately (see below).

### 1. Clone the repository and create the environment

```bash
git clone https://github.com/pedrobeltrao/vulcan-proteomics
cd vulcan-proteomics
mamba env create -f environment.yml      # exact recorded versions: env/environment.lock.yml
conda activate vulcan
. env/activate.sh                         # sets VULCAN_ROOT, VULCAN_SCRATCH, CUDA_HOME, CUDA arch flags
```

The code runs from the repository root; `pyproject.toml` installs only the pure-Python package (`pip install -e .` is optional). Compiled artefacts go to `tools/`, which git ignores.

### 2. Build the native components

```bash
bash scripts/m2/build_vulcan_core.sh      # PyO3 crate on sage-core (Sage v0.14.7): digestion, mods, masses
                                          #   -> tools/vulcan_core/vulcan_core.abi3.so  (cargo build --release)
bash scripts/m2/build_kernels.sh          # CUDA kernels (vulcan/gpu/kernels.cu) -> tools/vulcan_kernels/
                                          #   (also rebuilt automatically on first use if the source changed)
bash scripts/ingest/build_vulcan_tims.sh  # optional: Bruker timsTOF reader (timsrust) -> tools/vulcan_tims/
```

`maturin` is part of the environment, but the scripts call `cargo` directly and copy the resulting cdylib. Set `VULCAN_CORE_DIR` to use a crate built elsewhere.

### 3. Thermo RAW support (obtain RawFileReader separately)

Vulcan reads Thermo RAW files in-process through a small C# helper (`vulcan/io/dotnet/VulcanRaw`) loaded with pythonnet. The helper links against Thermo's **RawFileReader** assemblies (`ThermoFisher.CommonCore.Data.dll`, `ThermoFisher.CommonCore.RawFileReader.dll`). **These assemblies are not part of this repository and must not be redistributed.** They come under Thermo Fisher Scientific's own licence. Obtain them yourself and accept that licence, for example by installing [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser), which ships them:

```bash
mamba create -p conda_envs/trfp -c conda-forge -c bioconda thermorawfileparser
export VULCAN_TRFP_BIN=$PWD/conda_envs/trfp/bin   # directory that contains ThermoFisher.CommonCore.*.dll
bash scripts/ingest/build_vulcanraw.sh            # installs the .NET 8 SDK into tools/dotnet-sdk if missing,
                                                  # builds tools/vulcanraw/VulcanRaw.dll
```

Optional overrides: `VULCAN_RAW_DLL` and `VULCAN_DOTNET_ROOT`. mzML input does not need any of this: `python -m vulcan.io.convert run.mzML ...` uses a streaming mzML parser.

### 4. Check the installation

```bash
python -c "import vulcan, vulcan.hw; print(vulcan.hw.fingerprint())"
pytest -q --import-mode=importlib tests/stats tests/a2 tests/a3 tests/a6 tests/report   # CPU-only, no data needed
pytest -q --import-mode=importlib tests/gpu tests/ms1                                   # needs GPU + built kernels + vulcan_core
```

## Quickstart

All commands run from the repository root after `. env/activate.sh`. Comparator tools (Sage, Comet, DIA-NN, mokapot, xiSEARCH) are only needed for the benchmark scripts, not for Vulcan itself.

**1. Build a proteome index** (trypsin, 2 missed cleavages, length 7–50, carbamidomethyl-C fixed, Met-ox variable, paired pseudo-reverse decoys; content-hashed, built once and reused across runs):

```bash
python -m vulcan.index.build --fasta human.fasta contaminants.fasta --out work/index_human
#   optional: --entrap-seeds 2001 2002 2003  (shuffled entrapment families for empirical FDP estimation)
```

**2a. One file, step by step** (RAW, mzML or Bruker .d → Arrow spectrum store → GPU search):

```bash
python -m vulcan.io.convert run.raw  work/stores/      # or run.mzML; writes work/stores/<run>.vspec.arrow
python -m vulcan.gpu.run --index work/index_human --stores work/stores/run.vspec.arrow --out work/psms
#   -> work/psms/<run>.psms.parquet (label-free PSM features, schema in vulcan/gpu/PSM_SCHEMA.md) + stats.json
```

**2b. A cohort of Thermo RAW files** in one long-lived process. Decode workers stage each RAW to local disk while the GPU index stays resident:

```bash
ls /data/project/*.raw > raws.txt
python -m vulcan.cohort --index work/index_human --raw-list raws.txt --out work/cohort \
       --stage $VULCAN_SCRATCH/vulcan_stage --workers 24 --skip-failed
#   -> work/cohort/<run>.psms.parquet and work/cohort/cohort_stats.json (per-run timings, core-s, GPU use)
```

Useful flags (see `python -m vulcan.cohort --help`): `--ms1` (MS1 features and precursor hypotheses), `--ms1-backend gpu`, `--family-top2`, `--search-threads 2`, `--no-psms` (throughput only).

**3. Production scoring and FDR** (frozen v1 `static_gbm_pred`, paired target-decoy competition and picked protein groups) runs inside the cohort process:

```bash
python -m vulcan.cohort --index work/index_human --raw-list raws.txt --out work/cohort --score-lib PREDLIB_DIR
#   adds a `score` column and writes <run>.peptides.parquet and <run>.proteins.parquet (q-values)
python -m vulcan.cohort ... --chimeric --score-lib PREDLIB_DIR   # chimeric search + v2_sep scorer (candidate)
```

`PREDLIB_DIR` is a consolidated AlphaPeptDeep predicted library (fragment intensities and RT) for the index peptides. The only builder in this repository is the research script `scripts/scoring_speed/build_library.py` (with `vulcan/predict/`). That script assumes the development index and cache layout (`scripts/m3/m3_common.py`) and needs `peptdeep==1.5.1`, so a packaged one-command path does not exist yet. The default model files are under `results/2026-10-27_M3_scoring/models/` and `results/2026-10-28_E7_chimeric_v2/models/`. Their sha256 values are pinned in `config/a2.yaml` / `config/a3.yaml`.

Other workflows (PTM and open search, cross-linking, LiP, CF-MS, LFQ/TMT, repository-scale L6 FDR) live in `vulcan/{ptm,xl,lip,cf,quant,stats}` and are driven by the milestone scripts in `scripts/` (`m4` PTM, `m5`/`m5lc` XL/LiP/CF, `m6`/`m6d` quantification, `e9` repository scale). Each `results/<date>_<milestone>/NOTES.md` gives the exact commands used.

## Reproducing figures and reports

```bash
python reports/figures/make_figures.py          # fig1-fig6 (PDF + PNG) and reports/figures/SOURCES.md from committed results/
python reports/render/render_pdfs.py            # BRIEFING / TECHNICAL_REPORT / MANUSCRIPT PDFs (needs pandoc + tectonic)
```

Each plotted number is traced to a file and a JSON key in `reports/figures/SOURCES.md`. The figure script reads only committed development results, not raw data.

## Datasets

All datasets are public. Sources, versions, sha256 checksums and licences are listed in [DATA.md](DATA.md), and per-file manifests are in `data/manifests/` (`*.tsv` with URLs and sizes, `*.sha256`). No raw data is included. The fetch scripts re-download and verify it:

```bash
python scripts/data/fetch_tier_s.py   # Tier S benchmark suite -> data/raw/<accession>/
python scripts/data/fetch_dev.py      # development sets
python scripts/data/fetch_fasta.py    # UniProt proteomes / contaminants -> data/fasta/
python scripts/data/verify_data.py    # every manifest row vs disk (presence, size, sha256)
# set VULCAN_CONTACT_EMAIL to identify yourself in the HTTP User-Agent sent to PRIDE/UniProt
```

Main PRIDE accessions: PXD001468 (D0, HEK293), PXD001077 (D1, *P. furiosus* entrapment), PXD000138 / PXD009449 / PXD007058 (synthetic phosphopeptide pools), PXD014337 and PXD029252 (cross-linking libraries), PXD016766 (TMT), PXD031616 (LiP-MS), PXD059608 and PXD015406 (co-fractionation), PXD003881 and PXD014415 (LFQ / match-between-runs truth), plus ProteoBench DDA/DIA modules and a 1,000-run human DDA stream from 142 PRIDE projects (`data/manifests/e9_d8.tsv`). The sealed test sets (T1–T5) are listed by checksum only (`data/manifests/dev_sealed.tsv`, `sealed.sha256`).

## Repository layout

```
vulcan/            engine: io (RAW/mzML/Bruker -> Arrow store), index (L1, sage-core crate in rust/), gpu (CUDA kernels,
                   search), ms1, rescore (features, scorers), predict (AlphaPeptDeep/Koina), stats (TDC, picked protein,
                   FLR, entrapment, global FDR), ptm, xl, lip, cf, quant, chimeric, encoders, gates (pre-registered tests),
                   cohort.py (batch driver), hw.py (hardware fingerprint, thread pools)
bench/             benchmark and sealed-evaluation drivers (T1/T2/T5, comparators, one-look guard)
scripts/           milestone scripts that produced every result (data fetch, baselines, m2-m7, e7-e12, a2-a6 amendments)
tests/             unit and integration tests (CPU-only and GPU)
config/            frozen pre-registration configs (A2-A6, T2)
results/           per-milestone NOTES.md, metrics.json, small tables and logs (dev data only); SEALED_LEDGER.md
reports/           briefing, technical report, manuscript (MD + PDF), figures + make_figures.py, render scripts
reviews/           independent design/amendment reviews and check scripts
data/manifests/    dataset manifests and checksums (no raw data)
lit/               literature notes behind LITERATURE.md
env/               generic activate.sh, pinned environment lock file, ENVIRONMENT.md (development-node history)
```

## Notes on this export

- This is a cleaned public snapshot with fresh git history. Commit hashes, tags (such as `prereg-v0.4.4` and `milestone-M7-2026-10-28`) and `git_sha` files in the research record refer to the private development history, which is not included. The sealed-look machinery (`bench/sealed_run.py`) checks those tags, so it cannot run as shipped.
- Cluster-specific absolute paths in the research record were replaced with placeholders: `$VULCAN_ROOT` (repository root), `$VULCAN_SCRATCH` (node-local scratch), `$VULCAN_LOOK_ROOT` (sealed-look output root) and `gpu-node` (the development host). `config/a5.yaml` and `config/t2.yaml` contain such placeholders; substitute them (for example with `envsubst`) before use. Scientific content and numbers are unchanged.
- The `results/2026-10-28_E12_speed_abundance_ptm/` result files come from an unmerged development branch. The reports and the figure script cite them, but the corresponding code changes are not part of this snapshot.
- The A5 freeze hashes in `config/a5.yaml` refer to the frozen bytes of the development repository, and the pinned files are unchanged in this export.

## Limitations

- **No held-out test-set comparison.** All numbers are from development data, and some of those sets were used for tuning.
- **Error control is not gate-validated.** The frozen v1 scorer failed the pre-registered null-spectrum gates. Later re-tests were not evaluable (A2), not adopted (A4) or infeasible on the available node (A6). FDR figures are entrapment estimates.
- Speed was measured on a shared node (load 14–110). Comparators were often run without the benchmark lock, and Vulcan's core-s figures exclude GPU time.
- Production scoring depends on a predicted-spectrum library whose build is tied to the development layout (see Quickstart step 3). Phospho site localization is not run-to-run reproducible on about 0.7% of PSMs. The DIA arm was not done (there is no DIA result), and MSFragger and DIA-NN were not run as comparators.
- Limitations for each workflow are listed in [REPORT.md](REPORT.md) and in the technical report's inconsistencies appendix.

## Licence

The code and documentation in this repository are released under the [MIT License](LICENSE) (Copyright (c) 2026 Vulcan AI agent team).

Third-party components and their licences:

| Component | Use | Licence |
|---|---|---|
| [Sage](https://github.com/lazear/sage) (sage-core, v0.14.7 / 99407db6) | Cargo git dependency of `vulcan_core` (digestion, FASTA, modifications, ion series); comparator | MIT |
| [timsrust](https://github.com/MannLabs/timsrust), timsrust_pyo3 | Bruker reader (`vulcan_tims`, optional) | Apache-2.0 |
| PyO3, rust-numpy, rayon, twox-hash | Rust crates | MIT / Apache-2.0, BSD-2-Clause, MIT |
| Thermo RawFileReader | Thermo RAW reading (**not redistributed**; obtain separately) | Thermo Fisher Scientific licence |
| [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser) | Source of the RawFileReader assemblies; RAW→mzML comparator path | Apache-2.0 (bundled RawFileReader: Thermo licence) |
| pythonnet | .NET interop | MIT |
| AlphaPeptDeep (peptdeep), alphabase | Predicted spectra/RT features and pretrained weights (downloaded at runtime) | Apache-2.0 |
| koinapy / Koina (Prosit models) | IT-CID/TMT predicted features (optional) | Apache-2.0 (models: see Koina) |
| mokapot, pyteomics, psm-utils | Comparator rescoring, mzML parsing, PSM I/O | Apache-2.0 |
| PyTorch, NumPy, SciPy, pandas, polars, pyarrow, numba, scikit-learn, matplotlib | Python stack | BSD / Apache-2.0 / MIT-style |
| DIA-NN, Comet, xiSEARCH/xiFDR, MSstatsLiP, FlashLFQ, directLFQ | Benchmark comparators only (not included) | their own licences (DIA-NN Academia: free for academic use only) |

Public datasets keep their original licences (mostly CC0 or the EBI terms of use; see [DATA.md](DATA.md)).
