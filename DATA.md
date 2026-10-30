# DATA.md - datasets, FASTAs, ground truth

Retrieved 2026-10-26, completed and verified 2026-10-27 on the 32-core node (UniProt release 2026_03). Completeness is checked by `python3 scripts/data/verify_data.py` (latest run: `results/2026-10-27_data_verify/`). All paths are relative to the project root. Raw data are **not** in git
(`data/raw/`, `data/sealed/`, `data/fasta/`, `data/cache/`, `*.raw`, `*.d.zip` are gitignored); everything here is re-creatable by script:

```
python3 scripts/data/verify_data.py        # every manifest row vs disk: presence, size, sha256 (cache keyed by path+size+mtime); sealed files by stat only
python3 scripts/data/check_pride_sha1.py  # optional: readable PRIDE files vs the SHA-1 published by PRIDE (old projects publish none)
python3 scripts/data/fetch_tier_s.py      # manifest data/manifests/tier_s.tsv -> data/raw/<acc>/, data/truth/<acc>/ (Tier S = D0,D1,D2(PXD007058),D3,D5,D6 + T3 truth)
python3 scripts/data/fetch_dev.py         # manifest data/manifests/dev_sealed.tsv, role=dev    -> data/raw/<acc>/, data/truth/<acc>/   (D2 pools, D3b, D4, D7)
python3 scripts/data/fetch_sealed.py      # manifest data/manifests/dev_sealed.tsv, role=sealed -> data/sealed/<acc>/ (chmod 000), data/truth/sealed/<acc>/   (T1, T2, T3, T5)
python3 scripts/data/build_dev_sealed_manifest.py   # regenerates dev_sealed.tsv (seeded pool selection lives in this file)
python3 scripts/data/make_pool_truth.py   # per-pool FASTAs for the PXD000138 pools (after fetch)
python3 scripts/data/fetch_fasta.py       # data/fasta/*.fasta, provenance in data/manifests/fasta.tsv
(cd data && sha256sum -c manifests/tier_s.sha256 && sha256sum -c manifests/dev_sealed.sha256)   # verify (paths relative to data/)
```

All fetch scripts take `--jobs N` (parallel downloads), `--only <ids>`, `--check`; they resume `.part` files with HTTP Range, discard any file whose size or PRIDE SHA-1 is wrong, and are idempotent (a second run only skips). Per-file sha256: `data/manifests/tier_s.sha256`, `dev_sealed.sha256` (dev raw + all truth), `derived_truth.sha256` (per-pool FASTAs written by `make_pool_truth.py`), `sealed.sha256` (sealed raw only; the files
are mode 000, so `sha256sum -c` cannot read them; `python3 bench/sealed_run.py --verify` is the sanctioned check) and `fasta.sha256`.
Per-file URL, size, role, set id and purpose: `data/manifests/tier_s.tsv`, `dev_sealed.tsv` (columns `id, accession, filename, url,
expected_size, purpose, role, set_id, member`; `member` = single file extracted from a remote PRIDE SEARCH zip by HTTP Range). URLs were resolved
from the PRIDE archive API v3 `/projects/<ACC>/files` on 2026-10-26; ProteoBench files come from `proteobench.cubimed.rub.de`. The fetcher resumes
partial downloads, verifies sizes, is idempotent, and uses at most one connection per process (run one process per set with `--only D4` etc. to parallelise).
PRIDE FTP gave 1-4 MB/s per connection; ProteoBench ~65 MB/s.

Licences (from the PRIDE project metadata unless stated): CC0: PXD007058, PXD014337, PXD029252, PXD059608, PXD016766, PXD031616, PXD052597.
"EBI terms of use" (free to download and use; cite the original papers): PXD001077, PXD001468, PXD000138, PXD009449. ProteoBench files: open access as
published by ProteoBench (DDA-QExactive module; DIA-plasma data are also on jPOST JPST003358, Distler et al., Nat Commun 2025). UniProt: CC BY 4.0.

## Roles (DESIGN.md section 6.6 is authoritative; the manifest `role` / `set_id` columns must match it)

| role | meaning | where | git |
|---|---|---|---|
| `dev` | may gate merges; numbers are development curves, never claims | `data/raw/<acc>/` (+ `data/truth/<acc>/`) | raw ignored; small truth committable |
| `sealed` | claims-grade; read only through `bench/sealed_run.py` (clean tree + `milestone-*` tag + ledger line in `results/SEALED_LEDGER.md`); sealed entrapment seeds 2004-2006 | `data/sealed/<acc>/` (mode 000); truth stays readable in `data/truth/sealed/<acc>/` | ignored |
| `train` | training sets: listed in the training manifests, not in `tier_s.tsv` / `dev_sealed.tsv` | n/a | n/a |

Sealed rules: never open or inspect sealed raw contents beyond size and checksum (the fetcher hashes each file once, right after download, then
`chmod 000`); the fetch scripts write sealed sets straight into `data/sealed/`. T4 (stratified PRIDE sample) is selected from `data/manifests/t4_frame.tsv`
and streamed through `/tmp`, so it has no entry here.

### Dev / sealed assignment (matches DESIGN.md 6.6)

| Set | role | Data in this repo | Manifest |
|---|---|---|---|
| D0 | dev | PXD001468, 2 HEK files (S0) | tier_s |
| D1 | dev | PXD001077 (S1) | tier_s |
| D2 | dev | PXD000138 dev pools + PXD009449 dev pools + **PXD007058 (S2)**; the PXD007058 truth table is **unavailable** (see S2 row), only submitter reference exports exist | tier_s + dev_sealed |
| D3 | dev | PXD014337 (S3) | tier_s |
| D3b | dev | PXD029252 `20210409_QExHFX2_RSLC3_DSSO_mainlib_rep2.raw` (+ xiSEARCH mzid) | dev_sealed |
| D4 | dev | PXD016766, 12 files | dev_sealed |
| D5 | dev | ProteoBench DIA-plasma PYE (S5) | tier_s |
| D6 | dev | PXD059608 CF fractions (S6) | tier_s |
| D7 | dev | PXD031616 DDA LiP/TrP pair | dev_sealed |
| D8 | dev | E9 repository-scale stream: 1,000 human DDA OT-HCD runs from 142 PRIDE projects, streamed through `/tmp` and deleted (added 2026-10-27, see "D8" below; dataset D8, not decision D8) | e9_d8 |
| D9 | dev | ProteoBench DDA-Astral HYE A/B x 3 + 3 human-only DDA-Astral runs (MBR negative control) + HYE FASTA; M6 LFQ dev benchmark (added 2026-10-27, see "D9" below) | d9 |
| D11 | dev | PXD003881 IonStar human + E. coli spike-in, Orbitrap Fusion DDA: group A (1x) x 4 vs group E (3x) x 4; M6d second LFQ benchmark (added 2026-10-27, see "D11" below; id skips D10 to avoid confusion with decision D10) | d11 |
| D12 | dev | PXD014415 (Lim 2019) Lumos DDA: 6 human+yeast donor runs + 3 human-only acceptor runs (yeast absent); M6d independent MBR false-transfer truth (added 2026-10-27, see "D12" below) | d12 |
| D7b | dev | PXD031616, 6 more DDA LiP/TrP runs (3 replicates per digest x condition; M5 LiP; see "D7 extension and D10cf" below) | d7b |
| D10cf | dev | PXD015406 HEK293T SEC fractions 10-38 (control series), Orbitrap Fusion DDA, IT-CID MS2; M5 CF-MS / E10 (see "D7 extension and D10cf" below; id avoids decision D10) | d10cf |
| T1 | sealed | ProteoBench DDA-QExactive, 6 files | dev_sealed |
| T2 | sealed | held-out PXD009449 pools + held-out PXD000138 pools | dev_sealed |
| T3 | sealed | PXD029252 `20210203_QExHFX3_RSLC10__DSSO_mainlib_1_5HEK_ctrl.raw` (old S4 raw; moved to `data/sealed/` on 2026-10-26; confirmed as the one-DSSO+HEK file of 6.6) | dev_sealed (raw), tier_s (truth) |
| T4 | sealed | stratified PRIDE sample, streamed | `t4_frame.tsv` |
| T5 | sealed | PXD052597 yeast DDA, 2 files | dev_sealed |

## Tier S datasets (benchmark suite, DESIGN.md section 6; S-ids are the original labels, role/set id in the `Role` column)

| ID | Role / set | Accession | Files (full names in manifest) | Source | Size | Instrument | Purpose | Notes |
|---|---|---|---|---|---|---|---|---|
| S0 | dev D0 | PXD001468 (Chick 2015, HEK293) | `b1906_293T_proteinID_01A_QE3_122212.raw`, `b1937_293T_proteinID_01B_QE3_122212.raw`, community SDRF | PRIDE FTP/HTTPS `ftp.pride.ebi.ac.uk/pride/data/archive/2015/06/PXD001468/` | 3.47 + 3.37 GB (+26 kB) | Q Exactive | Raw-speed and ID parity vs Sage (dev loop); open search | **Chosen files:** 01A is the first file in the project (the file used in the Sage/MSFragger open-search references); 01B is the matching replicate injection. Files are 3.1-3.9 GB, not 1-3 GB (no smaller ones exist). Licence: EBI ToU |
| S1 | dev D1 | PXD001077 (Vaudel 2012, *P. furiosus*) | `Velos005137.raw` | PRIDE | 0.33 GB | LTQ Orbitrap Velos | DDA ID; FDR by foreign-proteome entrapment (Pfu sample searched against Pfu+human) | Licence: EBI ToU |
| S2 | dev D2 | PXD007058 (Ferries 2017 JPR, synthetic phosphopeptide pools) | `SF_200217_pPeptideLibrary_pool{1,2,3,4}_HCDOT_rep1.raw` + submitter `pool{N}_HCDOT_rep1_txt.zip` (PD/Mascot-ptmRS exports; pool4 file is named `pool4_HCDOTrep1_txt.zip`) | PRIDE | 0.86 GB total | Orbitrap Fusion | Phosphosite localisation FLR; decoy-amino-acid FLR calibration | CC0. **Now dev (D2); truth table unavailable.** The library peptide/site list is in the JPR paper SI (ACS, not open access; no PMC copy) and is not on PRIDE. The `*_txt.zip` exports are reference identifications only, not exhaustive truth. Not obtainable without ACS access (a human could add `data/truth/PXD007058/library_sites.tsv`). Because D2 also contains PXD000138 and PXD009449 pools, FLR truth for development comes from those; for PXD007058 a proxy truth (sites supported by >=2 search engines in the submitter files + unique pool design) must be declared in RESEARCH_TREE.md before use |
| S3 | dev D3 | PXD014337 (Beveridge 2020 Nat Commun, Cas9 95-peptide XL library) | `XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw`, community SDRF | PRIDE | 0.96 GB | Q Exactive HF-X | Non-cleavable (DSS) XL; real FDR from library group design | CC0. Truth: `data/truth/IMP-X-FDR@92b2e29/peplib_Cas9_Beveridge.xlsx` (see Truth) |
| S4 | sealed T3 | PXD029252 (Matzinger 2022 Nat Commun, ribosomal 141-peptide library + HEK) | `20210203_QExHFX3_RSLC10__DSSO_mainlib_1_5HEK_ctrl.raw` (note double underscore) | PRIDE | 1.35 GB | Q Exactive HF-X | Cleavable (DSSO) XL in human background; inter vs intra real FDR | CC0. **The raw file is now sealed T3 and lives in `data/sealed/PXD029252/` (chmod 000)**; its truth/reference is in `data/truth/sealed/PXD029252/` (xiSEARCH mzid, library FASTAs, SDRF) and `data/truth/sealed/IMP-X-FDR@92b2e29/main_library.xlsx` |
| S5 | dev D5 | ProteoBench DIA-plasma PYE (Distler 2025; JPST003358) | `A9_G_DIA_nLC_tTOF_R{1,2}.d.zip`, `B9_G_DIA_nLC_tTOF_R{1,2}.d.zip` (do not rename) + `ProteoBenchFASTA_Distler_PYE.zip` | `https://proteobench.cubimed.rub.de/raws/DIA-plasma/`, `/fasta/` | 4.37 GB | timsTOF (diaPASEF, 30-min nanoLC) | Bruker DIA quant; expected ratios A:B human 1:1, yeast 1:3 (log2 -1.585), E. coli 2:1 (log2 +1) | Open access. Unzip `.d.zip` to `data/cache/` before use. The PYE search FASTA is stored under `data/truth/ProteoBench-DIA-plasma/` (small) |
| S6 | dev D6 | PXD059608 (thyroid CF-DIA, Nat Commun 2026) | `B20240606lyumg_thyroid_complex_FTC238_NO2_30minDIA_{55,56,57,58,59}.raw`, community SDRF | PRIDE | 1.06 GB | Q Exactive HF | Fraction-aware MBR guard (transfers must stay within fraction index) | CC0. **Chosen:** the 5 consecutive fractions of one 60-fraction series (FTC238, replicate NO2) with the smallest total size; they are late (high-fraction) fractions, so not representative of the full elution range. Swap to fractions 20-24 if mid-series behaviour is needed (edit manifest) |

Total Tier S download: about 15.8 GB (incl. T3, 1.35 GB) (S0 is larger than the ~10 GB estimate in lit/04 because the lit/04 total excluded S0).

## Dev and sealed sets added 2026-10-26 (manifest `data/manifests/dev_sealed.tsv`)

Selection seed: 11 (DESIGN 6.1, "pilot subsample, project folds"); the code is `scripts/data/build_dev_sealed_manifest.py`.
New download: dev 8.2 GB + sealed 27.2 GB raw (+ about 0.2 GB truth) = about 36 GB, under the 60 GB cap.

### D2 (dev): PTM / FLR pools

* **PXD000138** (Marx et al. 2013, Nat Biotechnol 31:557; 96 synthetic tryptic phospho/unmodified peptide libraries, Orbitrap Velos, HCD files `1.raw`..`96.raw`
  of 0.38-0.65 GB). Pool N = plate well A1..H12 row-major (N=5 -> A5). `random.Random(11).sample(range(1,97), 8)` = [58, 72, 60, 66, 76, 25, 24, 61].
  **Dev pools: 58 (E10), 72 (F12), 60 (E12), 66 (F6). Sealed pools (T2): 76 (G4), 25 (C1), 24 (B12), 61 (F1).** The other 88 pools are unused, not sealed. Library-entry sizes vary 660-32400 aa; sealed pool 61 (F1) is one of 12 small pools (780 aa), kept as drawn (selection rule fixed before any results).
  Truth: submitter phosphoRS exports `phosphoRS_HCD_<N>.csv` (reference IDs with site probabilities) and the library FASTA
  `009606_IPI_v3_72_with_P_Lib_062111.fasta` (IPI human + 96 pool entries `IPI:A1..H12` holding the designed peptide concatenations). The FASTA contains
  sealed pools, so it lives in `data/truth/sealed/`; `make_pool_truth.py` writes per-pool FASTAs (`data/truth/PXD000138/pool_<well>.fasta` for dev,
  `data/truth/sealed/PXD000138/` for sealed). **Gap:** exact per-pool peptide lists and site assignments are only in the Nat Biotechnol SI (not retrieved);
  library peptides have internal K/R by design, so the concatenation does not mark peptide boundaries.
* **PXD009449** (Zolg et al. 2018 Mol Cell Proteomics, ProteomeTools PTM kits; 25 pools = one synthetic-peptide kit per modification, Orbitrap Fusion Lumos).
  Only the HCD file `..._200fmol_3xHCD_R*.raw` (0.52-0.64 GB) of each chosen pool is fetched. Pools: 21 modified + 4 `*_Unmod` (Kmod/Rmod/Ymod/Pmod). The `*_Unmod` pools hold the
  unmodified counterparts of every modified pool of the same residue class and are excluded from both roles to avoid sequence leakage across the split.
  `Ymod_Phospho` is pinned to dev (phospho FLR use case); the other 20 modified pools are shuffled with `Random(11)`.
  **Dev pools: Ymod_Phospho (BA4), Kmod_Dimethyl (GA2), Rmod_Citrullin (GE1), Kmod_Formyl (GC2). Sealed pools (T2): Kmod_Acetyl (GG1), Kmod_Succinyl (GF2),
  Kmod_Ubiquitinyl (GD3, GlyGly), Kmod_Crotonyl (GA3).** Unused: the other 13 pools. Truth: per-pool MaxQuant `combined/txt/evidence.txt` pulled out of the submitter
  archive `SEARCH_<pool>.zip` (reference identifications at 1% FDR, filter `Raw file` to the 3xHCD file) and `sequences.fasta` (peptide packets per residue class; sealed folder).
  **Gap:** the per-pool synthesized-peptide table is in the MCP paper SI; the PMC copy is behind a browser proof-of-work page and the MCP site is paywalled, and
  it is not on PRIDE, Zenodo or PROSPECT (PROSPECT, zenodo 6602020, covers ProteomeTools parts 1-3 only). So truth here = reference IDs, not the full design.
* PXD007058 (S2) stays in D2; truth table unavailable (see the Tier S table).

### D3b (dev), D4 (dev), D7 (dev)

| Set | Data | Why this choice | Size | Instrument |
|---|---|---|---|---|
| D3b | PXD029252 `20210409_QExHFX2_RSLC3_DSSO_mainlib_rep2.raw` (+ `..._KK_1FDR_linkandCSM_Xi1.7.6.5.mzid`, xiSEARCH 1% reference) | DSSO main library, different file, instrument (QExHFX2 vs QExHFX3) and acquisition day than T3; has a submitter xiSEARCH result (no HEK background, unlike T3) | 1.43 GB (+41 MB) | Q Exactive HF-X |
| D4 | PXD016766 (Sage TMT benchmark data, Eclipse) `dq_00095..00106_11cell_90min_sps_{A1,A3,A5,A7,A9,A11,B1,B3,B5,B7,B9,B11}.raw` | 12 files of the 90 min SPS-MS3 TMT11 series (2.6 GB, the closest match to the "12 files, ~2.5 GB" brief; the project has 48 raw files / 16.7 GB: 180 min sps, 90 min rts-hs, 90 min hrMS2, 90 min sps). I could not verify which 12 files the Sage paper used | 2.6 GB | Orbitrap Eclipse |
| D7 | PXD031616 DDA LiP/TrP pair, 6 files `LM{2481,2486,2492,2490,2480,2487}_DDA.raw` + SDRF | See below | 8.1 GB | Orbitrap Fusion Lumos |

**D7 choice and why.** PXD031616 (Picotti lab, ETH Zurich; high-throughput LiP-MS protocol, 2022; yeast BY4742 lysate +/- human alpha-synuclein monomer M or fibrils F at 5 pmol/ug (1) or 20 pmol/ug (2), each
split into a LiP sample (Proteinase K, then trypsin) and a trypsin-only control (TrP)) has 32 `LM24xx_DDA.raw` DDA runs (1.3-1.55 GB; spectral-library DDA runs of the same samples as the
32 DIA `LM24xx.raw`, which are 6+ GB and not fetched). It is a true DDA LiP/TrP pair on one instrument, CC0, small enough, and its yeast background suits entrapment. Alternatives
rejected: LiP-Quant PXD015446 (1.95 TB; mostly DIA, only 10 `*_DDA.raw` runs), CF-MS sets such as PXD042664/CFdb (not evaluated in depth: more selection work, no paired design). The SDRF does not distinguish LiP from TrP; the
mapping below comes from the submitter MSStats reports (`R.FileName` vs `R.Condition`): LiP: F1 = 2480, 2494, 2497, 2511; F2 = 2483, 2495, 2499, 2508; M1 = 2481, 2492, 2498, 2510; M2 = 2482, 2493, 2496, 2509.
TrP: F1 = 2487, 2489, 2502, 2504; F2 = 2484, 2491, 2503, 2507; M1 = 2486, 2490, 2500, 2505; M2 = 2485, 2488, 2501, 2506 (replicate = list order). Fetched: M1 rep1 and rep2 (LiP 2481, 2492; TrP 2486, 2490) and F1 rep1
(LiP 2480; TrP 2487). Caveat: the Picotti lab (ETH) is therefore in the dev sets, so T5 (yeast) is from another lab. Truth: none (the LiP "ground truth" is the alpha-synuclein
spike-in structural change; no peptide-level table).

### Sealed sets

| Set | Data | Size | Notes |
|---|---|---|---|
| T1 | ProteoBench DDA-QExactive module: `LFQ_Orbitrap_DDA_Condition_{A,B}_Sample_Alpha_0{1,2,3}.raw` from `proteobench.cubimed.rub.de/raws/DDA/` (do not rename) + `ProteoBenchFASTA_MixedSpecies_HYE.zip` (truth/sealed) | 21.7 GB | Q Exactive HF-X; the server also offers the pre-bundled `all_data_LFQ_Quant_DDA_QExactive.tar.gz` (19 GB) and mzML, not fetched |
| T2 | PXD000138 pools 76, 25, 24, 61 (1.9 GB) + PXD009449 pools Kmod_Acetyl, Kmod_Succinyl, Kmod_Ubiquitinyl, Kmod_Crotonyl (2.2 GB), disjoint from D2 by pool; truth tables in `data/truth/sealed/` | 4.1 GB | 8 pools; see D2 for the selection |
| T3 | PXD029252 `20210203_QExHFX3_RSLC10__DSSO_mainlib_1_5HEK_ctrl.raw` | 1.35 GB | truth in `data/truth/sealed/PXD029252/` |
| T5 | PXD052597 `HFX_19152_EL_298-1_DDA_191223.raw`, `..._298-2_...` (+ SDRF, `Sample_info.xlsx` in `data/truth/sealed/PXD052597/`) | 3.1 GB | see below |

**T5 choice.** Required: yeast, DDA, Orbitrap Exploris / Q Exactive HF or newer, 2020+, a lab not in the dev sets, 1-2 runs of at most 3 GB. Candidate list came from the PRIDE API
(`/search/projects?keyword=yeast|Saccharomyces|cerevisiae`, 1452 projects, filtered to organism = S. cerevisiae only, submitted 2020+, instrument HF/HF-X/Exploris/Eclipse/Astral/Ascend, label-free/untargeted title).
Chosen **PXD052597** (PI Emmanuel Levy, submitted from the Weizmann Institute (Merbl lab): "Impact of Folded-State Protein Filaments on the Yeast Proteome"; Q Exactive HF-X, DDA, acquired 2023-12, label-free, S-trap whole-cell digest, CC0, 95 runs of
1.4-1.8 GB): label-free whole-yeast-proteome DDA, no isotope labelling, distinct lab/instrument from every dev set, not listed in `predictor_training_accessions.tsv`. Runs 298-1 and 298-2 are assumed to be two replicates of one sample (name pattern `EL_<sample>-<rep>`; `Sample_info.xlsx` not parsed). Caveats: the strains express
recombinant YFP-tagged homomers, so the protein DB needs those few sequences or accepts a handful of unmatched spectra; PXD068196 (Exploris 480, Technion) was rejected because its runs are dimethyl-labelled.

**Predictor leakage (DESIGN 6.6, owned by the leakage-check agent).** `data/manifests/predictor_training_accessions.tsv` lists PXD000138, PXD009449 (12 rows) and PXD029252 among third-party
predictor training sources, i.e. T2 and T3 (and D2, D3b) overlap with training data of some predictors; PXD052597, PXD031616 and PXD016766 do not appear. Resolve in the leakage check before the prereg tag.


## D9 (dev, added 2026-10-27 for M6): multi-species DDA LFQ benchmark that does not overlap T1

Manifest `data/manifests/d9.tsv`, hashes `data/manifests/d9.sha256`, fetch `python3 scripts/data/fetch_d9.py [--jobs 6] [--check]`
(the generic fetcher; resumable, size-verified, idempotent). 16 files, 37.9 GB (raw 22.4 GB + ProteoBench mzML 11.2 GB, used only by the
Sage/FlashLFQ comparators, + FASTA). Source: `https://proteobench.cubimed.rub.de/raws/DDA-astral/` (+ `mzml/`), `/raws/DDA-astral-entrapment/`,
`/fasta/`. Licence: open access as published by ProteoBench (no PRIDE licence: the data are not on PRIDE, see below). Downloaded at ~150 MB/s.

| What | Files | Design | Instrument |
|---|---|---|---|
| HYE benchmark (ProteoBench module "LFQ, precursor ion (Astral)", DDA) | `LFQ_Astral_DDA_15min_50ng_Condition_{A,B}_REP{1,2,3}.raw` (2.56-2.62 GB each; do not rename) | Van Puyvelde 2022 HYE composition, 50 ng, ~16 min uPAC gradient; expected log2(A/B): human 0, yeast +1, E. coli -2 (ProteoBench docs) | Orbitrap Astral (MS1 Orbitrap 180k, MS2 Astral HCD 30%, 2 Th isolation) |
| MBR negative control | `LFQ_Astral_DDA_15min_50ng_Human_0{1,2,3}.raw` (2.29-2.31 GB) from ProteoBench `DDA-astral-entrapment` | human-only digest, same instrument and method name pattern; any yeast/E. coli quantity transferred into these runs is false (Lim 2019 design) | Orbitrap Astral |
| Comparator inputs | `LFQ_Astral_DDA_15min_50ng_Condition_{A,B}_REP{1,2,3}.mzML` (ProteoBench-provided, 1.86-1.88 GB) | read by Sage and FlashLFQ only; Vulcan reads the .raw | |
| FASTA (D9-truth) | `ProteoBenchFASTA_MixedSpecies_HYE.zip` -> `data/truth/ProteoBench-DDA-Astral/` | 31,889 entries (20,537 human, 6,722 yeast, 4,401 E. coli, `Cont_` contaminants); species from the `_HUMAN/_YEAST/_ECOLI` suffix | |

**Why this set and not another.** T1 (ProteoBench DDA-QExactive = PXD028735 "Alpha" QE HF-X runs) is sealed, so D9 must not share raw files, lab
acquisition or instrument with it. Options checked on 2026-10-27:
- *PXD028735 other DDA runs*: same samples, same lab and (for QE) the same instrument series as T1, so tuning on them would leak into T1. Rejected.
- *ProteoBench DDA-Astral* (chosen): a different instrument (Orbitrap Astral), gradient (15 min uPAC) and acquisition; same HYE design, so the
  ProteoBench epsilon definition applies unchanged. Small (6 x 2.6 GB). Its companion `DDA-astral-entrapment` human-only runs give a true MBR
  false-transfer control on the same instrument and method. It is not in DESIGN §6.6's sealed list: T6 = "ProteoBench DIA-Astral plus the Astral
  entrapment module" refers to the DIA modules (`scripts/data/select_t4.py` maps T6 to PXD046444 = DIA-Astral; LITERATURE/lit/04 list only the
  DIA entrapment module, HeLa DIA). **If the T6 owner intends the DDA-astral-entrapment human runs to be part of T6, D9's MBR control must be
  replaced before T6 is opened** (recorded as an open issue).
- *PXD014777 (HYE ddaPASEF)*, *PXD003881 (IonStar)*, *PXD014415 (Lim 2019)*: larger (40-204 GB) or not HYE; PXD014415 remains the best external
  MBR false-transfer set if D9's control has to go.

**PRIDE verification.** PRIDE API v3 `/search/projects?keyword=ProteoBench` returns only PXD062685 (diaPASEF) and PXD085732; the ProteoBench docs
call the DDA-Astral data "a not-yet-released Astral DDA dataset". So D9 has no PRIDE accession and no PRIDE SHA-1; integrity rests on the
ProteoBench server's Content-Length (checked) and our sha256 (`d9.sha256`). Same situation as T1/D5 (ProteoBench-server files).

**Overlap and leakage.** No raw file is shared with T1. The HYE FASTA zip is byte-identical to T1's (`cae44ed2e3cb9ce5`): a shared search
database, not shared data. The D9 data are unreleased, so they cannot be in any predictor's training set (`predictor_training_accessions.tsv`
has no ProteoBench/Astral entry). D9 is not in the M3 static-scorer training set (13 PRIDE runs, `train.tsv`).

**D4 submitter reference (added 2026-10-27, M6).** Manifest `data/manifests/d4_ref.tsv` (+ `.sha256`), fetch
`python3 scripts/data/fetch_tier_s.py --manifest d4_ref --roles dev`: the 12 submitter SEQUEST/LDA mzIdentML files of the D4 raws
(`<id>_dq_000NN_11cell_90min_sps_<well>.mzIdentML`, 34-40 MB each, PRIDE category SEARCH, CC0) -> `data/truth/PXD016766/`. They carry
identifications (`passThreshold`) but no reporter intensities. D4 design (PRIDE sample protocol): one TMT11 plex of 11 human cell lines (RKO,
A549, U87 MG, HCT116, HEK293T, HeLa, MCF7, U2OS, SUM159, PANC1, Jurkat; channel assignment not deposited), basic-pH RP fractionated, 12 of 24
pooled fractions analysed; there are no replicate channels.

## Ground truth / reference tables

`data/truth/` is gitignored (decision D5): every file below is re-fetchable by script and its sha256 is in the manifest named in the last column
(`tier_s.sha256`, `dev_sealed.sha256`, `derived_truth.sha256`; full 64-hex digests there, 16 shown). Sealed truth stays readable under `data/truth/sealed/`.
The ProteoBench DIA-plasma PYE FASTA (`ProteoBenchFASTA_Distler_PYE.zip`) is in `data/truth/ProteoBench-DIA-plasma/`; T1's HYE FASTA is in `data/truth/sealed/ProteoBench-DDA-QExactive/`.
The IMP-X-FDR repository (`github.com/vbc-proteomics-org/imp-x-fdr`, commit `92b2e295d980ad4472e9b5a23146b82cd91900f1`) has no LICENSE file; it is fetched, not redistributed.

| Path (under `data/`) | Role / set | Source | Licence | sha256 (manifest) |
|---|---|---|---|---|
| `truth/PXD007058/pool1_HCDOT_rep1_txt.zip` | dev D2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/09/PXD007058/pool1_HCDOT_rep1_txt.zip | CC0 | `0f79fd56dc5b7110` (tier_s.sha256) |
| `truth/PXD007058/pool2_HCDOT_rep1_txt.zip` | dev D2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/09/PXD007058/pool2_HCDOT_rep1_txt.zip | CC0 | `ac8540020bf188ef` (tier_s.sha256) |
| `truth/PXD007058/pool3_HCDOT_rep1_txt.zip` | dev D2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/09/PXD007058/pool3_HCDOT_rep1_txt.zip | CC0 | `d90d2e89434f0041` (tier_s.sha256) |
| `truth/PXD007058/pool4_HCDOTrep1_txt.zip` | dev D2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/09/PXD007058/pool4_HCDOTrep1_txt.zip | CC0 | `db41a93a35259b6c` (tier_s.sha256) |
| `truth/PXD014337/PXD014337_community_annotated.sdrf.tsv` | dev D3 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/07/PXD01...munity_annotated.sdrf.tsv | CC0 | `622764f5c2324a40` (tier_s.sha256) |
| `truth/IMP-X-FDR@92b2e29/peplib_Cas9_Beveridge.xlsx` | dev D3 | https://raw.githubusercontent.com/vbc-proteomics-org/imp-x-f...eplib_Cas9_Beveridge.xlsx | none declared (no LICENSE file) | `cdfd03c9d877adf1` (tier_s.sha256) |
| `truth/sealed/PXD029252/20210203_QExHFX3_RSLC10__DSSO_mainlib_1_5HEK_ctrl_KK_1FDR_linkandCSM_Xi1.7.6.5.mzid` | sealed T3 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/05/PXD02...linkandCSM_Xi1.7.6.5.mzid | CC0 | `b3f4d071d300bd12` (tier_s.sha256) |
| `truth/sealed/PXD029252/Ribosome_EColi_SG_AND_human_abundant500excl-ribosomal_671sequences.fasta` | sealed T3 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/05/PXD02...osomal_671sequences.fasta | CC0 | `3ca8b27e55dc3ea9` (tier_s.sha256) |
| `truth/sealed/PXD029252/Ribosome_for_XL_all_contained_Proteins_from_shotgun_comb_rever.fasta` | sealed T3 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/05/PXD02..._shotgun_comb_rever.fasta | CC0 | `5b0cf8d5074da00d` (tier_s.sha256) |
| `truth/sealed/PXD029252/PXD029252_community_annotated.sdrf.tsv` | sealed T3 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/05/PXD02...munity_annotated.sdrf.tsv | CC0 | `b980a1d2f63c2e73` (tier_s.sha256) |
| `truth/sealed/IMP-X-FDR@92b2e29/main_library.xlsx` | sealed T3 | https://raw.githubusercontent.com/vbc-proteomics-org/imp-x-f...braries/main_library.xlsx | none declared (no LICENSE file) | `e109fd68d9db4e7c` (tier_s.sha256) |
| `truth/ProteoBench-DIA-plasma/ProteoBenchFASTA_Distler_PYE.zip` | dev D5 | https://proteobench.cubimed.rub.de/fasta/ProteoBenchFASTA_Distler_PYE.zip | open (ProteoBench) | `6c8bdef083f85678` (tier_s.sha256) |
| `truth/PXD000138/phosphoRS_HCD_58.csv` | dev D2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2013/05/PXD000138/phosphoRS_HCD_58.csv | EBI ToU | `1b50012f5d0f09e1` (dev_sealed.sha256) |
| `truth/PXD000138/phosphoRS_HCD_72.csv` | dev D2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2013/05/PXD000138/phosphoRS_HCD_72.csv | EBI ToU | `309139a5868f844e` (dev_sealed.sha256) |
| `truth/PXD000138/phosphoRS_HCD_60.csv` | dev D2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2013/05/PXD000138/phosphoRS_HCD_60.csv | EBI ToU | `db69c732897d0160` (dev_sealed.sha256) |
| `truth/PXD000138/phosphoRS_HCD_66.csv` | dev D2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2013/05/PXD000138/phosphoRS_HCD_66.csv | EBI ToU | `6f77c566fdede842` (dev_sealed.sha256) |
| `truth/PXD009449/Ymod_Phospho_evidence.txt` | dev D2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/06/PXD009449/SEARCH_Ymod_Phospho.zip (member `Ymod_Phospho/combined/txt/evidence.txt`) | EBI ToU | `42d1aa2b83dde09e` (dev_sealed.sha256) |
| `truth/PXD009449/Kmod_Dimethyl_evidence.txt` | dev D2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/06/PXD009449/SEARCH_Kmod_Dimethyl.zip (member `Kmod_Dimethyl/combined/txt/evidence.txt`) | EBI ToU | `5262ef3859cd60e5` (dev_sealed.sha256) |
| `truth/PXD009449/Rmod_Citrullin_evidence.txt` | dev D2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/06/PXD009449/SEARCH_Rmod_Citrullin.zip (member `Rmod_Citrullin/combined/txt/evidence.txt`) | EBI ToU | `371a8c4b7481eaa5` (dev_sealed.sha256) |
| `truth/PXD009449/Kmod_Formyl_evidence.txt` | dev D2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/06/PXD009449/SEARCH_Kmod_Formyl.zip (member `Kmod_Formyl/combined/txt/evidence.txt`) | EBI ToU | `96bb9dcc73be1a6d` (dev_sealed.sha256) |
| `truth/PXD029252/20210409_QExHFX2_RSLC3_DSSO_mainlib_rep2_KK_1FDR_linkandCSM_Xi1.7.6.5.mzid` | dev D3b | https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/05/PXD02...linkandCSM_Xi1.7.6.5.mzid | CC0 | `6315f828e207363f` (dev_sealed.sha256) |
| `truth/PXD031616/PXD031616_community_annotated.sdrf.tsv` | dev D7 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/01/PXD03...munity_annotated.sdrf.tsv | CC0 | `a39cb16ac8ae1d20` (dev_sealed.sha256) |
| `truth/sealed/PXD052597/PXD052597_community_annotated.sdrf.tsv` | sealed T5 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD05...munity_annotated.sdrf.tsv | CC0 | `a34e6a8bd954a70e` (dev_sealed.sha256) |
| `truth/sealed/PXD052597/Sample_info.xlsx` | sealed T5 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD052597/Sample_info.xlsx | CC0 | `325149ef6f5a3be0` (dev_sealed.sha256) |
| `truth/sealed/ProteoBench-DDA-QExactive/ProteoBenchFASTA_MixedSpecies_HYE.zip` | sealed T1 | https://proteobench.cubimed.rub.de/fasta/ProteoBenchFASTA_MixedSpecies_HYE.zip | open (ProteoBench) | `cae44ed2e3cb9ce5` (dev_sealed.sha256) |
| `truth/sealed/PXD000138/phosphoRS_HCD_76.csv` | sealed T2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2013/05/PXD000138/phosphoRS_HCD_76.csv | EBI ToU | `8cb9c1db34c9ea52` (dev_sealed.sha256) |
| `truth/sealed/PXD000138/phosphoRS_HCD_25.csv` | sealed T2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2013/05/PXD000138/phosphoRS_HCD_25.csv | EBI ToU | `7e8b622194fa7747` (dev_sealed.sha256) |
| `truth/sealed/PXD000138/phosphoRS_HCD_24.csv` | sealed T2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2013/05/PXD000138/phosphoRS_HCD_24.csv | EBI ToU | `0511dec652309e82` (dev_sealed.sha256) |
| `truth/sealed/PXD000138/phosphoRS_HCD_61.csv` | sealed T2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2013/05/PXD000138/phosphoRS_HCD_61.csv | EBI ToU | `29b3f218e0f64ad1` (dev_sealed.sha256) |
| `truth/sealed/PXD009449/Kmod_Acetyl_evidence.txt` | sealed T2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/06/PXD009449/SEARCH_Kmod_Acetyl.zip (member `Kmod_Acetyl/combined/txt/evidence.txt`) | EBI ToU | `6cfc2e263d2679ed` (dev_sealed.sha256) |
| `truth/sealed/PXD009449/Kmod_Succinyl_evidence.txt` | sealed T2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/06/PXD009449/SEARCH_Kmod_Succinyl.zip (member `Kmod_Succinyl/combined/txt/evidence.txt`) | EBI ToU | `a8a1c09d8e50765d` (dev_sealed.sha256) |
| `truth/sealed/PXD009449/Kmod_Ubiquitinyl_evidence.txt` | sealed T2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/06/PXD009449/SEARCH_Kmod_GlyGly.zip (member `Kmod_GlyGly/combined/txt/evidence.txt`) | EBI ToU | `f5ee7480da509c09` (dev_sealed.sha256) |
| `truth/sealed/PXD009449/Kmod_Crotonyl_evidence.txt` | sealed T2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/06/PXD009449/SEARCH_Kmod_Crotonyl.zip (member `Kmod_Crotonyl/combined/txt/evidence.txt`) | EBI ToU | `a96191830398311a` (dev_sealed.sha256) |
| `truth/sealed/PXD000138/009606_IPI_v3_72_with_P_Lib_062111.fasta` | sealed T2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2013/05/PXD00...2_with_P_Lib_062111.fasta | EBI ToU | `b5f66625ed4f7a99` (dev_sealed.sha256) |
| `truth/sealed/PXD009449/sequences.fasta` | sealed T2 | https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/06/PXD009449/sequences.fasta | EBI ToU | `153231e6c360c9a2` (dev_sealed.sha256) |
| `truth/PXD000138/pool_E10.fasta` | dev D2 | derived by `scripts/data/make_pool_truth.py` from the PXD000138 library FASTA (deterministic) | EBI ToU | `4e0239db8e5dd9a5` (derived_truth.sha256) |
| `truth/PXD000138/pool_E12.fasta` | dev D2 | derived by `scripts/data/make_pool_truth.py` from the PXD000138 library FASTA (deterministic) | EBI ToU | `2d1d011505e09941` (derived_truth.sha256) |
| `truth/PXD000138/pool_F12.fasta` | dev D2 | derived by `scripts/data/make_pool_truth.py` from the PXD000138 library FASTA (deterministic) | EBI ToU | `50a6c227468367ad` (derived_truth.sha256) |
| `truth/PXD000138/pool_F6.fasta` | dev D2 | derived by `scripts/data/make_pool_truth.py` from the PXD000138 library FASTA (deterministic) | EBI ToU | `a0bb5e04a7f43104` (derived_truth.sha256) |
| `truth/sealed/PXD000138/pool_B12.fasta` | sealed T2 | derived by `scripts/data/make_pool_truth.py` from the PXD000138 library FASTA (deterministic) | EBI ToU | `1d68a8c334fbc601` (derived_truth.sha256) |
| `truth/sealed/PXD000138/pool_C1.fasta` | sealed T2 | derived by `scripts/data/make_pool_truth.py` from the PXD000138 library FASTA (deterministic) | EBI ToU | `8396bd0814005170` (derived_truth.sha256) |
| `truth/sealed/PXD000138/pool_F1.fasta` | sealed T2 | derived by `scripts/data/make_pool_truth.py` from the PXD000138 library FASTA (deterministic) | EBI ToU | `20a71245f4d8855e` (derived_truth.sha256) |
| `truth/sealed/PXD000138/pool_G4.fasta` | sealed T2 | derived by `scripts/data/make_pool_truth.py` from the PXD000138 library FASTA (deterministic) | EBI ToU | `514222fa01a0c10b` (derived_truth.sha256) |

Truth caveat (from lit/04): the XL library-to-group assignment is defined by the xlsx tables; real FDR is computed per
the IMP-X-FDR method (Beveridge 2020; Matzinger 2022). All of `data/truth/` is gitignored (D5); reproducibility rests on the URLs and sha256 values in the manifests.

## FASTA databases (UniProt release 2026_03, released 2026-09-02, via REST `uniprotkb/stream?format=fasta`)

Canonical sequences only (no isoforms). Provenance (URL, release, retrieval date, counts, sha256) is in
`data/manifests/fasta.tsv`; the fetcher keeps existing files so the release is frozen once fetched. UniProt "current" moves on;
old releases are only available from the UniProt FTP archive, so keep a copy of `data/fasta/` if bit-exact reproduction
is needed.

| Name (`data/fasta/<name>.fasta`) | Content | Sequences | sha256 |
|---|---|---|---|
| `human_sp_canonical` | Human UP000005640, Swiss-Prot reviewed canonical (2026_03) | 20416 | `dffdff3807ec6eb7e752035b5475aaefa26d1b36e513e3344198d521a2169c42` |
| `human_proteome_canonical` | Human UP000005640 all entries (Swiss-Prot+TrEMBL), canonical (2026_03) | 147520 | `4477c2a858646a958fc84fca351b21a4f479ce4744cf3c054bc889fa69759369` |
| `human_opg` | Human UP000005640 one-protein-per-gene (reference_proteomes FTP) (2026_03) | 20652 | `9d749eb81cf7085246ce3aad7e8c23b3a9858e631fa1ddbe6cfacb439ce7af4e` |
| `pfu_UP000001013` | Pyrococcus furiosus DSM 3638 (2026_03) | 2044 | `060541ed86c636b76f7cc3f5bc5d7398b8844adee9adee3cdb5557ec91994586` |
| `ecoli_K12_UP000000625` | Escherichia coli K12 (2026_03) | 4403 | `a8052a607a98d2444c9d5dd685bd1370150df8c5570c3116811544bbf3d093a2` |
| `yeast_UP000002311` | Saccharomyces cerevisiae S288C (2026_03) | 6067 | `1aaa81860c8bda68e81e4290703d871c7c8bd12a928038772fbc91b7e17b342a` |
| `arath_UP000006548` | Arabidopsis thaliana (entrapment foreign proteome) (2026_03) | 39272 | `63edb51070a0f51cfce737f0ea84b186f2b511a5b7820a694139a04b30607464` |
| `contaminants_frankenfield2022` | Universal protein contaminants, Frankenfield et al. J Proteome Res 2022 (HaoGroup-ProtContLib @bddf6da) (git:bddf6da9334f4300caed75f7a87c096dd96b6bca) | 381 | `c0782146be6dad6ecf53e23e5a21a7ed4fea9b4802cb01e8d00091404ec66926` |

Primary human DB for searches: `human_sp_canonical` (Swiss-Prot reviewed, UP000005640) unless an experiment states otherwise;
`human_opg` is the one-protein-per-gene reference proteome. Arabidopsis is the foreign proteome for entrapment
(I/L-equivalent overlap with the target DB must be removed at entrapment-construction time, DESIGN.md section 5).
Contaminants: Frankenfield et al., J Proteome Res 2022, "Universal protein contaminant" FASTA from
`github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics`, commit
`bddf6da9334f4300caed75f7a87c096dd96b6bca`, file `0602_Universal Contaminants.fasta` (repo has no LICENSE file; the paper is open access).

## Per-file raw inventory (verified 2026-10-27)

Role and set id follow DESIGN 6.6 and the manifest columns. Sealed rows (mode 000) are size-checked by stat only; their sha256 was recorded right after download in `sealed.sha256` and is re-verifiable only through `bench/sealed_run.py --verify`.
Dataset-level licences: PRIDE metadata, re-read 2026-10-27 (CC0 for PXD007058, PXD014337, PXD016766, PXD029252, PXD031616, PXD052597, PXD059608; EBI terms of use for PXD000138, PXD001077, PXD001468, PXD009449).
PRIDE SHA-1 cross-check (`results/2026-10-27_data_verify/pride_sha1.tsv`): 37 files match, 0 mismatch; 28 files (all of PXD001468, PXD001077, PXD000138 and the 2025 SDRFs/metadata) have no checksum on PRIDE and rest on size plus our sha256 only.

| Path (under `data/`) | Set | Role | Size (B) | Licence | sha256 (manifest) |
|---|---|---|---|---|---|
| `raw/PXD001468/b1906_293T_proteinID_01A_QE3_122212.raw` | D0 | dev | 3470679015 | EBI ToU | `e00bf1518b5fc26b` (tier_s.sha256) |
| `raw/PXD001468/b1937_293T_proteinID_01B_QE3_122212.raw` | D0 | dev | 3368964863 | EBI ToU | `2b70f3f288897047` (tier_s.sha256) |
| `raw/PXD001468/PXD001468_community_annotated.sdrf.tsv` | D0 | dev | 26349 | EBI ToU | `32f1a34b49bbd774` (tier_s.sha256) |
| `raw/PXD001077/Velos005137.raw` | D1 | dev | 334619250 | EBI ToU | `c4a43e3423f1c892` (tier_s.sha256) |
| `raw/PXD007058/SF_200217_pPeptideLibrary_pool1_HCDOT_rep1.raw` | D2 | dev | 202338736 | CC0 | `45953b4421b19c61` (tier_s.sha256) |
| `raw/PXD007058/SF_200217_pPeptideLibrary_pool2_HCDOT_rep1.raw` | D2 | dev | 218563894 | CC0 | `20c69f3f27aa6eef` (tier_s.sha256) |
| `raw/PXD007058/SF_200217_pPeptideLibrary_pool3_HCDOT_rep1.raw` | D2 | dev | 218366958 | CC0 | `47cf724503fbad96` (tier_s.sha256) |
| `raw/PXD007058/SF_200217_pPeptideLibrary_pool4_HCDOT_rep1.raw` | D2 | dev | 210734906 | CC0 | `30cd5ca4ee5c4063` (tier_s.sha256) |
| `raw/PXD014337/XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw` | D3 | dev | 960697429 | CC0 | `128bf16054a40ba1` (tier_s.sha256) |
| `raw/ProteoBench-DIA-plasma/A9_G_DIA_nLC_tTOF_R1.d.zip` | D5 | dev | 1112972541 | open (ProteoBench) | `ce83c4cfa6730cb1` (tier_s.sha256) |
| `raw/ProteoBench-DIA-plasma/A9_G_DIA_nLC_tTOF_R2.d.zip` | D5 | dev | 1148906243 | open (ProteoBench) | `7c9b48fb9cafd9ce` (tier_s.sha256) |
| `raw/ProteoBench-DIA-plasma/B9_G_DIA_nLC_tTOF_R1.d.zip` | D5 | dev | 1138659226 | open (ProteoBench) | `035cbab61ac80121` (tier_s.sha256) |
| `raw/ProteoBench-DIA-plasma/B9_G_DIA_nLC_tTOF_R2.d.zip` | D5 | dev | 963587638 | open (ProteoBench) | `becf6a039fc27dd2` (tier_s.sha256) |
| `raw/PXD059608/B20240606lyumg_thyroid_complex_FTC238_NO2_30minDIA_55.raw` | D6 | dev | 207201600 | CC0 | `df9c6fd9f87b6ff3` (tier_s.sha256) |
| `raw/PXD059608/B20240606lyumg_thyroid_complex_FTC238_NO2_30minDIA_56.raw` | D6 | dev | 226705766 | CC0 | `e0ab1632a2f42cd2` (tier_s.sha256) |
| `raw/PXD059608/B20240606lyumg_thyroid_complex_FTC238_NO2_30minDIA_57.raw` | D6 | dev | 235499737 | CC0 | `ac2fadf557c6765d` (tier_s.sha256) |
| `raw/PXD059608/B20240606lyumg_thyroid_complex_FTC238_NO2_30minDIA_58.raw` | D6 | dev | 185638012 | CC0 | `0bb37793402a8b7d` (tier_s.sha256) |
| `raw/PXD059608/B20240606lyumg_thyroid_complex_FTC238_NO2_30minDIA_59.raw` | D6 | dev | 208162800 | CC0 | `e8553da0699d77b2` (tier_s.sha256) |
| `raw/PXD059608/PXD059608_community_annotated.sdrf.tsv` | D6 | dev | 352269 | CC0 | `a4a931479d7a5c9d` (tier_s.sha256) |
| `raw/PXD000138/58.raw` | D2 | dev | 417126672 | EBI ToU | `1e329ff2fa2a2884` (dev_sealed.sha256) |
| `raw/PXD000138/72.raw` | D2 | dev | 495419780 | EBI ToU | `1e49865e525e424f` (dev_sealed.sha256) |
| `raw/PXD000138/60.raw` | D2 | dev | 460032922 | EBI ToU | `bc2ed7556826eeae` (dev_sealed.sha256) |
| `raw/PXD000138/66.raw` | D2 | dev | 495419462 | EBI ToU | `20f508c4dccc459c` (dev_sealed.sha256) |
| `raw/PXD009449/02330a_BA4_4087_1_PTM_TrainKit_Ymod_Phospho_200fmol_3xHCD_R1.raw` | D2 | dev | 636373246 | EBI ToU | `e9d88487a7cf563b` (dev_sealed.sha256) |
| `raw/PXD009449/02330a_GA2_3990_09_PTM_TrainKit_Kmod_Dimethyl_200fmol_3xHCD_R1.raw` | D2 | dev | 552477520 | EBI ToU | `633a1f7391289a2e` (dev_sealed.sha256) |
| `raw/PXD009449/02330a_GE1_3990_05_PTM_TrainKit_Rmod_Citrullin_200fmol_3xHCD_R1.raw` | D2 | dev | 518683564 | EBI ToU | `38b8684887755691` (dev_sealed.sha256) |
| `raw/PXD009449/02330a_GC2_3990_11_PTM_TrainKit_Kmod_Formyl_200fmol_3xHCD_R1.raw` | D2 | dev | 521663368 | EBI ToU | `6cb4049303f1f203` (dev_sealed.sha256) |
| `raw/PXD029252/20210409_QExHFX2_RSLC3_DSSO_mainlib_rep2.raw` | D3b | dev | 1425918096 | CC0 | `a983232bf1fb3c95` (dev_sealed.sha256) |
| `raw/PXD016766/dq_00095_11cell_90min_sps_A1.raw` | D4 | dev | 209342012 | CC0 | `36cb2e1a95ab1d32` (dev_sealed.sha256) |
| `raw/PXD016766/dq_00096_11cell_90min_sps_A3.raw` | D4 | dev | 203822033 | CC0 | `7b24441bf351bd92` (dev_sealed.sha256) |
| `raw/PXD016766/dq_00097_11cell_90min_sps_A5.raw` | D4 | dev | 221882856 | CC0 | `b3de450ecd39da54` (dev_sealed.sha256) |
| `raw/PXD016766/dq_00098_11cell_90min_sps_A7.raw` | D4 | dev | 234723755 | CC0 | `d3d344f37a7a6182` (dev_sealed.sha256) |
| `raw/PXD016766/dq_00099_11cell_90min_sps_A9.raw` | D4 | dev | 235081101 | CC0 | `35dbc827dd6e6e1f` (dev_sealed.sha256) |
| `raw/PXD016766/dq_00100_11cell_90min_sps_A11.raw` | D4 | dev | 231970426 | CC0 | `8ca650a9eeac1896` (dev_sealed.sha256) |
| `raw/PXD016766/dq_00101_11cell_90min_sps_B1.raw` | D4 | dev | 229691147 | CC0 | `3f1d3e7c1f3792c1` (dev_sealed.sha256) |
| `raw/PXD016766/dq_00102_11cell_90min_sps_B3.raw` | D4 | dev | 209226080 | CC0 | `7971ebbb3868d043` (dev_sealed.sha256) |
| `raw/PXD016766/dq_00103_11cell_90min_sps_B5.raw` | D4 | dev | 209090667 | CC0 | `f74cf5c8a708a732` (dev_sealed.sha256) |
| `raw/PXD016766/dq_00104_11cell_90min_sps_B7.raw` | D4 | dev | 208561029 | CC0 | `27f0a5d775f9c7b4` (dev_sealed.sha256) |
| `raw/PXD016766/dq_00105_11cell_90min_sps_B9.raw` | D4 | dev | 215470314 | CC0 | `214e7d895515687f` (dev_sealed.sha256) |
| `raw/PXD016766/dq_00106_11cell_90min_sps_B11.raw` | D4 | dev | 232375511 | CC0 | `9c8828d1aaca20d8` (dev_sealed.sha256) |
| `raw/PXD031616/LM2481_DDA.raw` | D7 | dev | 1365148534 | CC0 | `ee7b251ad2b89289` (dev_sealed.sha256) |
| `raw/PXD031616/LM2486_DDA.raw` | D7 | dev | 1279637792 | CC0 | `35d3d49a617ffd18` (dev_sealed.sha256) |
| `raw/PXD031616/LM2492_DDA.raw` | D7 | dev | 1348489012 | CC0 | `21996f72145cfd56` (dev_sealed.sha256) |
| `raw/PXD031616/LM2490_DDA.raw` | D7 | dev | 1270724976 | CC0 | `9a26f15f2d427956` (dev_sealed.sha256) |
| `raw/PXD031616/LM2480_DDA.raw` | D7 | dev | 1477148524 | CC0 | `5dfcaddb638e03a4` (dev_sealed.sha256) |
| `raw/PXD031616/LM2487_DDA.raw` | D7 | dev | 1355740750 | CC0 | `e65fe3a135d8a46b` (dev_sealed.sha256) |
| `sealed/PXD052597/HFX_19152_EL_298-1_DDA_191223.raw` | T5 | sealed | 1532412808 | CC0 | `78ff7f484b842cf9` (sealed.sha256) |
| `sealed/PXD052597/HFX_19152_EL_298-2_DDA_191223.raw` | T5 | sealed | 1558908720 | CC0 | `a7e371f3282f0cf6` (sealed.sha256) |
| `sealed/ProteoBench-DDA-QExactive/LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.raw` | T1 | sealed | 3459277876 | open (ProteoBench) | `bc7c31ee9188941a` (sealed.sha256) |
| `sealed/ProteoBench-DDA-QExactive/LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_02.raw` | T1 | sealed | 3612100302 | open (ProteoBench) | `df3ebbf7c47237b3` (sealed.sha256) |
| `sealed/ProteoBench-DDA-QExactive/LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_03.raw` | T1 | sealed | 3701529661 | open (ProteoBench) | `c1f0d2174502ae1a` (sealed.sha256) |
| `sealed/ProteoBench-DDA-QExactive/LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_01.raw` | T1 | sealed | 3501153502 | open (ProteoBench) | `6f6f5391c2c57773` (sealed.sha256) |
| `sealed/ProteoBench-DDA-QExactive/LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_02.raw` | T1 | sealed | 3904127239 | open (ProteoBench) | `190b9dcf3bbe246b` (sealed.sha256) |
| `sealed/ProteoBench-DDA-QExactive/LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_03.raw` | T1 | sealed | 3693116392 | open (ProteoBench) | `e74af37fa94756ad` (sealed.sha256) |
| `sealed/PXD000138/76.raw` | T2 | sealed | 483033470 | EBI ToU | `85042f86c0ac1446` (sealed.sha256) |
| `sealed/PXD000138/25.raw` | T2 | sealed | 436181432 | EBI ToU | `a5b6168885f5ef7f` (sealed.sha256) |
| `sealed/PXD000138/24.raw` | T2 | sealed | 439844190 | EBI ToU | `9331e13d83e13975` (sealed.sha256) |
| `sealed/PXD000138/61.raw` | T2 | sealed | 439875376 | EBI ToU | `5a41ca7aef3fe00d` (sealed.sha256) |
| `sealed/PXD009449/02330a_GG1_3990_07_PTM_TrainKit_Kmod_Acetyl_200fmol_3xHCD_R1.raw` | T2 | sealed | 540363056 | EBI ToU | `d96a7ef92102b817` (sealed.sha256) |
| `sealed/PXD009449/02330a_GF2_3991_03_PTM_TrainKit_Kmod_Succinyl_200fmol_3xHCD_R1.raw` | T2 | sealed | 536698096 | EBI ToU | `25ac25b7663d1b3d` (sealed.sha256) |
| `sealed/PXD009449/02330a_GD3_3991_09_PTM_TrainKit_Kmod_Ubiquitinyl_200fmol_3xHCD_R1.raw` | T2 | sealed | 571151246 | EBI ToU | `676216923e64a48b` (sealed.sha256) |
| `sealed/PXD009449/02330a_GA3_3991_06_PTM_TrainKit_Kmod_Crotonyl_200fmol_3xHCD_R1.raw` | T2 | sealed | 545493856 | EBI ToU | `03f74feda24b3e60` (sealed.sha256) |
| `sealed/PXD029252/20210203_QExHFX3_RSLC10__DSSO_mainlib_1_5HEK_ctrl.raw` | T3 | sealed | 1350506891 | CC0 | `c217d01665d4e7ac` (sealed.sha256) |

## Open issues and how they were resolved (2026-10-27)

* **PXD036062 (Messner et al., Cell 2023, yeast KO proteomes, 4,699 strains, SWATH/DIA)**: resolved as to why it 404s. It is a ProteomeXchange dataset hosted by **MassIVE**, not PRIDE
  (data availability statement: "deposited to the ProteomeXchange Consortium via the MassIVE repository, PXD036062"), so the PRIDE API cannot see it. Not downloaded, since it is not in any dev or sealed set and is a possible M8 cohort-scale set only; the MSV id was not retrieved.
* **PXD056598 (PYE plasma, Distler et al., Nat Commun 16:8774, 2025, 1,116 runs, 12 sites)**: **unresolved**. PRIDE API and ProteomeCentral return nothing, and the paper text was not reachable here, so I could not confirm the repository. The subset we actually use (D5) comes from ProteoBench (`proteobench.cubimed.rub.de/raws/DIA-plasma/`, also jPOST JPST003358) and is complete. Do not plan around PXD056598.
* **PXD080610 (ProteoBench low-input, per lit/02)**: **unresolved**; 404 on PRIDE and no ProteomeCentral record, possibly unreleased or mis-cited. It is not in any role.
* **InfinDIA accession to condition mapping**: partly verified against the DIA-NN README table (github.com/vdemichev/DiaNN, InfinDIA benchmark) and PRIDE metadata.

  | Accession | DIA-NN README condition (pairs/min, millions) | PRIDE check | Verdict |
  |---|---|---|---|
  | PXD005573 | QE HF-X 65 min K-GG human + E. coli (440); QE HF 30 min (315) | Bruderer DIA-optimisation project, instrument Q Exactive HF, organisms human + mouse, keyword "Human dia proteomics"; no K-GG or E. coli in metadata | not confirmed; README label doubtful (the 65 min K-GG rows may come from another project) |
  | PXD044991 | Astral 40 SPD (215; unspecific digest 210) | "One-Tip", PRIDE instrument tag QE HF but file names `20230721_AST1_EVO5_..._1h_40SPD_*.raw` and keyword Astral | confirmed Astral (PRIDE tag is wrong) |
  | PXD017703 | timsTOF 60 SPD (51) | diaPASEF, timsTOF Pro, yeast + human, Evosep files present | confirmed |
  | PXD034222 | timsTOF 70 min (29) | Oenothera/senescent fibroblasts DIA, timsTOF Pro, `.d.zip` | confirmed |
  | MSV000093613 | Astral 30 min phospho (275) | MassIVE, not checked | open |

  Treat PXD005573 conditions as unverified and MSV000093613 as unchecked until someone inspects the raw files. None of these is a dev or sealed set.

## Known gaps in the data (not resolvable without a human)

* PXD007058 library site table (JPR SI, ACS paywall), PXD009449 per-pool synthesized peptide table (MCP SI), PXD000138 per-pool peptide lists (Nat Biotechnol SI): truth is limited to submitter reference identifications and library FASTAs (see the D2 notes).
* IMP-X-FDR has no licence; PXD052597 `Sample_info.xlsx` is unparsed (298-1/298-2 replicate relation assumed).
* T4 is streamed and has no entry in the verifier (`t4_selection.tsv` is a selection list only).

## Fixed seeds (copied verbatim from DESIGN.md section 6.1)

| Purpose | Seeds |
|---|---|
| General | `SEEDS = [11, 23, 47, 83, 131]` |
| Decoy shuffles | `DECOY_SEEDS = [1001, 1002, 1003, 1004, 1005]` |
| Entrapment, dev | `ENTRAP_SEEDS_DEV = [2001, 2002, 2003]` |
| Entrapment, sealed | `ENTRAP_SEEDS_SEALED = [2004, 2005, 2006]` |
| Calibration decoys | `CALIB_DECOY_SEED = 3001` |
| Pilot subsample, project folds | 11 |

Replacement shuffles are hash(peptide, seed), so they do not depend on DB order. Every entrapment peptide must be unique against all targets and all decoys of all decoy seeds.
The default pseudo-reverse decoy is deterministic (DESIGN 6.2). Dev/sealed pool selection in this file used seed 11.

## Held-out and evaluation sets

The frozen dev/sealed split is the table at the top of this file (DESIGN 6.6); no further placeholder remains here.

- Evaluation datasets excluded from training by accession: all dev and sealed accessions above (PXD001468, PXD001077, PXD007058, PXD014337, PXD029252, PXD059608, PXD000138, PXD009449, PXD016766, PXD031616, PXD052597,
  ProteoBench DIA-plasma and DDA-QExactive), plus the T4 frame accessions drawn by `scripts/data/select_t4.py` (`data/manifests/t4_frame.tsv`) and the Tier M/L accessions listed in DESIGN.md section 6.
- Project-level train/val/test fold assignment uses seed 11 and is written by the training agent into the training manifests (not here).

## Training data (role = `train`; M3 static-first scorer, added 2026-10-27)

Static-first scorer training runs (DESIGN §4 L3). Selected by `scripts/m3/select_train.py` (seed 11) from the eligible T4 frame (`t4_frame.tsv`) with hard exclusion of: every T4 *selection* accession (`t4_selection.tsv`), every dev/sealed accession (`tier_s.tsv`, `dev_sealed*.tsv`), every predictor-training accession (`predictor_training_accessions.tsv`), yeast (memorization arm M / T5) and Arabidopsis (arm-C entrapment proteome). Instrument filter: Q Exactive family / Exploris (OT-HCD); label-free, no enrichment; filenames suggesting phospho/TMT/QC/blank/DIA excluded; raw 0.3–1.3 GB; one run per project.
Fetched by `scripts/m3/fetch_train.py` into `data/raw/train/<acc>/` (resumable, size-verified); sha256 in `data/manifests/train.sha256`. Project-disjoint fold map (3 folds, round-robin over the seeded project order): `data/manifests/train_folds.tsv` (committed). Licences: per PRIDE project metadata (not re-checked per project; research use only).
**Rule for later T4 redraws:** these 13 projects are training data and must be excluded from any future T4 sample (they are in the frame, not in the current selection).
Excluded after ingest: TR10 `MN-WT1-Frac8.raw` (PXD019362): MS2 are ion-trap CID (frame metadata said Q Exactive HF), so it is not OT-HCD; kept on disk, not used for training.

| id | role | accession | species | fold | file | size (B) | sha256 (train.sha256) |
|---|---|---|---|---|---|---|---|
| TR00 | train | PXD068465 | human | 0 | `DAI-1333_hCells_WT-3_STrapMicro_PSRTandC_90m_20240617_QE.raw` | 1050266312 | `124ca91d3c36b0c4` |
| TR01 | train | PXD026663 | human | 1 | `KA121LQ_CSF_65.raw` | 645067665 | `3c59e653cb2d5f2c` |
| TR02 | train | PXD074418 | human | 2 | `30min_LPS_2.raw` | 1283744578 | `a152b50a73ca40ad` |
| TR03 | train | PXD044783 | human | 0 | `PT10402-292_190921044254.raw` | 1278946759 | `dd723e84e30e923a` |
| TR04 | train | PXD009743 | human | 1 | `Ctrl_JH_repli1_SCX4.raw` | 916293263 | `b795910656aa92e5` |
| TR05 | train | PXD047664 | human | 2 | `P22_081_Singapur_Smim26_SampleII_02.raw` | 580640937 | `628ffe86e9a4c449` |
| TR06 | train | PXD047868 | human | 0 | `No_guide_control_2_promoter_reporter.raw` | 1275063819 | `21179a3fe4c11932` |
| TR07 | train | PXD045252 | human | 1 | `c3.raw` | 515478016 | `11f95f81982abbff` |
| TR08 | train | PXD031264 | human | 2 | `211221_LL_17183_1.raw` | 547614046 | `013ecaf7ea935cf3` |
| TR09 | train | PXD059835 | mouse | 0 | `210521-MClass-21_0844-14-BAZ.raw` | 1109102212 | `f77adca474d26521` |
| TR10 | train | PXD019362 | mouse | 1 | `MN-WT1-Frac8.raw` | 1023762710 | `9c2e302232f69968` |
| TR11 | train | PXD019012 | ecoli | 2 | `XS_Shura_High_pH_Frac12_api05_4.raw` | 1158541821 | `fe940c8408007acb` |
| TR12 | train | PXD029006 | ecoli | 0 | `005_ECOLI72_F20.raw` | 600901192 | `f861aa0145c1b975` |
| FASTA | train | UniProt Swiss-Prot mouse (organism 10090, reviewed, canonical; release of 2026-10-27) | mouse | - | `fasta/mouse_sp_canonical.fasta` | - | `313a03a01c364337` |

## DN (dev, added 2026-10-27 for M3e / Amendment A1.3): wrong-species G2 null sources

Purpose: null spectra for G2 / tail-G2 / G1 / G3 / G5 (DESIGN §6.4, §6.0 A1.3). They are never searched against their own proteome, are never used for training or selection, and appear only as null-provenance spectra inside the G2 mixtures (scripts/m3e/g2mix.py).
Selection (`scripts/m3e/select_nulls.py`, seed 11, fixed before fetching; log `data/manifests/g2null_select.log.json`):
* PRIDE frame `t4_frame.tsv`; single organism from an allow-list of prokaryotes and unicellular eukaryotes far from human; DDA, LFQ, no enrichment, Thermo Q Exactive/Exploris (OT-HCD).
* Hard-excluded: the T4 selection, all dev/sealed/tier-S accessions, the 13 training projects, predictor-training accessions, and the E9/D8 stream. Organisms already used (P. furiosus D1, S. cerevisiae D7), training species (human, mouse, E. coli) and Arabidopsis (arm C) are also excluded.
* 4 projects × 3 runs; files 0.3–1.5 GB.
* Together with D1 (Pfu) and D7 (yeast, 6 runs) this gives 6 null projects (A1.3 requires ≥ 3).
**Rule for T4 redraws:** these 4 projects are dev data and must be excluded from any future T4 sample. They are non-human, so they cannot enter the human E9/D8 stream.
Fetch: `python3 scripts/m3e/fetch_nulls.py`. RAWs go to `data/raw/<acc>/` (manifest `data/manifests/g2null.tsv`, sha256 `data/manifests/g2null.sha256`). Homolog-filter proteomes go to `data/fasta/dn_<upid>.fasta` (sha256 `data/manifests/g2null_fasta.sha256`; UniProt REST, 2026-10-27).
Correction made before any DN file was fetched: the free-text proteome lookup returned *Synechococcus* WH7803 for PXD021787. The project strain is PCC 7002, so the lookup was replaced by UP000001688. Its UniProtKB proteome query is empty, so the sequences come from UniParc `upid:UP000001688` (3,191). PXD023811's strain is not specified, so USA300 ∪ NCTC 8325 is used.
Licences: per PRIDE project metadata (not re-checked per project; research use only).

| id | role | accession | organism | file | size (B) | sha256 (g2null.sha256) |
|---|---|---|---|---|---|---|
| DN00 | dev | PXD015648 | T. brucei | `20170221_QEP2_FBU_FB_collabJanzen_2T1p99xC4x48_2.raw` | 689185737 | `03d9fd7c84950804` |
| DN01 | dev | PXD015648 | T. brucei | `20170221_QEP2_FBU_FB_collabJanzen_2T1p99xC4x72_4.raw` | 791651730 | `5ee89a3d2c3190d8` |
| DN02 | dev | PXD015648 | T. brucei | `20170221_QEP2_FBU_FB_collabJanzen_2T1p99xC9x72_3.raw` | 734444296 | `d6da5273938ae117` |
| DN03 | dev | PXD023811 | S. aureus | `374-A.raw` | 1412798221 | `356cac1ed7953cce` |
| DN04 | dev | PXD023811 | S. aureus | `536-A.raw` | 1295044249 | `d6c0d0343a4e643c` |
| DN05 | dev | PXD023811 | S. aureus | `860-A2.raw` | 1415610074 | `5280b83299b51758` |
| DN06 | dev | PXD021787 | Synechococcus sp. PCC 7002 | `TL190218_KDahlgren_CameronLab_0p5ug_GFP-APEX-B.raw` | 457252077 | `912d0c5dc86fd05c` |
| DN07 | dev | PXD021787 | Synechococcus sp. PCC 7002 | `TL190621_KDahlgren_CameronLab_250ng_GFP-APEX-A.raw` | 1023430429 | `08d97c268f98d984` |
| DN08 | dev | PXD021787 | Synechococcus sp. PCC 7002 | `TL190621_KDahlgren_CameronLab_250ng_psbU-APEX-B.raw` | 772146199 | `9f0bfbc03c2dc7c7` |
| DN09 | dev | PXD055392 | B. cereus ATCC 14579 | `Q15827_MS23-053_6.raw` | 1113700662 | `fa67d5435b20e477` |
| DN10 | dev | PXD055392 | B. cereus ATCC 14579 | `Q15836_MS23-053_11.raw` | 1131961908 | `8dfaeafc194c3235` |
| DN11 | dev | PXD055392 | B. cereus ATCC 14579 | `Q15837_MS23-053_12.raw` | 1034421410 | `dcbe61743907d0eb` |

## D11 (dev, added 2026-10-27 for M6d): second DDA LFQ benchmark, human + E. coli spike-in (IonStar, Orbitrap Fusion)

The ids D11 and D12 skip "D10" so that they are not confused with decision D10 in RESEARCH_TREE.md. Manifest `data/manifests/d11.tsv`,
hashes `data/manifests/d11.sha256`, fetch `python3 scripts/data/fetch_d11.py [--jobs 5] [--check]` (the generic fetcher; resumable,
size-verified, PRIDE SHA-1 checked on arrival). 9 files, 11.92 GB (8 raw, 1.47-1.50 GB each, + 15 MB FASTA) in `data/raw/PXD003881/` and
`data/truth/PXD003881/`. Source: PRIDE `https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/05/PXD003881/`. Licence: EBI terms of use.

| What | Files | Design | Instrument |
|---|---|---|---|
| Condition A (1x) | `B03_{10,11,20,21}_150304_human_ecoli_A_3ul_3um_column_95_HCD_OT_2hrs_30B_9B.raw` | E. coli (DH5α) lysate spiked into human Panc-1 lysate at 3% of total protein, then alkylated/precipitated/digested; 2 µg peptide per injection | Orbitrap Fusion, DDA, MS1 Orbitrap 120k, MS2 Orbitrap 15k HCD 35%; 100 cm x 75 µm column, 2.5 h gradient (PRIDE text; file names say "2hrs") |
| Condition E (3x) | `B03_{05,06,15,16}_150304_human_ecoli_E_3ul_3um_column_95_HCD_OT_2hrs_30B_9B.raw` | same, E. coli at 9% of total protein | same method, same batch (B03, 2015-03-04) |
| FASTA (D11-truth) | `Human_ecoli_trypsin_1501v_uniprot_sprot.fasta` -> `data/truth/PXD003881/` | submitter database (human + E. coli SwissProt 2015-01); reference only | |

**Expected ratios (log2 E/A; note the direction: E is the high-spike condition).** Source: the PRIDE PXD003881 project description
(sample processing protocol, verified via API v3 on 2026-10-27): "DH5α E. Coli lysate were spiked in human pancreatic cancer cells (Panc-1)
lysate at five different levels: 1 fold (3% E.Coli), 1.5 fold (4.5% E.Coli), 2 fold (6% E.Coli), 2.5 fold (7.5% E.Coli) and 3 fold (9%
E.Coli)", 2 µg peptides per injection; group letters A-E = levels 1x-3x in that order (file names `human_ecoli_<A..E>`). With equal total
load: **E. coli log2(E/A) = log2(3) = +1.585; human log2(91/97) = -0.092**. If intensities are normalised so that the human median is 0 (the
usual LFQ normalisation), the expected E. coli value is log2(3 x 97/91) = **+1.677**. State which convention a metric uses. Paper: Shen et
al. 2018 PNAS 115:E4767 (doi:10.1073/pnas.1800541115; the full text was not machine-readable here, so the PRIDE text is the cited source).
lit/04 describes the design as "human proteins spiked in E. coli", which is reversed; PRIDE says E. coli into human.

**Why A vs E and these runs.** A vs E is the largest fold change in the set (3x, close to D9's yeast 2x / E. coli 4x). The project has 4 runs
per level (20 raw, 40.4 GB, + mzXML/mgf/mzid). We take all 4 runs of A and E (8 runs, 11.9 GB), not a 3+3 subset, so no run choice is needed.
Runs were acquired end-over-end alternating (A1 B1 C1 D1 E1 E2 D2 ... per PRIDE), so A and E are interleaved in time, not blocked.
The other levels (B, C, D: 1.5x, 2x, 2.5x) can be added later by editing the manifest (same file pattern, B03_02/09/12/19 = B, 03/08/13/18 =
C, 04/07/14/17 = D).

**Candidates checked and rejected (2026-10-27, PRIDE API v3 keyword search: HYE, three species, mixed proteome, LFQbench, two-proteome,
hybrid proteome, species mixture, E. coli spike-in, benchmark label-free).**
- *PXD007683* (O'Connell 2018 JPR, human + yeast 10/5/3.3%, Lumos): the 11 LFQ raws (`a05191..a05202`) have no run-to-group mapping in PRIDE
  or the community SDRF, so the expected ratio per run cannot be stated from metadata. Rejected.
- *PXD017812/PXD017865/PXD017915* (UHR-IonStar, Wang 2021 Anal Chem, HYE, Lumos 240k MS1): SDRF gives groups A-G but the spike levels are only in
  paper Table S1 (not obtainable here). Larger (47-77 GB). Rejected; a good HYE upgrade if a human fetches Table S1.
- *PXD045500 / PXD039576* (Matzinger 2024, HYE 8:1:1, Exploris 480 FAIMS): a single HYE mix, no second condition. Rejected.
- *PXD002952* (LFQbench, TripleTOF), *PXD053241* (CHIMERYS; re-analyses PXD028735), *PXD070049/PXD071205* (Astral/timsTOF/ZenoTOF): wrong
  instrument class or contain T1 data. PXD028735 is excluded outright (contains sealed T1).

**PRIDE verification.** Project PXD003881 (API v3): title "Quantitative Proteomics Benchmark Dataset to evaluate label-free quantitative
methods ... E coli proteomes spiked-in Human proteins at 5 different levels (N=20)", instrument Orbitrap Fusion, organisms E. coli + human,
submitted 2016-03-31, public 2018-05-16. All 100 project files carry a PRIDE SHA-1; every D11 file was checked against it on arrival and again
after download: **9/9 SHA-1 match, 0 mismatch**, 9/9 sizes match PRIDE, sha256 in `d11.sha256` (downloaded 2026-10-27, ~15-30 MB/s).

**Overlap and leakage.** PXD003881 is not in `dev_sealed`, `dev_sealed_extra`, `tier_s`, `t4_selection`, `train`, `e9_d8`, `g2null`, `d9` or
`predictor_training_accessions.tsv` (grep, 2026-10-27). It is in the T4 crawl frame (`t4_frame` / `all_projects_classified.tsv`) but was not
drawn into T4. Different lab (Qu lab, Buffalo), instrument (Fusion) and species pair from T1 and D9. Absence from the predictor list means "not
listed", not proven absence from every predictor's training data.

**FASTA.** Use the ProteoBench HYE FASTA (`data/truth/ProteoBench-DDA-Astral/ProteoBenchFASTA_MixedSpecies_HYE.zip`, unzipped copy at
`/tmp/vulcan_m6/fasta/ProteoBenchFASTA_MixedSpecies_HYE.fasta`): 20,537 human + 6,722 yeast + 4,401 E. coli SwissProt + contaminants. It covers
both species (DH5α is an E. coli K-12 derivative; the SwissProt E. coli set is K-12). The yeast entries act as an extra entrapment proteome
(yeast is absent from every D11 run). The submitter FASTA is kept only for reference.

## D12 (dev, added 2026-10-27 for M6d): independent MBR false-transfer truth, human-only vs human + yeast (Lim 2019, Lumos)

Manifest `data/manifests/d12.tsv`, hashes `data/manifests/d12.sha256`, fetch `python3 scripts/data/fetch_d12.py [--jobs 5] [--check]`.
10 files, 2.83 GB (9 raw, 0.30-0.32 GB each, + community SDRF) in `data/raw/PXD014415/` and `data/truth/PXD014415/`. Source: PRIDE
`https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/09/PXD014415/`. Licence: CC0. Not ProteoBench, not D9, not T1.

| Role | Files | Sample | Notes |
|---|---|---|---|
| Donor (yeast present) | `a{11869..11874}_human_yaeast_90min_hrMS2.raw` (6; "yaeast" is the submitter's spelling) | HY: human HCT116 digest + yeast digest at 10% of human mass (3 µg yeast into 30 µg human) | first 6 runs of the HY block |
| Acceptor (yeast ABSENT) | `a{11866,11867,11868}_human_90min_hrMS2.raw` (3) | H: the same human digest, water instead of yeast | last 3 runs of the H block, i.e. acquired immediately before the first donor |
| SDRF (D12-truth) | `PXD014415_community_annotated.sdrf.tsv` -> `data/truth/PXD014415/` | H vs HY annotation | |

**Design and source.** PRIDE sample protocol: "the human sample was split into 2 aliquots of 30ug each. 3ug of yeast peptides was added to one
of the human samples, while an equal volume of HPLC grade water was added to the other sample." PRIDE data protocol and the paper (Lim, Paulo,
Gygi, J Proteome Res 2019, doi:10.1021/acs.jproteome.9b00492, PMC7346880, Methods and Figure 1): "the 20 analyses of the human only sample (H)
were queued first followed by the 20 analyses of the human + 10% yeast spike-in sample (HY)" on the same column, to prevent carry-over. File
numbering matches: `a11849..a11868` = H (20), `a11869..a11888` = HY (20). **Yeast is absent from the 3 acceptor runs**; any yeast-specific
peptide quantity transferred into them is a false transfer. Because H was run before HY, the acceptors cannot contain yeast carry-over.
Human abundance per unit total protein differs slightly: log2(HY/H) for human = log2(30/33) = -0.138 (yeast is ~9.1% of HY total).

**Instrument / method.** Orbitrap Fusion Lumos + EASY-nLC 1200, 35 cm Accucore C18 column, about 1 µg on column, 90 min gradient 4-35% ACN,
DDA Top Speed 3 s; MS1 Orbitrap 120k (350-1,400 Th), MS2 Orbitrap 7.5k HCD 30%, 1.6 Th isolation, dynamic exclusion 120 s. Same method for
all 40 runs.

**Why these runs.** Acceptors right next to the donors in time give the smallest RT drift between H and HY, which is the realistic (hard) case
for MBR: transfers are easiest to make there. 6 donors x 3 acceptors fits the brief and is 2.8 GB. More acceptors (any of `a11849..a11865`,
0.30 GB each) can be added by editing the manifest, for example to study drift with run distance.

**Caveat for scoring.** Some tryptic peptides are identical in human and yeast (conserved proteins; also check I/L equivalence). A transfer of
such a peptide into an H run is not false. A false-transfer count must use yeast peptides absent from the human (and contaminant) part of the
search FASTA.

**PRIDE verification.** Project PXD014415 (API v3): title "Evaluating False Transfer Rates from the Match-Between-Runs algorithm with a
two-proteome model", Orbitrap Fusion Lumos, human + S. cerevisiae, CC0, submitted 2019-06-27, public 2019-09-25. All files carry a PRIDE
SHA-1 (except the SDRF); each D12 raw was checked against it on arrival and again after download:
**9/9 raw SHA-1 match, 0 mismatch** (SDRF has no PRIDE SHA-1), 10/10 sizes match PRIDE, sha256 in `d12.sha256`. The submitter MaxQuant outputs
(`H_v_HY_MBR_*_LFQ_*.zip`, 1.6 GB each) were not downloaded.

**Overlap and leakage.** PXD014415 is not in `dev_sealed`, `dev_sealed_extra`, `tier_s`, `t4_selection`, `train`, `e9_d8`, `g2null`, `d9` or
`predictor_training_accessions.tsv` (grep, 2026-10-27); it is in the T4 crawl frame but not drawn. It is the external MBR set named in the D9
section and in lit/02. D12 makes the D9 MBR control replaceable if T6 claims the DDA-astral-entrapment runs.

**FASTA.** The ProteoBench HYE FASTA (see D11) contains human + yeast SwissProt, so it is sufficient; its E. coli entries act as an
entrapment proteome (no E. coli in any D12 run).

## DN (dev, added 2026-10-27 for M3e / Amendment A1.3): wrong-species G2 null sources

Purpose: null spectra for G2 / tail-G2 / G1 / G3 / G5 (DESIGN §6.4, §6.0 A1.3). They are never searched against their own proteome, are never used for training or selection, and appear only as null-provenance spectra inside the G2 mixtures (scripts/m3e/g2mix.py).
Selection (`scripts/m3e/select_nulls.py`, seed 11, fixed before fetching; log `data/manifests/g2null_select.log.json`):
* PRIDE frame `t4_frame.tsv`; single organism from an allow-list of prokaryotes and unicellular eukaryotes far from human; DDA, LFQ, no enrichment, Thermo Q Exactive/Exploris (OT-HCD).
* Hard-excluded: the T4 selection, all dev/sealed/tier-S accessions, the 13 training projects, predictor-training accessions, and the E9/D8 stream. Organisms already used (P. furiosus D1, S. cerevisiae D7), training species (human, mouse, E. coli) and Arabidopsis (arm C) are also excluded.
* 4 projects × 3 runs; files 0.3–1.5 GB.
* Together with D1 (Pfu) and D7 (yeast, 6 runs) this gives 6 null projects (A1.3 requires ≥ 3).
**Rule for T4 redraws:** these 4 projects are dev data and must be excluded from any future T4 sample. They are non-human, so they cannot enter the human E9/D8 stream.
Fetch: `python3 scripts/m3e/fetch_nulls.py`. RAWs go to `data/raw/<acc>/` (manifest `data/manifests/g2null.tsv`, sha256 `data/manifests/g2null.sha256`). Homolog-filter proteomes go to `data/fasta/dn_<upid>.fasta` (sha256 `data/manifests/g2null_fasta.sha256`; UniProt REST, 2026-10-27).
Correction made before any DN file was fetched: the free-text proteome lookup returned *Synechococcus* WH7803 for PXD021787. The project strain is PCC 7002, so the lookup was replaced by UP000001688. Its UniProtKB proteome query is empty, so the sequences come from UniParc `upid:UP000001688` (3,191). PXD023811's strain is not specified, so USA300 ∪ NCTC 8325 is used.
Licences: per PRIDE project metadata (not re-checked per project; research use only).

| id | role | accession | organism | file | size (B) | sha256 (g2null.sha256) |
|---|---|---|---|---|---|---|
| DN00 | dev | PXD015648 | T. brucei | `20170221_QEP2_FBU_FB_collabJanzen_2T1p99xC4x48_2.raw` | 689185737 | `03d9fd7c84950804` |
| DN01 | dev | PXD015648 | T. brucei | `20170221_QEP2_FBU_FB_collabJanzen_2T1p99xC4x72_4.raw` | 791651730 | `5ee89a3d2c3190d8` |
| DN02 | dev | PXD015648 | T. brucei | `20170221_QEP2_FBU_FB_collabJanzen_2T1p99xC9x72_3.raw` | 734444296 | `d6da5273938ae117` |
| DN03 | dev | PXD023811 | S. aureus | `374-A.raw` | 1412798221 | `356cac1ed7953cce` |
| DN04 | dev | PXD023811 | S. aureus | `536-A.raw` | 1295044249 | `d6c0d0343a4e643c` |
| DN05 | dev | PXD023811 | S. aureus | `860-A2.raw` | 1415610074 | `5280b83299b51758` |
| DN06 | dev | PXD021787 | Synechococcus sp. PCC 7002 | `TL190218_KDahlgren_CameronLab_0p5ug_GFP-APEX-B.raw` | 457252077 | `912d0c5dc86fd05c` |
| DN07 | dev | PXD021787 | Synechococcus sp. PCC 7002 | `TL190621_KDahlgren_CameronLab_250ng_GFP-APEX-A.raw` | 1023430429 | `08d97c268f98d984` |
| DN08 | dev | PXD021787 | Synechococcus sp. PCC 7002 | `TL190621_KDahlgren_CameronLab_250ng_psbU-APEX-B.raw` | 772146199 | `9f0bfbc03c2dc7c7` |
| DN09 | dev | PXD055392 | B. cereus ATCC 14579 | `Q15827_MS23-053_6.raw` | 1113700662 | `fa67d5435b20e477` |
| DN10 | dev | PXD055392 | B. cereus ATCC 14579 | `Q15836_MS23-053_11.raw` | 1131961908 | `8dfaeafc194c3235` |
| DN11 | dev | PXD055392 | B. cereus ATCC 14579 | `Q15837_MS23-053_12.raw` | 1034421410 | `dcbe61743907d0eb` |

## D8 (dev, added 2026-10-27 for E9 / tree 3d): repository-scale human DDA stream

The dataset id is D8; it is unrelated to decision D8 in RESEARCH_TREE.md. Manifest `data/manifests/e9_d8.tsv`; selection log
`e9_d8.log.json`; sha256 of every downloaded RAW in `e9_d8.sha256` (path `e9/<id>__<file>.raw`, appended by the streamer as each file
lands). Selector `python3 scripts/e9/select_e9.py` (stdlib, metadata only; seed 11). Streamer `python3 scripts/e9/stream.py --start A --end B`
downloads to `/tmp/vulcan_e9/raw`, searches, scores, and deletes. Only PSM/feature/score tables are kept, on `/tmp`; they are not committed.
Licences follow each PRIDE project's metadata (not re-checked per project; research use only).

- **Frame.** The PRIDE crawl frame `t4_frame.tsv`, restricted to eligible, *Homo sapiens*, DDA, label-free, no enrichment, Thermo,
  single-instrument Q Exactive / Q Exactive Plus / HF / HF-X / Exploris 240/480 projects. That is the OT-HCD instrument class of the default
  scorer (D7). Titles and keywords suggesting HLA/XL/glyco/SILAC/DIA/phospho/TMT/single-cell/AP-MS/targeted were dropped. 1,666 candidate
  projects.
- **Hard exclusions (project level).** Every accession in `t4_selection.tsv` (sealed T4 sample), `tier_s.tsv`, `dev_sealed.tsv`,
  `dev_sealed_extra.tsv` (all D*/T*), `train.tsv` (the 13 M3 training projects) and `predictor_training_accessions.tsv`. That left 1,654
  eligible projects. A re-check after selection found 0 overlaps with each manifest and with the dev/sealed accession list above.
- **Files.** PRIDE category RAW, `.raw`, 0.1-2.0 GB, with names not suggesting QC/blank/wash/standard/DIA/phospho/TMT. At most 8 per
  project, drawn with `random.Random(11 + accession digits)`. 24 projects in the seeded order had no usable file and were skipped (logged).
- **Order (nested N).** Projects are shuffled with seed 11. Runs are ordered round-robin over projects (the first run of every project,
  then the second, and so on). N = 1, 10, 100, 1000 are therefore nested prefixes: N=10 spans 10 projects, N=100 spans 100, and N=1000
  spans all 142.
- **Size.** 1,000 runs, 1.08 TB (N=100: 111 GB). Instruments: QE 509, QE HF 330, QE Plus 78, Exploris 480 59, HF-X 16, Exploris 240 8.
  Submission years 2014-2026.
- **Integrity.** Each file is size-checked and checked against its PRIDE SHA-1. A mismatch is deleted and re-downloaded from scratch; 7% of
  the resumed proxy transfers were corrupt. 350 of the 1,000 files have no PRIDE SHA-1 and cannot be verified. sha256 of the final copy
  is in `e9_d8.sha256`.
- **Completed.** 982 of 1,000 runs, from 140 projects. PXD074511 (8 runs) cannot be opened by RawFileReader. PXD072218 (8 runs) gives
  0 PSMs. E9_0052 and E9_0735 fail to decode. See results/2026-10-27_E9_scale/NOTES.md.
- **Role.** dev: development curves only, never claims. A1.1 (DESIGN §6.0) also uses its first 10 runs by manifest order as the
  not-used-for-selection confirmation set for the frozen default scorer.

## D7 extension and D10cf (dev, added 2026-10-27 for M5 LiP / CF-MS, E10)

**D7 extension (manifest `data/manifests/d7b.tsv`, sha256 `d7b.sha256`, fetched with `fetch_tier_s.py --manifest d7b --roles dev`).** Six more
PXD031616 DDA runs, so that D7 has 3 replicates per digest and condition (MSstatsLiP needs replicates): LiP F1 rep 2/3 `LM2494`, `LM2497`; TrP F1
rep 2/3 `LM2489`, `LM2502`; LiP M1 rep 3 `LM2498`; TrP M1 rep 3 `LM2500` (`_DDA.raw`, 1.30-1.49 GB each, 8.3 GB). Same project, licence (CC0), design
mapping (submitter MSStats reports) and caveats as D7 above. Design table: `scripts/m5lc/common.py:D7_DESIGN`.

**D10cf (manifest `data/manifests/d10cf.tsv`, sha256 `d10cf.sha256`, fetched with `fetch_tier_s.py --manifest d10cf --roles dev`).** DDA CF-MS dev
set for E10. The id is D10cf (not D10) to avoid confusion with decision D10.
- **Data.** PXD015406 (Mallam et al. 2019 Cell Rep "Systematic discovery of endogenous human ribonucleoprotein complexes", RNA DIF-FRAC, Marcotte
  lab), HEK293T native lysate, size-exclusion chromatography. Control (no RNase) series `Control_HEK_SEC_020118_<n>{b,c}.raw`, fractions
  **10-38** (29 consecutive fractions, 29.1 GB, 0.72-1.23 GB each; the series has fractions 10-70, 58 GB). Orbitrap Fusion, DDA: **Orbitrap MS1,
  ion-trap CID MS2** (verified from the filter strings). Licence CC0 (PRIDE metadata). Fractions 10-14 are void (no peptide IDs).
- **Why this set.** PRIDE API keyword search (co-fractionation / size exclusion / SEC / complexome / elution profiles / biochemical fractionation,
  2026-10-27). Rejected: Heusel 2019 PXD007038/PXD010288 and PXD036711 (SWATH/DIA), PXD044083 (DIA), PXD059608 (DIA, = D6), PXD033459 and
  PXD001220 (SEC DDA but 1.3-1.9 GB/fraction: > 30 GB for any 30-fraction window; PXD001220 not CC0), PXD066287 (purified proteasome SEC, too few
  proteins for CORUM), Wan 2015 PXD002322 (EBI terms; LTQ-Orbitrap Velos), complexome-profiling BN-PAGE sets (mitochondria only), PXD023211 (TMT).
  PXD015406 is human (CORUM applies), CC0, DDA, a single contiguous SEC series and fits the size budget.
- **Caveats.** IT-CID MS2: no frozen static model for this class (tree 3c), so M5 scores it with hyperscore. Only 29 of 61 fractions, one replicate
  series (RNase series and the PXD014820 Lumos replicate not fetched). No spiked foreign proteome: the E10 null uses paired entrapment peptides.

**M5 reference tables (`data/truth/M5/`, gitignored; `python3 scripts/data/fetch_m5_truth.py`; sha256 `data/manifests/m5_truth.sha256`).**
`omnipath_corum.tsv` (CORUM complexes through the OmniPath web service, rows whose sources include CORUM; CORUM's own download server was unreachable
from the node; CORUM licence: free for academic use), `yeast_uniprot_features.tsv` (UniProt UP000002311 initiator-Met / signal / transit /
propeptide / chain / peptide features, CC BY 4.0; half-tryptic exclusions), `P37840.fasta` (human alpha-synuclein, the D7 spike).

## DN2 (role = dev-null, Amendment A2 fresh null set; walk 2026-10-27/28, PI decision D16)
- Selector: `scripts/a2/select_dn2.py` (tag `prereg-v0.4.2` + the D16 fix d6d199c), seed 23, run through
  `scripts/a2_ops/run_dn2_walk.py` from a clean worktree at 8fb2330 (results/2026-10-27_A2_dn2/DEVIATIONS.md items 1-20).
  Run 1 (2650cb8) admitted 0 because of the reference-proteome query bug (item 11); this is run 2 (the single D16 re-walk,
  resumed after infrastructure failures with identical decisions).
- Manifest `data/manifests/dn2.tsv` (18 RAWs, 6 projects), checksums `data/manifests/dn2.sha256` (sha256), log
  `data/manifests/dn2_select.log.json`; per-file PRIDE SHA-1 status in `results/2026-10-27_A2_dn2/selection_log.tsv`
  (15 verified OK; 3 files of PXD062500 have no PRIDE checksum: size + sha256 only, identical across 2 downloads).
- **6 projects admitted, below the 10-project minimum; the stopping criteria were not met** (candidates exhausted).
  DN2 is the A2 null set as drawn; step 4 was not run (D16 point 4). Own proteomes: UniProt reference proteomes (upid
  per row). Persistent copies: `data/cache/a2/` (gitignored; manifests `results/2026-10-27_A2_dn2/cache_*.sha256`).
- License: PRIDE public data (EMBL-EBI terms of use; data as deposited by the submitters). No sealed data touched.

| set | role | accession | organism | own proteome | file | size (B) | sha256 | PRIDE SHA-1 |
|---|---|---|---|---|---|---|---|---|
| DN2 | dev-null | PXD006818 | Streptomyces lividans tk24 | UP000028682 | SL203_set1_rep1.raw | 580829183 | 2369e190a7295d7c3bf83e43dff0a95fef9602ed9780feff2fb54a5acb0d1e04 | OK |
| DN2 | dev-null | PXD006818 | Streptomyces lividans tk24 | UP000028682 | SL203_set1_rep6.raw | 623186987 | 1f7dbc2c707644688e2dfcabd9c1fa8f2dcdd5e838d9783ffd097ce3da68a414 | OK |
| DN2 | dev-null | PXD006818 | Streptomyces lividans tk24 | UP000028682 | SL322_set1_rep5.raw | 632534029 | 1584d3e7d4ddc1da374d261cee4a0cb190d7692881f1e3f5079c198b9a94298f | OK |
| DN2 | dev-null | PXD018953 | Caballeronia mineralivorans | UP000035963 | Q06711_MS18-007_Biotite_GB1.raw | 568422658 | 1db14e5367e20983d9d6c6dcac896d5ee46fcfa9b2e0d6d843bc9bc1e02767f3 | OK |
| DN2 | dev-null | PXD018953 | Caballeronia mineralivorans | UP000035963 | Q06717_MS18-007_Biotite_GBB3.raw | 557052061 | 94ca424a7d9698efd820c4969d48846e70f59feea189f31af8ed44543671911e | OK |
| DN2 | dev-null | PXD018953 | Caballeronia mineralivorans | UP000035963 | Q06722_MS18-007_Biotite_MBB4.raw | 610791957 | b2f23a3ee390cab783252ba7c137e8da948a5b1874eebc4b5bb126ca99b33e26 | OK |
| DN2 | dev-null | PXD044078 | Bacillus pumilus | UP000001355 | Q12238_MS21-015_Bp2.raw | 1174258877 | 199c59ad7928d72dfbd9005c861d1c720f974cd0851b24796dd295d973f635c8 | OK |
| DN2 | dev-null | PXD044078 | Bacillus pumilus | UP000001355 | Q12239_MS21-015_Bp3.raw | 1184320832 | 49173346782c8d30e0cd6108f26a9248a303d8d70e95017a65abb309ff935b2f | OK |
| DN2 | dev-null | PXD044078 | Bacillus pumilus | UP000001355 | Q12240_MS21-015_Bp4.raw | 1175302694 | edd809ff0585e22a7ee3146342fef81494fbef85f4c7cfdd16f4ce222a073095 | OK |
| DN2 | dev-null | PXD011591 | Rhodobacter capsulatus | UP000002361 | 1010_1.raw | 1018857264 | 933d77c9bff7197a4362d6884b6bf940b404a159235ce027c6ab5fe869e9d0bd | OK |
| DN2 | dev-null | PXD011591 | Rhodobacter capsulatus | UP000002361 | 101_G.raw | 1152273587 | 8c526889384e24cff8bc6287b28751b2d1cb07c1fe5cbeff0a647b1958393214 | OK |
| DN2 | dev-null | PXD011591 | Rhodobacter capsulatus | UP000002361 | D2_G.raw | 973523843 | 5c5bba4f387c490bce6e67620175ffeb3edcf947b7f83186a2165b83bb492339 | OK |
| DN2 | dev-null | PXD062500 | Deinococcus radiodurans | UP000002524 | HFX_BIT_LF_DR1_20231107.raw | 1305413384 | 66aec183bd6d40a5346365f4600bbe67b0d009d2a0e452d70b79f99d0ccd94e6 | no_pride_checksum |
| DN2 | dev-null | PXD062500 | Deinococcus radiodurans | UP000002524 | HFX_BIT_LF_DR2_20231107.raw | 1325258508 | 71dd199ede378b4c0f4a8e36b079c9bc489438323ecba238551ada3aa2fc4216 | no_pride_checksum |
| DN2 | dev-null | PXD062500 | Deinococcus radiodurans | UP000002524 | HFX_BIT_LF_pprI4_20231107.raw | 1411283940 | e2d65b8fbe6f112e9c90dffac9df5fba64373bb732ea7ba53244a11dd3d25e7d | no_pride_checksum |
| DN2 | dev-null | PXD038243 | Emiliania huxleyi | UP000013827 | RCC911_23_SubOpt_1.raw | 961473170 | 7f47b02d7aad12198772912e60c87b4462e660f461836cb7ef208d2d11910be2 | OK |
| DN2 | dev-null | PXD038243 | Emiliania huxleyi | UP000013827 | RCC911_28_SupraOpt_1.raw | 968316736 | 84534444b5522fc2a35a79a95f48db01c2e64951351c46506433c5f2ab187432 | OK |
| DN2 | dev-null | PXD038243 | Emiliania huxleyi | UP000013827 | RCC911_28_SupraOpt_2.raw | 967138500 | 9c1a6f086ba3a7e6f85e19bcfb5a97ba21bdb1e71bd8ce40bf22395567923b24 | OK |

## IT-CID training data (role = `train`; 3c IT-CID static-first scorer, added 2026-10-27)

Selected by `scripts/m3c/select_train_itcid.py` (seed 11, rule committed before any download) from the eligible T4 frame, human/mouse only, with hard exclusion of the union of every `data/manifests/*.tsv` accession except the frame (dev, sealed, tier S, T4 selection, M3 `train.tsv`, predictor-training incl. Prosit 2020 CID/TMT, DN, D8, D7b, D9-D12, D10cf) and a mechanical DN2 guard (no accession seen by the A2 DN2 walk, no ALLOW-reachable organism; 0/58 overlap). Strata: TMT (Fusion/Lumos/Eclipse, SPS-MS3, TMT 6/10/11) 6 projects; label-free (LTQ Orbitrap XL/Velos/Elite, Fusion/Lumos) 7 projects; one run per project. Acceptance after ingest (pre-specified): >= 95% MS2 ion-trap CID; TMT: MS3 >= 50% of MS2 and no TMTpro reporters; LF: no MS3. 20 walked candidates were rejected (17 with Orbitrap MS2, incl. RTS methods; 3 TMTpro by file name) and 2 passing candidates (TMT14, TMT15) were beyond the stratum target — `train_itcid_fetch.log.json`. Fetched by `scripts/m3c/fetch_train_itcid.py` into `data/raw/train_itcid/<acc>/`; sha256 `data/manifests/train_itcid.sha256`; fold map (3 project-disjoint folds, round-robin within stratum) `data/manifests/train_itcid_folds.tsv`.
**Rule for later T4 redraws:** these 13 projects are training data and must be excluded from any future T4 sample.

| id | role | accession | species | stratum | instrument | fold | file | sha256 (prefix) |
|---|---|---|---|---|---|---|---|---|
| TI00 | train | PXD006293 | human | TMT | Orbitrap Fusion | 0 | `m06973_PMQ_02.raw` | `94ec26cb0682c899` |
| TI01 | train | PXD042759 | human | TMT | Orbitrap Fusion Lumos | 1 | `270120-Urine-Set31.raw` | `641cf84306675cea` |
| TI02 | train | PXD009023 | human | TMT | Orbitrap Fusion | 2 | `Fraction_4.raw` | `f0a76de7ff513364` |
| TI03 | train | PXD043848 | mouse | TMT | Orbitrap Fusion | 0 | `m21032.raw` | `257d09f10c22ccea` |
| TI04 | train | PXD024175 | human | TMT | Orbitrap Fusion | 1 | `Set3_6.raw` | `4ebd7d3ecb8932f0` |
| TI05 | train | PXD062779 | human | TMT | Orbitrap Fusion Lumos | 2 | `Req111-SM-Nik-Zainal_A_24.raw` | `3419dde3e333c110` |
| TI06 | train | PXD013226 | mouse | LF | LTQ Orbitrap Elite | 0 | `141210_Christiano_plantaris_3hit_tech1_141216070912.raw` | `31ec30a5e2a32bb5` |
| TI07 | train | PXD005480 | human | LF | LTQ Orbitrap | 1 | `MH_2.raw` | `102e5db6fbc8962c` |
| TI08 | train | PXD044866 | mouse | LF | LTQ Orbitrap Velos | 2 | `20210527_3d_mito_totalprot1.raw` | `542495179e9237f5` |
| TI09 | train | PXD000383 | human | LF | LTQ Orbitrap Velos | 0 | `110804_110530_ProstateLMW_N_226_38.raw` | `4f8236686cd9ac33` |
| TI10 | train | PXD038703 | mouse | LF | LTQ Orbitrap Velos | 1 | `62598.raw` | `288cb01f11f2cda4` |
| TI11 | train | PXD002838 | mouse | LF | LTQ Orbitrap Velos | 2 | `ath017934.raw` | `86c3dbb011bb9d5d` |
| TI12 | train | PXD007725 | human | LF | LTQ Orbitrap Elite | 0 | `20150114_ZhongshanHospital_Af_1_240min_2ug.raw` | `fc410a04032d62df` |

## DN3c (dev, added 2026-10-27 for 3c): IT-CID wrong-species G2 null sources

Yeast (*S. cerevisiae*) ion-trap CID runs used ONLY as null-provenance spectra inside IT-CID G2-original/G5 mixtures (hosts = D10cf fractions); never trained on or used for selection. Selected by `scripts/m3c/select_nulls_itcid.py` (seed 23; same exclusions + DN2 guard; S. cerevisiae is excluded from DN2 by `config/a2.yaml`). Acceptance: both runs of a project >= 95% IT-CID MS2, no MS3. Pre-specified 3 projects; stopped at 2 because the next eligible project (PXD002607) stalled on PRIDE FTP (deviation logged in `null_itcid_fetch.log.json` before any null search). sha256 `data/manifests/null_itcid.sha256`.

| id | role | accession | file | sha256 (prefix) |
|---|---|---|---|---|
| DN3c_00 | dev (DN3c) | PXD008785 | `b023p018_OE3_TR01.raw` | `2f51388e47870acc` |
| DN3c_01 | dev (DN3c) | PXD008785 | `b023p018_MSS2_TR01.raw` | `b40221c1df5e9d7d` |
| DN3c_02 | dev (DN3c) | PXD010849 | `150225_S_ACVM_11_01.raw` | `4e5855cf3689ae40` |
| DN3c_03 | dev (DN3c) | PXD010849 | `150225_S_ACVM_12_01.raw` | `75597631b5864821` |

## DNIT (dev, added 2026-10-27, 3c follow-up): IT-CID fresh null attempt (A2 adapted to IT-CID) — power insufficient

Rules `results/2026-10-27_3c_itcid/followup/NULL_RULES.md` + `config/a2_itcid.yaml` (committed 2814f06 before any listing; dated notes 1–4). Seeded walk (seed 47) over free-living non-host-associated organisms on IT-CID-class instruments, excluding every manifest accession, every DN2-seen or DN2-reachable project and every ALLOW-list genus (3c DN2 guard). Downloads verified by PRIDE SHA-1 (API `checksum` or the project's `checksum.txt`) and size; RAW kept on node NVMe (`/tmp/vulcan_dnit/raw`, re-fetchable from PRIDE). 3 projects admitted, then stopped by the futility rule (|A| = 581, projection 3,099 vs 120,000 needed): **no gate run**. Used for nothing else; spent under P4. sha256 `data/manifests/dnit.sha256`.

| id | role | accession | organism | file | size | sha256 (prefix) |
|---|---|---|---|---|---|---|
| DNIT_00 | dev (DNIT) | PXD014351 | *Dinoroseobacter shibae* | `K3597.raw` | 592121722 | `34dcce036a17c6cc` |
| DNIT_01 | dev (DNIT) | PXD014351 | *Dinoroseobacter shibae* | `K3869.raw` | 509206920 | `340a743145fbe530` |
| DNIT_02 | dev (DNIT) | PXD014351 | *Dinoroseobacter shibae* | `K3876.raw` | 529181278 | `da418eb0b2f5f275` |
| DNIT_03 | dev (DNIT) | PXD008499 | *Nannochloropsis gaditana* | `FcoJavFernandez_sol29-16_M1_01092016.raw` | 1076739282 | `4b3877f4c7958ae2` |
| DNIT_04 | dev (DNIT) | PXD008499 | *Nannochloropsis gaditana* | `FcoJavFernandez_sol29-16_M2_01092016.raw` | 1036648828 | `c519d906db06d25b` |
| DNIT_05 | dev (DNIT) | PXD008499 | *Nannochloropsis gaditana* | `FcoJavFernandez_sol29-16_MB_08072016_160901222432.raw` | 938872964 | `9a006e6933e3a12d` |
| DNIT_06 | dev (DNIT) | PXD011268 | *Trichococcus patagoniensis* | `NS2_30oC_1b.raw` | 287767502 | `c1f01232c0a2e117` |
| DNIT_07 | dev (DNIT) | PXD011268 | *Trichococcus patagoniensis* | `NS3_30oC_1c.raw` | 237275074 | `82be503a9e394272` |
| DNIT_08 | dev (DNIT) | PXD011268 | *Trichococcus patagoniensis* | `NS7_30oC_2c.raw` | 225898538 | `96e7db5a3a351781` |
