# T4 selector dry run (metadata only), 2026-10-27

No raw file downloaded. Script: `scripts/data/select_t4.py` (seed 11). Cache: `data/cache/pride_meta/`.

## Frame

- PRIDE projects crawled: 41754 unique accessions (API header total_records: 41754; crawl 2026-10-26).
- **Metadata-eligible frame: 14031 projects** (`results/2026-10-27_T1driver_dryrun/t4_selector/draw/t4_frame.tsv`). Raw-file presence is verified only for selection candidates: 44 candidates checked, failures/skips {'fewer_than_3_usable_raw_files': 4}.
- Exclusion mode `draw`. Frame without the M7 exclusions (train, dev_manifest:*, e9_d8_stream): 14190 projects; with them: 14031. In `draw` mode allocation and within-cell shuffles use the former and M7-excluded candidates are skipped while drawing; in `frame` mode they use the latter.
- Exclusion reasons (a project can have several; counts per reason over all 41754 projects):

  - no_quant_or_phospho_evidence: 16379
  - dia: 6322
  - other_vendor: 4086
  - no_reference_proteome: 4084
  - topdown_maldi_meta: 2998
  - silac: 2897
  - organism_count_2: 2732
  - xl: 1707
  - ion_trap_only: 495
  - organism_count_3: 444
  - dda_dia_mixed: 375
  - mixed_vendor: 243
  - organism_count_4: 145
  - e9_d8_stream: 142
  - predictor_training: 74
  - organism_count_5: 59
  - organism_count_6: 35
  - organism_count_7: 33
  - organism_count_8: 28
  - dev_sealed: 16
  - train: 13
  - organism_count_9: 12
  - organism_count_11: 8
  - organism_count_13: 7
  - dev_manifest:tier_s: 6
  - organism_count_10: 6
  - dev_manifest:dev_sealed: 6
  - organism_count_24: 5
  - organism_count_0: 5
  - organism_count_26: 5
  - dev_manifest:g2null: 4
  - organism_count_12: 4
  - organism_count_14: 3
  - organism_count_15: 3
  - organism_count_16: 3
  - no_instrument: 3
  - organism_count_20: 2
  - organism_count_28: 2
  - dev_manifest:dev_sealed_extra: 2
  - organism_count_30: 2
  - organism_count_47: 1
  - organism_count_50: 1
  - organism_count_43: 1
  - organism_count_56: 1
  - organism_count_18: 1
  - organism_count_44: 1
  - organism_count_27: 1
  - organism_count_54: 1
  - dev_manifest:d4_ref: 1
  - organism_count_21: 1
  - organism_count_32: 1
  - organism_count_35: 1

### Frame composition

- year_bin: {'2015-2017': 1451, '2018-2020': 3572, '2021-2023': 5199, '<=2014': 267, '>=2024': 3542}
- family: {'astral': 22, 'bruker_qtof': 162, 'orbitrap': 13320, 'tims': 527}
- label: {'LFQ': 8786, 'TMT': 4869, 'UNK': 376}
- enrichment: {'none': 12461, 'phospho': 1570}
- acq_confidence: {'assumed': 1490, 'text': 12541}
- label_confidence: {'': 376, 'structured': 9069, 'text': 4586}
- cells (year x family x label x enrichment) non-empty: 53

## Selection

- Projects: 40 (target 40, allowed 30-50); runs: 300 (cap 300); bytes: 383.7 GB = 0.384 TB (cap 1 TB); files dropped for byte cap: 0.
- Cells hit: 26. Allocation rule: sqrt-weighted systematic PPS (see script header).
- selected by year_bin: {'2015-2017': 6, '2018-2020': 9, '2021-2023': 12, '<=2014': 3, '>=2024': 10}
- selected by family: {'bruker_qtof': 3, 'orbitrap': 34, 'tims': 3}
- selected by label: {'LFQ': 21, 'TMT': 15, 'UNK': 4}
- selected by enrichment: {'none': 28, 'phospho': 12}
- selected by vendor: {'bruker': 6, 'thermo': 34}
- selected acq_confidence: {'assumed': 4, 'text': 36} (non-'text' rows rely on defaults; 'assumed' = no DIA or DDA wording found)

### Estimated ingest time

Assumes 36630 spectra (MS2) per GB raw in total (40k/GB Thermo, 25k/GB Bruker; planning figure) -> 1.41e+07 spectra.
- at 2000 spectra/s: 2.0 h (single process)
- at 5000 spectra/s: 0.8 h (single process)

## DDA vs DIA classification and error modes

PRIDE has no machine-readable acquisition-mode field. Rules (in order): (1) strong DIA = DIA/SWATH/diaPASEF/
PRM/SRM terms or DIA software (DIA-NN, Spectronaut, OpenSWATH, EncyclopeDIA) in title, keywords, tags, experiment types, quantification
methods or software -> dropped; (2) DIA wording only in protocol text -> dropped unless DDA wording is also present, in which
case labelled `mixed` and dropped; (3) DDA wording (data-dependent, DDA, TopN, shotgun) -> DDA/text; (4) none of these -> DDA/assumed.
Error modes: (a) DDA-labelled projects that are DDA library + DIA quantification with the DIA only described in files (kept as
`assumed`/`text` DDA; file-name filter removes DIA/SWATH names only); (b) projects mentioning DIA as a comparison or in the
Introduction-style description (dropped = false negatives of the frame, harmless for validity but shrink it); (c) software lists
rarely record DIA tools so structured recall is low; (d) 'shotgun proteomics' experiment type is applied to DIA too, so it is not used as
DDA evidence; (e) projects with both DDA and DIA runs pass only if no DIA wording exists, hence residual DIA contamination possible
in `assumed` projects (quantified above via acq_confidence). Label errors: TMT/LFQ mentions in protocols can describe a comparison or
earlier work (TMT needs >=2 text mentions or one structured one; SILAC needs >=2); LFQ is often undeclared, so many true LFQ projects are
absent from the frame (`no_quant_or_phospho_evidence`). Phospho needs phospho wording and an enrichment term. Instrument lists are
submitter-chosen, may be incomplete (e.g. Astral listed as Exploris) and may mix vendors (dropped). Organism is matched by species
name to UniProt reference-proteome lineages (strain/synonym mismatches lose a few projects; multi-organism projects dropped).
Year is submissionDate, not acquisition date.

Open uncertainties: spectra per GB; whether the Bruker `.d.zip` ingest path exists; ground-truth DDA/DIA per selected project
should be confirmed from raw-file headers at ingest time (a cheap check: Thermo scan filter / tdf acquisition mode) and a project
failing it replaced from the same cell. Predictor-training exclusions: see `data/manifests/predictor_training_accessions.tsv`
(re-run `--stage frame`/`all` after it changes; ProteomeTools accessions are already excluded).
