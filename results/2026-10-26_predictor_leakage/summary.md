# Predictor training-set leakage check, 2026-10-26

Manifest: `data/manifests/predictor_training_accessions.tsv` (predictor, version, accession, source, confidence, role, note; built by
`scripts/data/build_predictor_manifest.py` + `scripts/data/resolve_im2deep.py`). Sources were papers, supplements, Zenodo/Koina/docs pages and the
PRIDE API (project titles/references). No raw data was fetched.

## Coverage and caveats
- Resolved with accessions: AlphaPeptDeep (paper Data availability list, 17 PXD + HLA sets), DeepLC 2021 (Table of train/val/test sets), MS2PIP docs table
  (HCD2019, iTRAQ, TMT, TTOF5600, Immuno-HCD, CID-TMT, timsTOF2023/2024), IM2Deep multiconformer (Supp Table 1: 30 projects; 21 accessions were garbled
  in the PDF text, the dropped digit was '3', each recovered accession was checked in PRIDE: timsTOF, 2022-2023, reference matches; one, `PXD04?82`, stays
  unresolved), Prosit 2019/2020/TMT/timsTOF/XL models, Prosit-PTM (2025).
- **Not resolvable**: MS2PIP HCD2021 training set ("combined dataset, 520,579 peptides", accessions only in Supp Table S1, which was not retrievable; ProteomeTools
  rows are kept as low-confidence); DeepLC v4 multitask model (6,543 LC setups, not enumerated in the docs); Prosit-PTM new ProteomeTools PTM data (no accession);
  MS2PIP HCD training MassIVE-KB (no PXD); AlphaPeptDeep gives a combined list without per-accession train/test roles (all rated "medium", treated as training).
- Roles matter: Prosit-XL and Prosit-PTM name PXD014337, PXD029252, PXD009449 and PXD007058 as **hold-out or external benchmarks, not training data**.

## Overlaps with our dev/sealed accessions (DESIGN 6.6)
| Our set | Accession | Predictor | Finding | Action |
|---|---|---|---|---|
| D2 / T2 (FLR, sealed) | PXD009449 (ProteomeTools 21 PTM, Zolg 2018) | AlphaPeptDeep (listed training/testing data, PTM transfer-learning); DeepLC 2021 (explicit train/val/test split, "ProteomeTools PTM") | **Overlap** | see below |
| D2 / T2 | PXD000138 (Marx 2013 synthetic phosphopeptide library) | AlphaPeptDeep (listed training/testing data) | **Overlap** | see below |
| D2 | PXD007058 (Ferries) | Prosit-PTM 2025 external test only | benchmark known to the authors, not training | note only |
| D3 / D3b / T3 | PXD014337, PXD029252 | Prosit-XL (Koina Prosit_2023_intensity_XL_CMS2, Prosit_2024_intensity_XL_NMS2) hold-out sets | not training; models were tuned/validated on them, unfair if the metric is the XL predictor gain | report as such; use XL rescoring only with the "known benchmark" caveat |
| D0, D1, D4 (PXD016766), D6 (PXD059608), T7 (PXD049412), T8 (PXD015446), ProteoBench PXD028735 (T1) and PXD046444 (T6), T5 | none found | none | clean | none |
| ProteomeTools PXD004732 / PXD010595 / PXD021013 | all Prosit 2019/2020/TMT/timsTOF models; AlphaPeptDeep; DeepLC (004732, 010595); possibly MS2PIP HCD2021 | training | not in our lists; already excluded from T4 | keep excluded |

**PXD009449 vs Prosit.** It is the ProteomeTools PTM deposit (Zolg et al. 2018), a separate PRIDE project from the three tryptic/non-tryptic parts (PXD004732,
PXD010595, PXD021013). It is not listed among the Prosit 2019/2020 intensity/iRT training data (those models only support Ox/CAM anyway), and Prosit-PTM (2025)
holds it out. So Prosit_2019_*/2020_* models did not train on it as far as published statements go (high confidence), while AlphaPeptDeep and DeepLC did.

## Recommendations
1. FLR sealed set T2 and dev set D2 (PXD009449, PXD000138): **do not use AlphaPeptDeep or DeepLC features (MS2, RT, CCS) in any FLR/localisation claim on these data**. Use Prosit_2019_intensity / Prosit_2020_intensity_HCD
   (neither saw them) or MS2PIP HCD2019 (training MassIVE-KB, no overlap) with the low-confidence caveat on HCD2021, and prefer Prosit-PTM only if its hold-out of PXD009449 is confirmed
   on Koina. For RT use Prosit_2019_irt (clean). AlphaPeptDeep/DeepLC may be used only as a clearly labelled "leaky" ablation.
   Alternative: retrain/fine-tune AlphaPeptDeep without those datasets (design allows "retrained").
2. PXD007058 and the XL sets (PXD014337, PXD029252): usable, with a note that Prosit-PTM/Prosit-XL authors benchmarked on them.
3. T4: 74 PRIDE projects are removed from the frame by predictor-training status (32 of them were otherwise eligible); AlphaPeptDeep's list contains several DIA/phospho projects. Low-confidence flags (PXD001189, PXD008034 evaluation sets of MS2PIP) remain flagged only; none was selected.
   If MS2PIP HCD2021 is used, treat ProteomeTools as trained (already excluded).
4. T5 (yeast, not chosen yet): avoid PXD003472 and PXD008812 (DeepLC yeast / Arabidopsis training), PXD000269, PXD001250, and any ProteomeTools project; verify the lab is not Kuster/Wilhelm/Mann/Heck.
5. TODO before the prereg tag: obtain MS2PIP Supp Table S1 and the DeepLC v4 setup list, then re-run `select_t4.py --stage all` (it reads the manifest); record results in DATA.md.
