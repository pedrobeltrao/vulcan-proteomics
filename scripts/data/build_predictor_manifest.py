#!/usr/bin/env python3
"""Writes data/manifests/predictor_training_accessions.tsv (hand-curated 2026-10-26 from papers/repos; see 'source').
roles: train = stated training data; train_or_test = listed as 'training and testing data' without per-accession role;
holdout/eval = stated NOT in training; app = application data only.  select_t4.py excludes confidence in {high, medium}
AND role in {train, train_or_test, unknown_role}.  Unresolved IM2Deep multiconformer accessions are filled by
scripts/data/resolve_im2deep.py (digit-gap candidates checked against crawled PRIDE metadata)."""
import csv, os
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
R = []
def add(pred, ver, accs, src, conf, role, note=""):
    for a in accs.split():
        R.append([pred, ver, a, src, conf, role, note])

APD = "Zeng et al. 2022 Nat Commun 13:7238 (doi:10.1038/s41467-022-34904-3), Data availability: 'training and testing data from ProteomeXchange'"
add("AlphaPeptDeep", "pretrained generic MS2/RT/CCS (peptdeep; Koina AlphaPept_*_generic)",
    "PXD010595 PXD004732 PXD021013 PXD019854 PXD019086 PXD004452 PXD014525 PXD017476 PXD019347 PXD021318 PXD026805 PXD026824 PXD029545 PXD000269 PXD001250",
    APD, "medium", "train_or_test", "paper does not map accession->train/test split")
add("AlphaPeptDeep", "pretrained generic MS2/RT/CCS", "PXD009449 PXD000138", APD, "medium", "train_or_test",
    "ProteomeTools 21-PTM synthetic peptides (PXD009449) and Marx 2013 synthetic phosphopeptide library (PXD000138): used for PTM transfer-learning/testing of the MS2 model; OVERLAPS our D2/T2 FLR sets")
add("AlphaPeptDeep", "pretrained generic (HLA models)", "MSV000084172 PXD004894 PXD022950", APD, "medium", "train_or_test", "HLA data; MassIVE id for MSV000084172")
DLC = "Bouwmeester et al. 2021 Nat Methods 18:1363 (bioRxiv 10.1101/2020.03.28.013003) Table of training/validation/test peptides"
add("DeepLC", "2021 paper models (Koina Deeplc_hela_hf = HeLa HF)", "PXD006932", DLC, "high", "train", "HeLa HF-X; Koina Deeplc_hela_hf")
add("DeepLC", "2021 paper models", "PXD010595 PXD004732 PXD009449 PXD000954 PXD013477 PXD005573 PXD010154 PXD003472 PXD008812", DLC, "high", "train",
    "dataset-specific train/validation/test splits; PXD009449 = ProteomeTools PTM; overlaps D2/T2")
add("DeepLC", "2021 paper (application only)", "PXD000561", DLC, "high", "app", "CD8 T-cell open-search application, not training")
add("DeepLC", "v4.x multitask pretrained model (6,543 LC setups)", "UNRESOLVED", "deeplc.readthedocs.io / github compomics/DeepLC: training setups not enumerated in docs", "low", "unknown_role", "must be checked in DeepLC repo/zenodo before use; treat as possibly containing any public RP data")
M23 = "MS2PIP docs 'Prediction models' (ms2pip.readthedocs.io); Declercq et al. 2023 NAR gkad335; Gabriels/Declercq Supp Table S1 (not retrievable)"
add("MS2PIP", "HCD2019", "PXD008034", M23, "high", "eval", "test set; training = MassIVE-KB (no PXD)")
add("MS2PIP", "HCD2021", "PXD008034", M23, "high", "eval", "evaluation set listed in docs; training = 'combined dataset (520,579 peptides, tryptic+chymotryptic)' whose accessions are in Supp Table S1 (not retrievable)")
add("MS2PIP", "HCD2021", "PXD004732 PXD010595 PXD021013", M23, "low", "unknown_role", "POSSIBLE: combined tryptic/chymotryptic HCD training set may include ProteomeTools; unverified")
add("MS2PIP", "iTRAQ / iTRAQphospho", "PXD001189", M23, "high", "eval", "")
add("MS2PIP", "TMT", "PXD009495", M23, "high", "eval", "training = Peng Lab TMT spectral library")
add("MS2PIP", "TTOF5600", "PXD000954 PXD001587", M23, "high", "train_or_test", "")
add("MS2PIP", "Immuno-HCD", "PXD005231 PXD020011", M23, "high", "train_or_test", "")
add("MS2PIP", "CID-TMT", "PXD041002 PXD005890", M23, "high", "train_or_test", "")
add("MS2PIP", "timsTOF2023", "PXD043026 PXD046535 PXD046543", M23, "high", "train_or_test", "")
add("MS2PIP", "timsTOF2024", "PXD043026 PXD046535 PXD046543 PXD038782", M23, "high", "train_or_test", "")
IM = "Devreese et al. 2025 Anal Chem / bioRxiv 10.1101/2025.02.18.638865, Supplementary Table 1 (30 PRIDE projects, 1248 runs)"
add("IM2Deep", "multiconformer (2025) and v1 base (Meier 2021 data)", "PXD019086 PXD042416 PXD042478 PXD040890 PXD049281 PXD042114 PXD040481 PXD040521 PXD048960", IM, "high", "train",
    "accession text fully legible in the PDF; the other 21 accessions are partially garbled in the PDF text and are resolved by scripts/data/resolve_im2deep.py")
add("IM2Deep", "multiconformer (2025)", "PXD046507", IM, "high", "eval", "independent demonstration set, not training")
add("IM2Deep", "multiconformer (2025)", "PXD056927", "same paper, main text", "medium", "eval", "independent dataset, species not in training (comparison with Meier model)")
P19 = "Gessulat et al. 2019 Nat Methods 16:509 (ProteomeTools Part I/II); PRIDE titles PXD004732, PXD010595, PXD010871"
add("Prosit", "Prosit_2019_intensity (Koina)", "PXD004732 PXD010595", P19, "high", "train", "ProteomeTools tryptic synthetic peptides, HCD")
add("Prosit", "Prosit_2019_irt (Koina)", "PXD004732 PXD010595", P19, "high", "train", "iRT model on same ProteomeTools data")
add("Prosit", "Prosit_2019_*", "PXD010871", P19, "medium", "derived", "PRIDE 'Prosit' project holding search results of the above data")
P20 = "Wilhelm et al. 2021 Nat Commun 12:3346 Data availability (PXD021013; PXD004732/PXD010595 incorporated from the earlier model)"
add("Prosit", "Prosit_2020_intensity_HCD / _CID (Koina)", "PXD021013 PXD004732 PXD010595", P20, "high", "train", "~30M spectra: 21M tryptic + 9M non-tryptic/HLA ProteomeTools")
add("Prosit", "Prosit_2020_intensity_TMT / Prosit_2020_irt_TMT", "PXD023119 PXD023120", "Gabriel et al. 2022 Anal Chem (Prosit-TMT); PRIDE titles 'ProteomeTools Part IV/V TMT'", "high", "train", "also reuses unmodified ProteomeTools parts (medium)")
add("Prosit", "Prosit_2020_intensity_TMT / Prosit_2020_irt_TMT", "PXD004732 PXD010595 PXD021013", "same", "medium", "train", "unmodified/non-tryptic peptides probably reused; not verified")
add("Prosit", "Prosit_2023_intensity_timsTOF (Koina)", "PXD019086 PXD043844", "Zenodo 8211811 / bioRxiv 10.1101/2023.07.17.549401 (fine-tuned on timsTOF ProteomeTools; base = Prosit 2020 HCD)", "high", "train", "also inherits Prosit_2020_intensity_HCD training data")
add("Prosit", "Prosit_2023_intensity_timsTOF (Koina)", "PXD004732 PXD010595 PXD021013", "same (base model)", "high", "train", "")
XL = "Kalhor et al. 2025 Nat Commun 16 (Prosit-XL), Methods"
add("Prosit", "Prosit_2023_intensity_XL_CMS2 / Prosit_2024_intensity_XL_NMS2", "PXD012546 PXD019926 PXD017711 PXD011861 PXD017620 PXD016554 PXD017695 PXD014675 PXD008550", XL, "high", "train", "cross-link training sets")
add("Prosit", "Prosit_2023_intensity_XL_CMS2 / Prosit_2024_intensity_XL_NMS2", "PXD029252 PXD014337", XL, "high", "eval", "stated hold-out/external benchmarks, NOT training; same as our D3/D3b/T3 (benchmarks known to the model authors)")
add("Prosit", "Prosit_2023_intensity_XL_CMS2 / Prosit_2024_intensity_XL_NMS2", "PXD017711", XL, "high", "train_and_eval", "also used as external set")
add("Prosit-PTM", "2025 (bioRxiv 10.1101/2025.11.07.687302; not on Koina's classic list)", "PXD009449", "Gabriel et al. 2025, 'previously published 21PTMs dataset (Zolg 2018) reserved as holdout'", "high", "eval", "held out; trained on new ProteomeTools PTM extension (no accession in preprint; HuggingFace Wilhelmlab/prospect-ptms-*)")
add("Prosit-PTM", "2025", "PXD007058", "same, Methods 'external phospho peptides PXD007058'", "high", "eval", "external test only; our D2 set is a benchmark of the model authors")
add("Prosit-PTM", "2025", "UNRESOLVED", "no accession for new ProteomeTools PTM training data in preprint", "low", "train", "ProteomeTools PTM extension; includes unmodified, TMT, dimethyl + modified peptides (3,945 raw files)")
os.makedirs(os.path.join(ROOT, "data", "manifests"), exist_ok=True)
out = os.path.join(ROOT, "data", "manifests", "predictor_training_accessions.tsv")
with open(out, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["predictor", "version", "accession", "source", "confidence", "role", "note"])
    w.writerows(R)
print(len(R), "rows ->", out)
