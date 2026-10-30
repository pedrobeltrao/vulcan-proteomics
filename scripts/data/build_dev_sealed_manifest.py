#!/usr/bin/env python3
"""Generate data/manifests/dev_sealed.tsv (DESIGN.md section 6.6 roles) with the seeded pool selection recorded here.

Selection (seed 11 = "pilot subsample / project folds" seed of DESIGN 6.1):
  * PXD000138 (Marx 2013, 96 synthetic phospho-libraries; HCD raw files are 1.raw..96.raw): rng=Random(11).sample(1..96, 8);
    first 4 -> dev (D2), last 4 -> sealed (T2).  Pool N is plate well A1..H12 row-major (N=5 -> A5).
  * PXD009449 (ProteomeTools PTM, 25 pools = modification kits, file tag ..._PTM_TrainKit_<Pool>_200fmol_3xHCD_R*.raw):
    the 4 *_Unmod pools are excluded from both roles (they hold the unmodified counterparts of every modified pool of
    the same residue class, so they would leak sequences across the split).  Of the 21 modified pools, Ymod_Phospho is
    pinned to dev (phosphosite localisation is the FLR use case); the other 20 are shuffled with Random(11): next 3 -> dev,
    following 4 -> sealed.  Only the 3xHCD raw file of each pool is fetched (HCD-only, 0.52-0.64 GB; 2xIT_2xHCD, ETD and
    survey files are not).  The remaining 13 pools are unused, not sealed.
Everything else is explicit below.  Re-running this script rewrites the manifest (url/size are re-resolved by the fetcher).
"""
import csv, json, os, random, re, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(ROOT, "data", "manifests", "dev_sealed.tsv")
COLS = ["id", "accession", "filename", "url", "expected_size", "purpose", "role", "set_id", "member"]
rows = []


def add(rid, acc, fn, purpose, role, set_id, url="", member="", size=""):
    rows.append(dict(id=rid, accession=acc, filename=fn, url=url, expected_size=size, purpose=purpose, role=role, set_id=set_id, member=member))


# ---- selection ------------------------------------------------------------------------------------------------
marx = random.Random(11).sample(range(1, 97), 8)
marx_dev, marx_sealed = marx[:4], marx[4:]
mods = sorted(["Kmod_Acetyl", "Kmod_Biotinyl", "Kmod_Butyryl", "Kmod_Crotonyl", "Kmod_Dimethyl", "Kmod_Formyl", "Kmod_Glutaryl",
               "Kmod_Hydroxyisobutyryl", "Kmod_Malonyl", "Kmod_Methyl", "Kmod_Propionyl", "Kmod_Succinyl", "Kmod_Trimethyl",
               "Kmod_Ubiquitinyl", "Pmod_Hydroxyproline", "Rmod_Citrullin", "Rmod_Dimethyl_asymm", "Rmod_Dimethyl_symm",
               "Rmod_Methyl", "Ymod_Nitrotyr"])
random.Random(11).shuffle(mods)
pt_dev, pt_sealed = ["Ymod_Phospho"] + mods[:3], mods[3:7]
# PRIDE file / SEARCH-archive names per pool (PXD009449 naming is irregular)
PT = {  # pool: (3xHCD raw filename, SEARCH zip stem)
    "Ymod_Phospho": ("02330a_BA4_4087_1_PTM_TrainKit_Ymod_Phospho_200fmol_3xHCD_R1.raw", "Ymod_Phospho"),
    "Kmod_Dimethyl": ("02330a_GA2_3990_09_PTM_TrainKit_Kmod_Dimethyl_200fmol_3xHCD_R1.raw", "Kmod_Dimethyl"),
    "Rmod_Citrullin": ("02330a_GE1_3990_05_PTM_TrainKit_Rmod_Citrullin_200fmol_3xHCD_R1.raw", "Rmod_Citrullin"),
    "Kmod_Formyl": ("02330a_GC2_3990_11_PTM_TrainKit_Kmod_Formyl_200fmol_3xHCD_R1.raw", "Kmod_Formyl"),
    "Kmod_Acetyl": ("02330a_GG1_3990_07_PTM_TrainKit_Kmod_Acetyl_200fmol_3xHCD_R1.raw", "Kmod_Acetyl"),
    "Kmod_Succinyl": ("02330a_GF2_3991_03_PTM_TrainKit_Kmod_Succinyl_200fmol_3xHCD_R1.raw", "Kmod_Succinyl"),
    "Kmod_Ubiquitinyl": ("02330a_GD3_3991_09_PTM_TrainKit_Kmod_Ubiquitinyl_200fmol_3xHCD_R1.raw", "Kmod_GlyGly"),
    "Kmod_Crotonyl": ("02330a_GA3_3991_06_PTM_TrainKit_Kmod_Crotonyl_200fmol_3xHCD_R1.raw", "Kmod_Crotonyl"),
}
PX9 = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/06/PXD009449/"
PX0 = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2013/05/PXD000138/"


def marx_rows(pools, role, sid):
    for n in pools:
        well = "ABCDEFGH"[(n - 1) // 12] + str((n - 1) % 12 + 1)
        add(sid, "PXD000138", f"{n}.raw", f"{sid} Marx synthetic phospho-library pool {n} (well {well}), HCD, Orbitrap Velos", role, sid)
        add(sid + "-truth", "PXD000138", f"phosphoRS_HCD_{n}.csv", f"{sid} pool {n} submitter phosphoRS PSMs (reference IDs; library peptide entry = IPI:{well} in the library FASTA)", role, sid)


def pt_rows(pools, role, sid):
    for p in pools:
        raw, zst = PT[p]
        add(sid, "PXD009449", raw, f"{sid} ProteomeTools PTM pool {p}, 3xHCD, Orbitrap Fusion Lumos", role, sid)
        add(sid + "-truth", "PXD009449", f"{p}_evidence.txt", f"{sid} pool {p}: MaxQuant evidence.txt (combined/txt) from submitter SEARCH_{zst}.zip; reference IDs, filter 'Raw file' to the 3xHCD file",
            role, sid, url=PX9 + f"SEARCH_{zst}.zip", member=f"{zst}/combined/txt/evidence.txt")


# ---- dev ------------------------------------------------------------------------------------------------------
marx_rows(marx_dev, "dev", "D2")
pt_rows(pt_dev, "dev", "D2")
add("D3b", "PXD029252", "20210409_QExHFX2_RSLC3_DSSO_mainlib_rep2.raw", "D3b DSSO main library rep2 (QExHFX2, 2021-04-09; T3 is a different file/instrument/day)", "dev", "D3b")
add("D3b-truth", "PXD029252", "20210409_QExHFX2_RSLC3_DSSO_mainlib_rep2_KK_1FDR_linkandCSM_Xi1.7.6.5.mzid", "D3b xiSEARCH 1% FDR link+CSM result (reference)", "dev", "D3b")
for i, tag in zip(range(95, 107), ["A1", "A3", "A5", "A7", "A9", "A11", "B1", "B3", "B5", "B7", "B9", "B11"]):
    add("D4", "PXD016766", f"dq_{i:05d}_11cell_90min_sps_{tag}.raw", f"D4 TMT11 Eclipse, 90 min SPS-MS3, well {tag} (12-file subset, 2.6 GB)", "dev", "D4")
extra = os.path.join(ROOT, "data", "manifests", "dev_sealed_extra.tsv")  # D7 / T5 rows chosen after PRIDE API search (hand-curated)
if os.path.exists(extra):
    rows += list(csv.DictReader(open(extra), delimiter="\t"))

# ---- sealed ---------------------------------------------------------------------------------------------------
PB = "https://proteobench.cubimed.rub.de/raws/DDA/"
for c in "AB":
    for r in "123":
        add("T1", "ProteoBench-DDA-QExactive", f"LFQ_Orbitrap_DDA_Condition_{c}_Sample_Alpha_0{r}.raw", "T1 ProteoBench DDA quant module (QExactive HF-X), do not rename", "sealed", "T1", url=PB + f"LFQ_Orbitrap_DDA_Condition_{c}_Sample_Alpha_0{r}.raw")
add("T1-truth", "ProteoBench-DDA-QExactive", "ProteoBenchFASTA_MixedSpecies_HYE.zip", "T1 HYE search FASTA (human, yeast, E. coli, contaminants); expected condition ratios are defined by the ProteoBench module", "sealed", "T1",
    url="https://proteobench.cubimed.rub.de/fasta/ProteoBenchFASTA_MixedSpecies_HYE.zip")
marx_rows(marx_sealed, "sealed", "T2")
pt_rows(pt_sealed, "sealed", "T2")
add("T2-truth", "PXD000138", "009606_IPI_v3_72_with_P_Lib_062111.fasta", "T2/D2 library FASTA: human IPI + 96 pool entries IPI:A1..H12 (the designed peptide concatenations); contains sealed pools, so kept under truth/sealed", "sealed", "T2")
add("T2-truth", "PXD009449", "sequences.fasta", "T2/D2 ProteomeTools PTM kit peptide-sequence FASTA (Kmod/Pmod/Rmod/Ymod/QC packets; all pools of a residue class share a packet)", "sealed", "T2")
add("T3", "PXD029252", "20210203_QExHFX3_RSLC10__DSSO_mainlib_1_5HEK_ctrl.raw", "T3 DSSO ribosomal XL library in HEK background (formerly Tier-S S4 raw; its truth rows stay in tier_s.tsv under data/truth/sealed/)", "sealed", "T3",
    url="https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/05/PXD029252/20210203_QExHFX3_RSLC10__DSSO_mainlib_1_5HEK_ctrl.raw", size="1350506891")
if __name__ == "__main__":
    if os.path.exists(OUT):  # keep already-resolved url / expected_size
        old = {(r["accession"], r["filename"]): r for r in csv.DictReader(open(OUT), delimiter="\t")}
        for r in rows:
            o = old.get((r["accession"], r["filename"]))
            if o and not r["url"]: r["url"] = o["url"]
            if o and not r["expected_size"] and o["url"] == r["url"]: r["expected_size"] = o["expected_size"]
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, COLS, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    print(f"wrote {OUT}: {len(rows)} rows; Marx dev={marx_dev} sealed={marx_sealed}; PTM dev={pt_dev} sealed={pt_sealed}")
