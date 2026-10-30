#!/usr/bin/env python3
"""Amendment A4 (draft): METADATA-ONLY count of the DN3 candidate pool, for sizing (A4.3, A4.7).

Reads only local files: the 2026-10-26 PRIDE crawl table data/cache/pride_meta/all_projects_classified.tsv (project
metadata, 41,754 rows; no file listing), data/manifests/*, config/a2_excluded_accessions.tsv, the DN2 / DNIT / A1
selection logs, and reviews/checks/a4_genus_lineage.json (UniProt taxonomy lineage per genus, fetched from
rest.uniprot.org/taxonomy while drafting -- not PRIDE).  No PRIDE file was listed, nothing was downloaded, no
spectrum was read.  The A4 selector rebuilds the lineage table at the tag and freezes it (A4.3.2); these numbers are
for sizing only.
    conda_envs/pe/bin/python reviews/checks/a4_pool_count.py   # writes reviews/checks/a4_pool_count.json
"""
from __future__ import annotations

import collections
import glob
import json
import os
import re

import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
HERE = os.path.dirname(os.path.abspath(__file__))
A1_INSTR = r"^(Q Exactive( HF| HF-X| Plus)?|Orbitrap Exploris (480|240))$"
EXT_INSTR = r"^(Orbitrap Fusion( Lumos)?|Orbitrap Eclipse|Orbitrap Ascend|Orbitrap Exploris 120)$"
DENY_GENERA = """Homo Mus Escherichia Arabidopsis Saccharomyces Pyrococcus Trypanosoma Staphylococcus Synechococcus Shigella
Salmonella Klebsiella Enterobacter Yersinia Acinetobacter Mycobacterium Mycobacteroides Neisseria Haemophilus
Bordetella Brucella Francisella Legionella Listeria Clostridioides Clostridium Bacteroides Phocaeicola Parabacteroides
Bifidobacterium Lactobacillus Lacticaseibacillus Lactiplantibacillus Limosilactobacillus Apilactobacillus Lactococcus
Leuconostoc Pediococcus Fructobacillus Helicobacter Campylobacter Vibrio Aliivibrio Chlamydia Chlamydophila Borrelia
Borreliella Treponema Mycoplasma Mycoplasmoides Metamycoplasma Ureaplasma Streptococcus Enterococcus Porphyromonas
Prevotella Tannerella Fusobacterium Aggregatibacter Moraxella Pasteurella Actinobacillus Coxiella Rickettsia Ehrlichia
Anaplasma Wolbachia Buchnera Sodalis Photorhabdus Xenorhabdus Serratia Proteus Citrobacter Cronobacter
Stenotrophomonas Achromobacter Burkholderia Aeromonas Cutibacterium Micrococcus Kocuria Schaalia Gardnerella
Akkermansia Faecalibacterium Roseburia Ruminococcus Eubacterium Blautia Veillonella Endozoicomonas Candida
Cryptococcus Malassezia Histoplasma Paracoccidioides Coccidioides Pneumocystis Blastomyces Talaromyces Emergomyces
Madurella Vairimorpha Anncaliia Nosema Encephalitozoon Durusdinium Cladocopium Symbiodinium Breviolum Naegleria
Acanthamoeba Candidatus""".split()
DENY_SPECIES_RE = (r"pseudomonas aeruginosa|aspergillus fumigatus|bacillus (?:cereus|anthracis|thuringiensis|"
                   r"weihenstephanensis|mycoides)|corynebacterium diphtheriae|pyrococcus furiosus")
DENY_LINEAGE = ("Metazoa", "Viruses", "Apicomplexa", "Kinetoplastea", "Metamonada", "Parabasalia", "Entamoeb",
                "Microsporidia")


def org_class(lineage: str) -> str:
    if not lineage or lineage == "ERR":
        return "unresolved"
    if any(k in lineage for k in DENY_LINEAGE):
        return "denied_lineage"
    if "Viridiplantae" in lineage:
        return "II_land_plant" if "Embryophyta" in lineage else "I_microbe"
    return "I_microbe"          # Bacteria, Archaea, Fungi, other eukaryotic microbes


def exclusions() -> set[str]:
    ex = set()
    for m in glob.glob(os.path.join(ROOT, "data/manifests/*.tsv")):
        if m.endswith("t4_frame.tsv"):
            continue
        try:
            d = pd.read_csv(m, sep="\t", dtype=str)
        except Exception:
            continue
        for c in d.columns:
            if c.lower() in ("accession", "project", "pxd"):
                ex |= set(d[c].dropna())
    ex |= set(pd.read_csv(os.path.join(ROOT, "config/a2_excluded_accessions.tsv"), sep="\t", header=None,
                          dtype=str)[0])
    ex |= {x["accession"] for x in json.load(open(os.path.join(ROOT, "data/manifests/dn2_select.log.json")))["log"]}
    ex |= {x["accession"] for x in json.load(open(os.path.join(ROOT, "data/manifests/dnit_select.log.json")))["log"]}
    ex |= {x[0] for x in json.load(open(os.path.join(ROOT, "data/manifests/g2null_select.log.json")))["tried"]}
    return {e for e in ex if str(e).startswith("PXD")}


def main():
    a = pd.read_csv(os.path.join(ROOT, "data/cache/pride_meta/all_projects_classified.tsv"), sep="\t",
                    dtype=str).fillna("")
    n0 = len(a)
    a = a[a.exclude_reasons.isin(["", "no_quant_or_phospho_evidence"])]
    a = a[(a.vendor == "thermo") & (a.acquisition == "DDA") & (a.enrichment == "none") & (a.label.isin(["LFQ", "UNK"]))]
    inst = a.instruments.str.split(";")
    a = a.assign(instr_tier=[("A1" if all(re.match(A1_INSTR, x) for x in s) else
                              "EXT" if all(re.match(A1_INSTR, x) or re.match(EXT_INSTR, x) for x in s) else "")
                             if s != [""] else "" for s in inst])
    a = a[a.instr_tier != ""]
    ex = exclusions()
    n_ex = int(a.accession.isin(ex).sum())
    a = a[~a.accession.isin(ex)]
    lin = json.load(open(os.path.join(HERE, "a4_genus_lineage.json")))
    a = a.assign(genus=a.organisms.str.split(" ").str[0].str.capitalize())
    a = a.assign(cls=[org_class(lin.get(g, "")) for g in a.genus])
    deny = a.genus.isin(DENY_GENERA) | a.organisms.str.lower().str.contains(DENY_SPECIES_RE)
    keep = a[~deny & a.cls.isin(["I_microbe", "II_land_plant"])]
    tiers = {"I (microbes, A1 instruments)": keep[(keep.cls == "I_microbe") & (keep.instr_tier == "A1")],
             "II (land plants, A1 instruments)": keep[(keep.cls == "II_land_plant") & (keep.instr_tier == "A1")],
             "III (microbes + plants, Fusion/Lumos/Eclipse/Ascend/Exploris 120)": keep[keep.instr_tier == "EXT"]}
    out = {"crawl_rows": n0, "after_metadata_filters": int(len(a) + n_ex), "excluded_by_prior_use": n_ex,
           "n_exclusion_accessions": len(ex),
           "class_counts_after_exclusion": a.cls.value_counts().to_dict(), "denied_by_genus_or_species": int(deny.sum()),
           "tiers": {k: {"projects": int(len(v)), "genera": int(v.genus.nunique()), "labs": int(v.lab.nunique()),
                         "label_LFQ": int((v.label == "LFQ").sum()), "label_UNK": int((v.label == "UNK").sum()),
                         "top_genera": collections.Counter(v.genus).most_common(15)} for k, v in tiers.items()}}
    json.dump(out, open(os.path.join(HERE, "a4_pool_count.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "tiers"}))
    for k, v in out["tiers"].items():
        print(k, {kk: vv for kk, vv in v.items() if kk != "top_genera"}, v["top_genera"][:8])


if __name__ == "__main__":
    main()
