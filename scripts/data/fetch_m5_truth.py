#!/usr/bin/env python3
"""M5 LiP/CF reference tables (dev; re-fetchable, gitignored; sha256 in data/manifests/m5_truth.sha256).
  data/truth/M5/omnipath_corum.tsv          CORUM complexes via the OmniPath web service (rows whose sources include
                                            CORUM; components = UniProt ids). mips.helmholtz-muenchen.de (CORUM's own
                                            download) was unreachable from the node on 2026-10-27.
  data/truth/M5/yeast_uniprot_features.tsv  UniProt UP000002311 initiator Met / signal / transit / propeptide / chain /
                                            peptide features (half-tryptic exclusions, DESIGN L5 LiP)
  data/truth/M5/P37840.fasta                human alpha-synuclein (D7 spike)
Content at a URL can change between releases; the committed sha256 identifies the copy used in results/2026-10-27_M5_lip_cf.
    python3 scripts/data/fetch_m5_truth.py [--check]"""
import hashlib, os, sys, urllib.request
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
D = os.path.join(ROOT, "data", "truth", "M5")
URLS = {
    "omnipath_corum.tsv": "https://omnipathdb.org/complexes?databases=CORUM&format=tsv&fields=references",
    "yeast_uniprot_features.tsv": "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000002311&format=tsv"
                                  "&fields=accession,length,ft_init_met,ft_signal,ft_transit,ft_propep,ft_chain,ft_peptide",
    "P37840.fasta": "https://rest.uniprot.org/uniprotkb/P37840.fasta",
}
os.makedirs(D, exist_ok=True)
for name, url in URLS.items():
    p = os.path.join(D, name)
    if not os.path.exists(p) and "--check" not in sys.argv:
        urllib.request.urlretrieve(url, p + ".part")
        os.replace(p + ".part", p)
    if os.path.exists(p):
        print(hashlib.sha256(open(p, "rb").read()).hexdigest() + "  truth/M5/" + name)
