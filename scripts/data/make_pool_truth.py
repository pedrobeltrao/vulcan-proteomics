#!/usr/bin/env python3
"""Split the PXD000138 library FASTA (IPI human + 96 pool entries IPI:A1..H12) into one FASTA per pool.

Pool N (HCD file N.raw) is plate well A1..H12 row-major (N=5 -> A5; confirmed by the submitter phosphoRS exports, whose
protein-group accession for 5.raw is A5).  Each pool entry is the concatenation of the designed library peptides (an
unmodified peptide plus its phosphorylated counterparts are searched as variable Phospho (STY) on the same sequence).
Dev pools go to data/truth/PXD000138/pool_<well>.fasta (committable, ~20 kB); sealed pools to data/truth/sealed/PXD000138/.
NOTE: the concatenation does not mark peptide boundaries (library peptides contain internal K/R by design), so the exact
per-pool peptide list needs the Nat Biotechnol 2013 SI; searching with the pool entry as protein, trypsin/P with up to 4
missed cleavages, as the submitters did, recovers it.
"""
import csv, hashlib, os, textwrap
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
FA = os.path.join(ROOT, "data/truth/sealed/PXD000138/009606_IPI_v3_72_with_P_Lib_062111.fasta")
entries, name = {}, None
for l in open(FA):
    if l.startswith(">"):
        name = l[1:].strip(); entries[name] = []
    elif name: entries[name].append(l.strip())
shas = {}
want = {"dev": [], "sealed": []}
for r in csv.DictReader(open(os.path.join(ROOT, "data/manifests/dev_sealed.tsv")), delimiter="\t"):
    if r["accession"] == "PXD000138" and r["filename"].endswith(".raw"):
        want[r["role"]].append(int(r["filename"][:-4]))
for role, pools in want.items():
    out = os.path.join(ROOT, "data/truth/sealed/PXD000138" if role == "sealed" else "data/truth/PXD000138")
    os.makedirs(out, exist_ok=True)
    for n in pools:
        well = "ABCDEFGH"[(n - 1) // 12] + str((n - 1) % 12 + 1)
        seq = "".join(entries[f"IPI:{well}"])
        txt = f">pool_{well} PXD000138 file {n}.raw library entry IPI:{well}\n" + "\n".join(textwrap.wrap(seq, 60)) + "\n"
        fp = os.path.join(out, f"pool_{well}.fasta")
        with open(fp, "w") as fh: fh.write(txt)
        shas[os.path.relpath(fp, os.path.join(ROOT, "data"))] = hashlib.sha256(txt.encode()).hexdigest()
        print(role, n, well, len(seq))
# derived truth is deterministic given the library FASTA (itself checksummed in dev_sealed.sha256); record its sha256 so verify_data.py can check it
with open(os.path.join(ROOT, "data/manifests/derived_truth.sha256"), "w") as fh:
    for k in sorted(shas): fh.write(f"{shas[k]}  {k}\n")
