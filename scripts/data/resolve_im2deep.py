#!/usr/bin/env python3
"""Resolve the 21 IM2Deep Supp Table 1 accessions whose digit '3' was dropped by PDF text extraction (pattern
'PXD0<gap>NNNN' etc.). Candidate = gap filled with '3'; kept only if PRIDE says the project used a timsTOF and the
submission year is 2022-2024 (paper searched PRIDE for timsTOF Pro/Pro2 data). Appends rows to the manifest."""
import csv, json, os, urllib.request
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
P = os.path.join(ROOT, "data", "manifests", "predictor_training_accessions.tsv")
cands = {"PXD039469": "Zila", "PXD043511": "Ries", "PXD035986": "Romero-Gavilan", "PXD036206": "Pahmeier", "PXD037288": "Lamsal",
         "PXD036127": "Ries", "PXD035987": "Ries", "PXD038824": "Filatov", "PXD037622": "de Jonckheere", "PXD043166": "Wang",
         "PXD035675": "Zhu", "PXD036970": "Kovarik", "PXD036746": "Scally", "PXD037945": "Huang", "PXD039646": "Triglia",
         "PXD037089": "Wang", "PXD043226": "Feng", "PXD038840": "Li", "PXD036191": "Xu", "PXD036482": "Klusch", "PXD043082": "-"}
rows = []
for a, who in cands.items():
    try:
        d = json.load(urllib.request.urlopen(urllib.request.Request("https://www.ebi.ac.uk/pride/ws/archive/v3/projects/" + a), timeout=60))
        inst = [i["name"] for i in d.get("instruments", [])]
        ref = " ".join(r.get("referenceLine", "") for r in d.get("references", []))
        ok = any("tims" in i.lower() for i in inst) and d["submissionDate"][:4] in ("2022", "2023", "2024")
        print(a, who, d["submissionDate"][:10], inst, who.split()[0] in ref, d["title"][:70], "OK" if ok else "NO")
        if ok:
            rows.append(["IM2Deep", "multiconformer (2025)", a, "Devreese 2025 Supp Table 1 (gap-filled accession; %s; verified timsTOF in PRIDE)" % who,
                         "medium", "train", "digit recovered by inference"])
    except Exception as e:
        print(a, "ERR", e)
if rows:
    have = {r["accession"] for r in csv.DictReader(open(P), delimiter="\t")}
    with open(P, "a", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        for r in rows:
            if r[2] not in have:
                w.writerow(r)
