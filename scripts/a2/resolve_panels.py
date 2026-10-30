#!/usr/bin/env python3
"""Amendment A2 (DESIGN §6.0b) Appendix: resolve Panel H / keratins / Panel B to UniProt reviewed entries and fetch
cRAP.  Run ONCE before the tag prereg-v0.4.2 (it reads UniProt and the GPM FTP only; no DN2 / PRIDE access).  The
output FASTA sha256 values are written into config/a2.yaml (panels.*) and must match at every later use
(vulcan.gates.a2_admission.check_panel_files).

Resolution rule (fixed here; A2.11 item 5 "panel resolution"):
  * each panel member is (panel, symbol, NCBI taxid): H = human-exposure non-keratin (9606), K = keratins (9606,
    listed but never trigger A2.3 step 3), B = animal-derived media (bovine 9913, porcine 9823, chicken 9031);
  * query UniProtKB REST: (gene_exact:<symbol>) AND (organism_id:<taxid>) AND (reviewed:true);
    if that returns no entry, fall back to the entry name: (id:<symbol>_<MNEMONIC>) AND (reviewed:true)
    (needed for chicken OVAL, whose gene is SERPINB14);
  * EVERY returned reviewed entry is kept (stricter: more sequences can only remove more own peptides and detect more
    panel material), canonical sequence AND all isoforms (format=fasta, includeIsoform=true);
  * a member that resolves to nothing is a hard error (no silent drop);
  * FASTA header first token = PNL|<panel>|<symbol>|<taxid>|<UniProt accession incl. isoform suffix>; entries are
    written sorted by (panel, symbol, taxid, accession) so the file is order-stable;
  * cRAP = the GPM common Repository of Adventitious Proteins (ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta), first
    header token rewritten to CRAP|<original token>.  A1's contaminant list (contaminants_frankenfield2022.fasta) is
    used IN ADDITION (stricter: the A2 "cRAP" set = GPM cRAP union Frankenfield 2022).

Outputs: data/fasta/a2_panels.fasta, data/fasta/a2_crap_gpm.fasta (gitignored, like every FASTA; sha256 in config),
config/a2_panels.tsv (member -> accessions; committed).

    python3 scripts/a2/resolve_panels.py
"""
from __future__ import annotations

import hashlib
import os
import sys
import time
import urllib.parse
import urllib.request

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
from vulcan.gates.a2_panels import MNEMONIC, PANEL_MEMBERS  # noqa: E402

FASTA_DIR = os.path.join(ROOT, "data", "fasta")
OUT_PANEL = os.path.join(FASTA_DIR, "a2_panels.fasta")
OUT_CRAP = os.path.join(FASTA_DIR, "a2_crap_gpm.fasta")
OUT_TSV = os.path.join(ROOT, "config", "a2_panels.tsv")
CRAP_URL = "ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta"
UP = "https://rest.uniprot.org/uniprotkb/stream?"


def get(url: str) -> str:
    for k in range(5):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "vulcan-a2"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read().decode()
        except Exception:  # noqa: BLE001
            if k == 4:
                raise
            time.sleep(3 * (k + 1))
    raise RuntimeError(url)


def fasta_records(txt: str):
    hdr, cur = None, []
    for line in txt.splitlines():
        if line.startswith(">"):
            if hdr is not None:
                yield hdr, "".join(cur)
            hdr, cur = line[1:].strip(), []
        elif line.strip():
            cur.append(line.strip())
    if hdr is not None:
        yield hdr, "".join(cur)


def uniprot_fasta(query: str) -> list[tuple[str, str]]:
    url = UP + urllib.parse.urlencode({"query": query, "format": "fasta", "includeIsoform": "true"})
    return list(fasta_records(get(url)))


def resolve(panel: str, sym: str, tax: int) -> tuple[str, list[tuple[str, str]]]:
    q = f"(gene_exact:{sym}) AND (organism_id:{tax}) AND (reviewed:true)"
    recs = uniprot_fasta(q)
    if not recs:
        q = f"(id:{sym}_{MNEMONIC[tax]}) AND (reviewed:true)"
        recs = uniprot_fasta(q)
    if not recs:
        raise SystemExit(f"resolve_panels: {panel}:{sym}:{tax} resolved to no reviewed UniProt entry")
    return q, recs


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    os.makedirs(FASTA_DIR, exist_ok=True)
    rows, entries = [], []
    for panel, sym, tax in PANEL_MEMBERS:
        q, recs = resolve(panel, sym, tax)
        accs = []
        for hdr, seq in recs:
            acc = hdr.split()[0].split("|")[1]          # sp|P12345-2|NAME_HUMAN ...
            accs.append(acc)
            entries.append(((panel, sym, tax, acc), f"PNL|{panel}|{sym}|{tax}|{acc} {hdr}", seq))
        rows.append((panel, sym, tax, q, ",".join(sorted(accs))))
        print(panel, sym, tax, len(accs), flush=True)
    entries.sort(key=lambda e: e[0])
    with open(OUT_PANEL, "w") as fh:
        for _, h, s in entries:
            fh.write(f">{h}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    crap = list(fasta_records(get(CRAP_URL)))
    with open(OUT_CRAP, "w") as fh:
        for h, s in crap:
            fh.write(f">CRAP|{h}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    with open(OUT_TSV, "w") as fh:
        fh.write("panel\tsymbol\ttaxid\tquery\taccessions\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print("panel fasta", OUT_PANEL, len(entries), "entries sha256", sha256(OUT_PANEL))
    print("crap fasta", OUT_CRAP, len(crap), "entries sha256", sha256(OUT_CRAP))
    print("retrieved", time.strftime("%Y-%m-%d"))


if __name__ == "__main__":
    main()
