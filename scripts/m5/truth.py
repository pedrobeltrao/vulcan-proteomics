"""M5 XL truth from the library designs (L4 / evaluation only; never imported by vulcan.xl).

D3  (PXD014337, Beveridge 2020): 95(96) synthetic Cas9 peptides in 12 groups, cross-linked within groups and pooled.
D3b (PXD029252, Matzinger 2022): 141(202 rows) ribosomal-protein peptides in groups, cross-linked within groups.
Group tables: IMP-X-FDR@92b2e29 libraries (peplib_Cas9_Beveridge.xlsx; main_library.xlsx re-fetched for dev from
the public URL into /tmp/vulcan_m5/truth, sha256 identical to the manifest; data/truth/sealed is not read).

A target-target link is TRUE iff both peptides map to library peptides that share a group. It is FALSE if
(a) the peptides map to library peptides of different groups only ("between-group"), or
(b) either peptide is not a library peptide: human background protein (D3b search DB has 500 human proteins that
    are absent from this sample), contaminant (D3), or a non-library peptide of a library protein ("non-library").
A peptide maps to a library peptide if the sequences are equal or one contains the other (I = L). The IMP-X-FDR
rule is group-based (residue positions are not checked); we follow it at every level. At PPI level a protein pair
is TRUE iff at least one library group contains peptides of both proteins (for self links: of that protein).
"""
from __future__ import annotations

import os
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
LIB = {
    "D3": os.path.join(ROOT, "data/truth/IMP-X-FDR@92b2e29/peplib_Cas9_Beveridge.xlsx"),
    "D3b": "/tmp/vulcan_m5/truth/main_library.xlsx",
}


def _norm(s: str) -> str:
    return s.replace("I", "L")


def read_library(which: str) -> list[tuple[str, int]]:
    """[(peptide (I->L), group)] from the xlsx (first column; 'GroupN' rows start a group)."""
    sys.path.insert(0, "/tmp/vulcan_m5/pylib")               # openpyxl (scratch install)
    import openpyxl
    wb = openpyxl.load_workbook(LIB[which], read_only=True)
    ws = wb.worksheets[0]
    out, g = [], None
    for row in ws.iter_rows(values_only=True):
        v = row[0]
        if v is None:
            continue
        v = str(v).strip()
        if v.lower().startswith("group"):
            g = int(v[5:])
            continue
        out.append((_norm(v), g))
    return out


class Truth:
    def __init__(self, which: str, prot_seqs: dict):
        self.which = which
        self.lib = read_library(which)
        self.groups_of_pep: dict[str, frozenset] = {}
        self.prot_seqs = {a: _norm(s) for a, s in prot_seqs.items()}
        # protein -> groups with a library peptide of that protein
        self.prot_groups: dict[str, set] = {}
        for p, g in self.lib:
            for a, s in self.prot_seqs.items():
                if p in s:
                    self.prot_groups.setdefault(a, set()).add(g)

    def groups(self, pep: str) -> frozenset:
        p = _norm(pep)
        if p not in self.groups_of_pep:
            self.groups_of_pep[p] = frozenset(g for q, g in self.lib if q == p or q in p or p in q)
        return self.groups_of_pep[p]

    def link_class(self, pa: str, pb: str) -> str:
        ga, gb = self.groups(pa), self.groups(pb)
        if not ga or not gb:
            return "non_library"
        return "true" if ga & gb else "between_group"

    def ppi_true(self, a: str, b: str) -> bool:
        ga, gb = self.prot_groups.get(a, set()), self.prot_groups.get(b, set())
        return bool(ga & gb)
