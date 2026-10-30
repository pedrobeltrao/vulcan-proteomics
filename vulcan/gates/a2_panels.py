"""Amendment A2 Appendix panels (DESIGN §6.0b), verbatim from reviews/amendment_A2_draft.md.

Member = (panel, gene symbol, NCBI taxid).  Panels:
  H  human exposure, non-keratin (desmosome / cornified envelope; sweat, saliva, tears; blood and serum)  -> triggers
  K  human keratins (listed under Panel H in the Appendix; never trigger A2.3 step 3)                    -> no trigger
  B  animal-derived media (bovine, porcine, chicken)                                                    -> triggers
Resolution to UniProt accessions: scripts/a2/resolve_panels.py (rule in its docstring); the resolved table is
config/a2_panels.tsv and the FASTA sha256 is in config/a2.yaml.
"""
from __future__ import annotations

HUMAN, BOVINE, PIG, CHICKEN = 9606, 9913, 9823, 9031
MNEMONIC = {HUMAN: "HUMAN", BOVINE: "BOVIN", PIG: "PIG", CHICKEN: "CHICK"}

PANEL_H_SKIN = ["DSP", "DSG1", "DSC1", "DSC3", "JUP", "PKP1", "PKP3", "EVPL", "PPL", "TGM1", "TGM3", "CASP14", "FLG",
                "FLG2", "HRNR", "LOR", "IVL", "SPRR1A", "SPRR1B", "SPRR2A", "KRTDAP", "SBSN", "DMKN", "CDSN", "S100A7",
                "S100A8", "S100A9", "ARG1", "SERPINB3", "SERPINB4", "SERPINB12", "CSTA", "CST6", "CALML5", "ALOX12B",
                "KLK5", "KLK7", "LYPD3"]
PANEL_H_SECRETION = ["DCD", "PIP", "SCGB1D2", "SCGB2A2", "LYZ", "LTF", "AZGP1", "ZG16B", "MUCL1", "LCN1", "STATH",
                     "HTN1", "AMY1A", "CST1", "CST4", "SMR3B"]
PANEL_H_BLOOD = ["ALB", "HBA1", "HBB", "HBD", "APOA1", "APOA2", "SERPINA1", "TF", "HP", "A2M", "FGA", "FGB", "FGG",
                 "IGHG1", "IGKC", "CP", "C3"]
PANEL_K = ["KRT1", "KRT2", "KRT5", "KRT6A", "KRT9", "KRT10", "KRT14", "KRT16", "KRT17", "KRT77"]
PANEL_B_BOVINE = ["ALB", "TF", "APOA1", "SERPINA1", "AHSG", "FETUB", "HBA", "HBB", "LALBA", "CSN1S1", "CSN2", "CSN3",
                  "LGB"]
PANEL_B_PIG = ["ALB", "TF", "HBA", "HBB"]
PANEL_B_CHICKEN = ["OVAL"]

PANEL_MEMBERS: list[tuple[str, str, int]] = (
    [("H", s, HUMAN) for s in PANEL_H_SKIN + PANEL_H_SECRETION + PANEL_H_BLOOD]
    + [("K", s, HUMAN) for s in PANEL_K]
    + [("B", s, BOVINE) for s in PANEL_B_BOVINE]
    + [("B", s, PIG) for s in PANEL_B_PIG]
    + [("B", s, CHICKEN) for s in PANEL_B_CHICKEN])

TRIGGER_PANELS = ("H", "B")       # A2.3 step 3: >= 2 detected members of these panels exclude the project


def member_id(panel: str, sym: str, tax: int) -> str:
    return f"{panel}|{sym}|{tax}"


def parse_panel_accession(token: str):
    """FASTA first token 'PNL|<panel>|<symbol>|<taxid>|<acc>' -> member id '<panel>|<symbol>|<taxid>' (None if the
    token is not a panel entry)."""
    if not token.startswith("PNL|"):
        return None
    parts = token.split("|")
    return "|".join(parts[1:4])
