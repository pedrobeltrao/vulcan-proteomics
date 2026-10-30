"""Half-tryptic calling for LiP-MS (DESIGN §4 L5, LiP bullet), applied to identified peptides only.

For every protein occurrence of a peptide [x, y) (0-based, end-exclusive) in protein P:
  N-end tryptic   x == 0, or P[x-1] in {K,R} and P[x] != 'P'
  C-end tryptic   y == len(P), or P[y-1] in {K,R} and (y == len(P) or P[y] != 'P')
A non-tryptic end that sits on a *biological* terminus is not a protease (PK) cut and is excluded:
  * protein N-terminus after initiator-Met removal (x == 1 and P[0] == 'M'), whether or not UniProt annotates it;
  * mature-chain boundaries from UniProt features: signal peptide end, transit peptide end, propeptide ends,
    chain / peptide starts and ends (`ft_signal`, `ft_transit`, `ft_propep`, `ft_chain`, `ft_peptide`, `ft_init_met`).
Peptide-level call over all occurrences (most tryptic occurrence wins, as for the search class):
  'tryptic'       both ends tryptic in some occurrence
  'bio_terminal'  no tryptic occurrence, but each non-tryptic end is a biological terminus (excluded from LiP cuts)
  'half_tryptic'  exactly one non-tryptic, non-biological end (in the best occurrence) -> one PK cut site
  'nonspecific'   both ends non-tryptic and non-biological (not in the semi-specific search space; should not occur)
Cut site: the P1 residue of the PK-cleaved bond, 1-based: N-side cut -> x (residue P[x-1]); C-side cut -> y.
Sequences are matched I/L-insensitively; the original residues are returned (`orig_seq`) for downstream tools.

`insource_flags`: a half-tryptic peptide co-eluting (|apex RT difference| <= tol) in the same run with an
identified peptide that contains it as a prefix or suffix (same protein locus) is flagged as a possible in-source
fragment of that parent.  A flag, never a filter (DESIGN: "flag in-source fragments").
"""
from __future__ import annotations

import re
from collections import defaultdict

import numpy as np
import pandas as pd

_RANGE = re.compile(r"(SIGNAL|TRANSIT|PROPEP|CHAIN|PEPTIDE|INIT_MET) (\?|<?\d+)\.\.(\?|>?\d+)")


def parse_uniprot_features(tsv_path: str) -> dict:
    """accession -> {'n': set of 0-based N-boundaries, 'c': set of 0-based end-exclusive C-boundaries}."""
    df = pd.read_csv(tsv_path, sep="\t", dtype=str).fillna("")
    out = {}
    for _, r in df.iterrows():
        nb, cb = set(), set()
        for col in r.index[1:]:
            for kind, a, b in _RANGE.findall(r[col]):
                if a.startswith("?") or b.startswith("?"):
                    continue
                s, e = int(a.lstrip("<")), int(b.lstrip(">"))   # 1-based inclusive
                if kind in ("SIGNAL", "TRANSIT", "INIT_MET"):
                    nb.add(e)            # mature protein starts after residue e (0-based start e)
                elif kind in ("PROPEP", "CHAIN", "PEPTIDE"):
                    nb.add(s - 1)
                    cb.add(e)
                    if kind == "PROPEP":
                        nb.add(e)        # the next chain starts after the propeptide
                        cb.add(s - 1)
        out[r["Entry"]] = {"n": nb, "c": cb}
    return out


class ProteinDB:
    def __init__(self, accessions, sequences, features: dict | None = None):
        self.acc = [a.split("|")[1] if a.count("|") >= 2 else a for a in accessions]
        self.raw_acc = list(accessions)
        self.seq = list(sequences)
        self.il = [s.replace("I", "L") for s in self.seq]
        self.idx = {a: i for i, a in enumerate(self.acc)}
        self.feat = features or {}

    def find(self, pep_il: str, prot_idx):
        out = []
        for i in prot_idx:
            s = self.il[i]
            st = s.find(pep_il)
            while st >= 0:
                out.append((i, st, st + len(pep_il)))
                st = s.find(pep_il, st + 1)
        return out


def _end_status(P, x, y, fe):
    n = len(P)
    nt = x == 0 or (P[x - 1] in "KR" and P[x] != "P")
    ct = y == n or (P[y - 1] in "KR" and P[y] != "P")
    nbio = (not nt) and ((x == 1 and P[0] == "M") or x in fe.get("n", ()))
    cbio = (not ct) and (y in fe.get("c", ()))
    return nt, ct, nbio, cbio


def call(pep_il: str, prot_idx, db: ProteinDB) -> dict:
    occ = db.find(pep_il, prot_idx)
    if not occ:
        return {"call": "unmapped", "protein": "", "start": -1, "end": -1, "cut_site": -1, "cut_side": "",
                "orig_seq": pep_il, "n_occ": 0, "bio_reason": ""}
    best, rank = None, 9
    for i, x, y in occ:
        P = db.il[i]
        fe = db.feat.get(db.acc[i], {})
        nt, ct, nbio, cbio = _end_status(P, x, y, fe)
        n_ok = nt or nbio
        c_ok = ct or cbio
        if nt and ct:
            r, c = 0, "tryptic"
        elif n_ok and c_ok:
            r, c = 1, "bio_terminal"
        elif n_ok or c_ok:
            r, c = 2, "half_tryptic"
        else:
            r, c = 3, "nonspecific"
        if r < rank:
            rank, best = r, (i, x, y, c, nt, ct, nbio, cbio)
    i, x, y, c, nt, ct, nbio, cbio = best
    side, site = "", -1
    if c == "half_tryptic":
        if not (nt or nbio):
            side, site = "N", x
        else:
            side, site = "C", y
    reason = ";".join(k for k, v in (("Nbio", nbio), ("Cbio", cbio)) if v)
    if nbio and x == 1 and db.il[i][0] == "M":
        reason = reason.replace("Nbio", "Nbio_MetRemoval")
    return {"call": c, "protein": db.acc[i], "start": x + 1, "end": y, "cut_site": site, "cut_side": side,
            "orig_seq": db.seq[i][x:y], "n_occ": len(occ), "bio_reason": reason}


def annotate(peptides: pd.DataFrame, db: ProteinDB) -> pd.DataFrame:
    """peptides: columns target_seq (I/L collapsed) and prot_idx (list of protein indices)."""
    rows = [call(s, p, db) for s, p in zip(peptides["target_seq"], peptides["prot_idx"])]
    return pd.concat([peptides.reset_index(drop=True), pd.DataFrame(rows)], axis=1)


def insource_flags(prec: pd.DataFrame, tol_min: float = 0.1) -> pd.Series:
    """prec: one row per (run, precursor) with columns run, target_seq, call, rt_apex (minutes).  Returns a boolean
    Series (index aligned with prec): half-tryptic rows co-eluting with a prefix/suffix parent of the same run."""
    flag = np.zeros(len(prec), bool)
    for run, g in prec.groupby("run"):
        by_seq = defaultdict(list)
        for s, rt in zip(g["target_seq"], g["rt_apex"]):
            if np.isfinite(rt):
                by_seq[s].append(rt)
        seqs = list(by_seq)
        pre, suf = defaultdict(list), defaultdict(list)
        for s in seqs:
            for k in range(7, len(s)):
                pre[s[:k]].append(s)
                suf[s[-k:]].append(s)
        for j, (i, s, c, rt) in enumerate(zip(g.index, g["target_seq"], g["call"], g["rt_apex"])):
            if c != "half_tryptic" or not np.isfinite(rt):
                continue
            parents = pre.get(s, []) + suf.get(s, [])
            if any(abs(prt - rt) <= tol_min for p in parents for prt in by_seq[p]):
                flag[prec.index.get_loc(i)] = True
    return pd.Series(flag, index=prec.index)
