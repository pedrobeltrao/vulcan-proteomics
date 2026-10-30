"""Wording lint for T1/T5 claim text (DESIGN §6.0d, Amendment A5.11.7).

    lint(text, kind, labels=None) -> list of violation strings (empty = pass)

kind:
  "claims"   claims.md written by bench/t1/claims_a5.py: every check, structure included.
  "report"   REPORT.md: banned phrases and "validated" everywhere; the structure checks (L-S in the abstract and the
             headline sentence, table captions, "at 1% FDR" tags, headline not starting with a number) apply once the
             file has a T1 section (a heading containing the token "T1").  Before the T1 look REPORT.md has no T1
             numbers, so there is nothing to label; the banned-phrase checks still bind (stricter than skipping).
  "release"  any file under release/: the same as "report", but the structure checks always apply when the text
             contains the token "T1" (release text has no sections to anchor on).

labels: {"rendered": [every label string the claim script rendered], "L-S": [the rendered L-S variant(s)]}.
  Banned phrases are exempt only inside a curly-quoted span “...” whose content is EXACTLY one rendered label
  (A5.11.7 "outside quoted label text").  Quoting arbitrary text does not exempt it.

Conventions of claims.md (and REPORT.md, which copies claims.md verbatim):
  * labels are written as  [L-X] “<label text>”
  * a table caption is the line directly above a markdown table and starts with "**Table"; a caption of a table with
    FDR-filtered numbers carries [L-G …] “…”, [L-E] “…” (or per-row L-E) and [L-M] or [L-Q] tags; a table without
    FDR-filtered numbers says "[no FDR-filtered numbers]" in its caption.
  * the headline is the first paragraph after a heading "## Headline" (claims.md) or the first paragraph of the T1
    section (REPORT.md); its first sentence starts with the L-S tag, never with a digit.
"""
from __future__ import annotations

import re

LQ, RQ = "“", "”"

# A5.11.7 banned phrases (case-insensitive), each with the exception the text allows
BANNED = [
    ("entrapment-validated", r"entrapment-validated", None),
    ("entrapment-verified", r"entrapment-verified", None),
    ("empirically verified", r"empirically[ -]verified", None),
    ("validated FDR", r"validated[ -]FDR", None),
    ("FDR control", r"FDR[ -]control(?!led)", None),       # "FDR-control claim" only inside quoted A5.9 labels
    ("controls FDR", r"controls[ -]FDR", None),
    ("controlled at", r"controlled at", None),
    ("FDR-controlled", r"FDR[ -]controlled", r"not FDR[ -]controlled"),
    ("gate-passed", r"gate[ -]passed", None),
    ("contaminated nulls", r"contaminated nulls?", None),
]
# "calibrated" is banned when applied to v1/v2_sep FDR: flagged in any sentence that also names FDR, q, a scorer or a
# T1/T5 list (stricter than deciding "applied to" by hand); "miscalibrated" (the E7b caveat) is not the word.
CALIBRATED = re.compile(r"(?<![a-z])calibrated", re.I)
CALIB_CONTEXT = re.compile(r"FDR|FDP|q-value|\bq\b|own q|v1|v2_sep|vulcan|static_gbm|scorer|T1|T5", re.I)
# "validated" only negated, or in the exact A5.9.3 pass phrases
VALIDATED = re.compile(r"[A-Za-z-]*validated", re.I)
VALIDATED_OK_BEFORE = re.compile(r"\bnot\s+$", re.I)
PASS_PHRASES = ["validated by A6 (test 2), after A1 FAILED"]
FDR_AT = re.compile(r"\bat (?:own q |matched (?:arm-P )?FDP̂? ?(?:= )?)?\d+(?:\.\d+)?\s?% ?(?:FDR|q|FDP)\b", re.I)


def _strip_quoted_labels(text: str, rendered) -> str:
    rendered = set(rendered or [])

    def rep(m):
        return " " * len(m.group(0)) if m.group(1) in rendered else m.group(0)
    return re.sub(LQ + "(.*?)" + RQ, rep, text, flags=re.S)


def _sentences(text):
    return re.split(r"(?<=[.;!?])\s+|\n\s*\n", text)


def banned(text: str, rendered=None) -> list:
    t = _strip_quoted_labels(text, rendered)
    out = []
    for name, pat, exc in BANNED:
        for m in re.finditer(pat, t, re.I):
            if exc:
                s = max(0, m.start() - 4)
                if re.match(exc, t[m.start():], re.I) or re.match(exc, t[s:], re.I):
                    continue
            out.append(f"banned phrase '{name}': ...{t[max(0, m.start() - 40):m.end() + 20]!r}...")
    for sent in _sentences(t):
        if CALIBRATED.search(sent) and CALIB_CONTEXT.search(sent):
            out.append(f"banned 'calibrated' applied to FDR/scorer: {sent.strip()[:120]!r}")
    for m in VALIDATED.finditer(t):
        w = m.group(0).lower()
        before = t[:m.start()]
        if w in ("validated", "gate-validated") and VALIDATED_OK_BEFORE.search(before):
            continue
        if any(t.startswith(p, m.start()) for p in PASS_PHRASES):
            continue
        out.append(f"'validated' outside the negated forms / A5.9.3 pass phrases: ...{t[max(0, m.start() - 40):m.end() + 20]!r}...")
    return out


def _section(text, heading_re):
    lines = text.split("\n")
    for i, l in enumerate(lines):
        if re.match(heading_re, l):
            lvl = len(l) - len(l.lstrip("#"))
            body = []
            for l2 in lines[i + 1:]:
                if l2.startswith("#") and len(l2) - len(l2.lstrip("#")) <= lvl:
                    break
                body.append(l2)
            return "\n".join(body)
    return None


def _first_paragraph(body):
    for p in re.split(r"\n\s*\n", body or ""):
        if p.strip():
            return p.strip()
    return ""


def _has_ls(text, ls_texts):
    return any((LQ + s + RQ) in text for s in (ls_texts or []))


def structure(text: str, kind: str, labels: dict) -> list:
    out = []
    rendered = set(labels.get("rendered", []))
    ls = labels.get("L-S", [])
    if kind == "claims":
        head = _first_paragraph(_section(text, r"^## Headline"))
        abstract = _section(text, r"^## Abstract")
    else:
        t1 = _section(text, r"^#{1,6} .*\bT1\b")
        head = _first_paragraph(t1)
        abstract = _section(text, r"^#{1,6} Abstract")
        if kind == "release" and t1 is None:
            head, abstract = text.strip(), text
    if not head:
        out.append("no headline paragraph")
    else:
        if re.match(r"[\s*_>\[(]*[0-9+\-−]", head) and not head.lstrip().startswith("[L-S]"):
            out.append("headline sentence starts with a number")
        if not head.lstrip().startswith("[L-S]") or not _has_ls(head, ls):
            out.append("headline does not start with the L-S label")
    if abstract is None or not _has_ls(abstract, ls):
        out.append("L-S label missing from the abstract")
    lines = text.split("\n")
    for i, l in enumerate(lines):
        if l.startswith("|") and (i == 0 or not lines[i - 1].startswith("|")):
            cap = lines[i - 1] if i > 0 else ""
            if not cap.startswith("**Table"):
                out.append(f"table at line {i + 1} has no caption")
                continue
            if "[no FDR-filtered numbers]" in cap:
                continue
            for tag in ("[L-G", "[L-E", ):
                if tag not in cap:
                    out.append(f"caption at line {i} lacks {tag}]")
            if "[L-M]" not in cap and "[L-Q]" not in cap:
                out.append(f"caption at line {i} lacks [L-M] or [L-Q]")
            for m in re.finditer(r"\[(L-[A-Z0-9]+)[^\]]*\]\s*" + LQ + "(.*?)" + RQ, cap):
                if m.group(2) not in rendered:
                    out.append(f"caption at line {i}: {m.group(1)} text is not a rendered label")
    for l in lines:
        if l.startswith("|") or l.startswith("**Table"):
            continue
        s = _strip_quoted_labels(l, rendered)
        for m in FDR_AT.finditer(s):
            if "[L-M]" not in l[m.end():] and "[L-Q]" not in l[m.end():]:
                out.append(f"'{m.group(0)}' not followed by an [L-M]/[L-Q] tag: {l.strip()[:120]!r}")
    return out


def lint(text: str, kind: str = "claims", labels: dict | None = None) -> list:
    labels = labels or {}
    out = banned(text, labels.get("rendered"))
    if kind == "claims":
        out += structure(text, kind, labels)
    elif kind == "report":
        if re.search(r"^#{1,6} .*\bT1\b", text, re.M):
            out += structure(text, kind, labels)
    elif kind == "release":
        if re.search(r"\bT1\b", text):
            out += structure(text, kind, labels)
    else:
        raise ValueError(kind)
    return out
