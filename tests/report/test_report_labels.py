"""A5.11.7 wording lint (bench/t1/wording_lint.py) over claims.md fixtures, REPORT.md and release/*.

REPORT.md and release/ files are linted when present.  The banned-phrase and "validated" checks always apply; the
structure checks (L-S in abstract and headline, table captions, "at 1% FDR" tags, headline not a number) apply to
REPORT.md once it has a T1 section (a heading containing "T1") and to release files that mention T1.  The rendered
labels that exempt quoted text are taken from the committed claims.json of the T1 look (results/*_T1_sealed/claims/
claims.json), and from config/a5.yaml's fixed label texts before that.  A failing lint blocks any release (A5.11.7)."""
import glob
import json
import os

import pytest
import yaml

from bench.t1 import wording_lint as WL

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def _labels():
    cfg = yaml.safe_load(open(os.path.join(ROOT, "config", "a5.yaml")))
    rendered = [v["text"] for v in cfg["labels"].values()]
    ls = [cfg["labels"]["L-S"]["text"]]
    for p in sorted(glob.glob(os.path.join(ROOT, "results", "*_T1_sealed", "**", "claims.json"), recursive=True)):
        d = json.load(open(p))
        rendered += list((d.get("labels_rendered") or {}).values())
        if (d.get("labels_rendered") or {}).get("L-S"):
            ls.append(d["labels_rendered"]["L-S"])
    return {"rendered": rendered, "L-S": ls}


@pytest.mark.parametrize("path", [p for p in [os.path.join(ROOT, "REPORT.md")] if os.path.exists(p)] +
                         sorted(glob.glob(os.path.join(ROOT, "release", "**", "*.*"), recursive=True)))
def test_report_and_release(path):
    kind = "report" if path.endswith("REPORT.md") else "release"
    v = WL.lint(open(path, errors="replace").read(), kind, _labels())
    assert not v, "\n".join(v[:20])


@pytest.mark.parametrize("phrase", ["entrapment-validated", "entrapment-verified", "empirically verified",
                                    "validated FDR", "FDR control", "controls FDR", "controlled at 1%",
                                    "FDR-controlled", "gate-passed", "contaminated nulls", "gate-validated",
                                    "validated", "v1 is calibrated at 1% FDR",
                                    "no FDR-control claim is made"])
def test_banned(phrase):
    assert WL.banned(f"The list of v1 is {phrase} here.")


@pytest.mark.parametrize("ok", ["PSM counts are labelled not FDR-controlled.", "not validated by null tests",
                                "Secondary arm; not gate-validated", "the comparator's rank ≥ 2 IDs were "
                                "miscalibrated on dev",
                                "validated by A6 (test 2), after A1 FAILED", "mass error was calibrated per run"])
def test_allowed(ok):
    assert not WL.banned(ok), WL.banned(ok)


def test_quoted_label_exemption_only_exact():
    lab = "No OT-HCD scorer has gate-validated FDR claims."
    assert not WL.banned(f"[L-D] “{lab}”", [lab])
    assert WL.banned(f"[L-D] “v1 has gate-validated FDR claims.”", [lab])


def _claims(head="[L-S] “S” P1 is met.", caption="**Table P1** [L-G v1] “G” [L-E] “E” "
            "[L-M] “M”", extra=""):
    return (f"# x\n\n## Headline\n\n{head}\n\n## Abstract\n\n[L-S] “S”\n\n## T\n\n{caption}\n| a |\n|---|\n"
            f"| 1 |\n\n{extra}\n")


LAB = {"rendered": ["S", "G", "E", "M", "Q"], "L-S": ["S"]}


def test_structure_ok():
    assert WL.lint(_claims(), "claims", LAB) == []


@pytest.mark.parametrize("kw,needle", [
    ({"head": "1,234 more IDs. [L-S] “S”"}, "number"),
    ({"head": "P1 is met."}, "L-S"),
    ({"caption": "**Table P1** [L-E] “E” [L-M] “M”"}, "[L-G"),
    ({"caption": "**Table P1** [L-G v1] “G” [L-E] “E”"}, "[L-M] or [L-Q]"),
    ({"caption": "Some text"}, "no caption"),
    ({"caption": "**Table P1** [L-G v1] “G” [L-E] “E” [L-M] “made up”"}, "not a rendered"),
    ({"extra": "v1 finds 10% more at 1% FDR than Sage."}, "[L-M]/[L-Q]"),
])
def test_structure_violations(kw, needle):
    v = WL.lint(_claims(**kw), "claims", LAB)
    assert any(needle in x for x in v), v


def test_fdr_phrase_with_tag_ok():
    assert WL.lint(_claims(extra="v1 finds more at 1% FDR [L-Q] “Q”."), "claims", LAB) == []


def test_report_without_t1_section_only_banned():
    assert WL.lint("# Report\n\n## Dev\n\nNo sealed numbers yet.\n", "report", LAB) == []
    assert WL.lint("# Report\n\n## Dev\n\nv1 controls FDR.\n", "report", LAB)
    assert WL.lint("# Report\n\n## T1 results\n\n12 more.\n", "report", LAB)
