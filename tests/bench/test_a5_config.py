"""config/a5.yaml pinned to the normative text (reviews/amendment_A5_draft.md rev 1, A5.11.5 "a test pins the label
texts to this document"): every label text and every A5.4 gate-table cell occurs verbatim in the draft, the A6 hook
(path, closed enum, absent rule) equals the A6 text (A6.2) and claims_a5.py, and the C9 record matches its file."""
import hashlib
import os
import subprocess

import yaml

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
CFG = yaml.safe_load(open(os.path.join(ROOT, "config", "a5.yaml")))
DRAFT = open(os.path.join(ROOT, "reviews", "amendment_A5_draft.md")).read()
A6 = open(os.path.join(ROOT, "reviews", "amendment_A6_draft.md")).read()

REQUIRED_LABELS = ["L-G_v1", "L-G_v2_sep", "C9_clause_nogo", "C9_clause_timebox", "C9_clause_a6", "L-E", "L-M", "L-Q",
                   "L-R", "L-S", "L-2", "L-D", "L-P", "L-S_last_clause_a6_verdict", "L-S_pass_last_clause",
                   "L-Q_pass_last_clause", "L-R_pass", "v2_sep_title", "v2_sep_title_a6_pass", "stratum1_title",
                   "void_title", "R4"]


def test_labels_verbatim():
    for k in REQUIRED_LABELS:
        assert k in CFG["labels"], k
    for k, v in CFG["labels"].items():
        assert v["text"] in DRAFT, k


def test_every_short_form_contains_failed():
    assert "FAILED" in CFG["labels"]["L-S"]["text"]
    assert "FAILED" in CFG["labels"]["L-S_pass_last_clause"]["text"]


def test_gate_tables_verbatim():
    for sc, n in (("v1", 12), ("v2_sep", 8)):
        rows = CFG["gates"][sc]["rows"]
        assert len(rows) == n
        for r in rows:
            assert f"| {r['gate']} | {r['data']} | {r['result']} | {r['source']} |" in DRAFT, r


def test_a6_hook_matches_a6_text_and_code():
    from bench.t1 import claims_a5 as C
    a6 = CFG["a6"]
    assert a6["status_path"] == C.A6_STATUS_PATH == "results/2026-10-28_A6_hn1/a6_status.json"
    assert f"`{a6['status_path']}`" in A6
    assert tuple(a6["enum"]) == C.A6_ENUM
    for e in a6["enum"]:
        assert f"`{e}`" in A6, e
    assert a6["absent_state"] == C.A6_ABSENT == "stopped: node time limit"


def test_c9_record():
    c9 = CFG["c9"]
    assert c9["verdict"] == "GO" and c9["commit"] == "3097ef5"
    blob = subprocess.run(["git", "-C", ROOT, "show", f"{c9['commit']}:{c9['verdict_file']}"], capture_output=True,
                          check=True).stdout
    assert hashlib.sha256(blob).hexdigest() == c9["verdict_file_sha256"]
    import json
    d = json.loads(blob)
    assert d["summary"]["f0.20/own/R3/kappa1.0/budget200"]["P_v1_minimum"] >= 0.8      # GO rule (A5.1.1)


def test_comparator_nonreproducibility_d33():
    v = CFG["comparator_nonreproducibility"]["values"]
    assert v["sage+pred+mokapot"] == {"arm_P_fdp_abs": 4.3e-4, "ratio_abs": 9.1e-4}
    assert v["sage+mokapot"] == {"ratio_abs": 2.3e-3}


def test_a3_1_8_record():
    r = CFG["a3_1_8_records"]
    assert hashlib.sha256(open(os.path.join(ROOT, r["path"]), "rb").read()).hexdigest() == r["sha256"]
