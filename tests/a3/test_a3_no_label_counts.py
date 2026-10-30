"""A2.9.4 extended to scripts/a3/ (A3.2.2 "Commit"): no code that runs before step 5 (step 4b and the A-syn library)
may read gate-DB labels or print/store per-label counts; every scripts/a3 file is classified; dev-only scripts never
touch DN2 paths; the step-5 driver checks its committed records before any statistic."""
import os
import re

import pytest

from test_no_label_counts import COUNT_TOKENS, GATE_LABEL_TOKENS, OWN_LABEL_TOKENS, _code_lines, _enclosing, _funcs  # noqa: F401,E501

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PRE_STEP5 = ["scripts/a3/step4b_chimeric.py", "vulcan/gates/a3_asyn.py", "scripts/a3/build_offsets.py"]
STEP5 = ["scripts/a3/gates_a3.py"]
DEV = ["scripts/a3/dev_g4.py", "scripts/a3/dev_armM.py", "scripts/a3/dev_mc_v2sep.py", "scripts/a3/dev_truncation.py",
       "scripts/a3/g0_freeze.py", "scripts/a3/dryrun_compare.py"]
GUARDS = ["scripts/a3/v1_closure_guard.py"]
# own-species certification (A2.2a) is recomputed with the UNCHANGED scripts/a2/own_search functions (which hold the
# only own-label access, guarded by a2_admission.own_labels); step 4b calls them by name only
OWN_CALLS = {("scripts/a3/step4b_chimeric.py", "certified_A")}


@pytest.mark.parametrize("path", PRE_STEP5)
def test_pre_step5_never_touches_gate_labels(path):
    tree, lines = _code_lines(path)
    funcs = _funcs(tree)
    bad = [(i, c.strip()) for i, c in lines if GATE_LABEL_TOKENS.search(c) or OWN_LABEL_TOKENS.search(c)]
    assert not bad, bad


@pytest.mark.parametrize("path", PRE_STEP5)
def test_pre_step5_has_no_per_label_count_names(path):
    _, lines = _code_lines(path)
    bad = [(i, c.strip()) for i, c in lines if COUNT_TOKENS.search(c)]
    assert not bad, bad


@pytest.mark.parametrize("path", PRE_STEP5)
def test_pre_step5_does_not_import_gate_statistics(path):
    import ast
    tree = ast.parse(open(os.path.join(ROOT, path)).read())
    mods = []
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            mods += [a.name for a in n.names]
        elif isinstance(n, ast.ImportFrom):
            mods += [n.module or ""] + [f"{n.module}.{a.name}" for a in n.names]
    bad = [m for m in mods if re.search(r"a2_stats|a3_stats|gates\.stats|gates_a1|gates_a2|gates_a3|seqclf|m3_data", m)]
    assert not bad, (path, bad)


def test_step4b_own_search_only_via_certification():
    src = open(os.path.join(ROOT, "scripts/a3/step4b_chimeric.py")).read()
    tree, lines = _code_lines("scripts/a3/step4b_chimeric.py")
    funcs = _funcs(tree)
    uses = [(i, _enclosing(funcs, i)) for i, c in lines if "OS." in c]
    assert uses and all(("scripts/a3/step4b_chimeric.py", f) in OWN_CALLS for _, f in uses), uses
    assert "certify_and_markers" in src


def test_every_scripts_a3_file_is_classified():
    files = sorted(f"scripts/a3/{f}" for f in os.listdir(os.path.join(ROOT, "scripts/a3")) if f.endswith(".py"))
    unclassified = [f for f in files if f not in PRE_STEP5 + STEP5 + DEV + GUARDS]
    assert not unclassified, "new scripts/a3 file must be classified for the A2.9.4 checks: %s" % unclassified


@pytest.mark.parametrize("path", DEV + ["reviews/checks/a3_power.py"])
def test_dev_and_power_scripts_never_touch_dn2(path):
    src = open(os.path.join(ROOT, path)).read()
    assert not re.search(r"vulcan_a2|VULCAN_A2_WORK|A2mix|A2_dn2|dn2\.tsv|data/sealed|sealed_run", src), path


def test_power_reads_no_data():
    src = open(os.path.join(ROOT, "reviews/checks/a3_power.py")).read()
    assert not re.search(r"labels\.arrow|load_store|parquet|data/raw|/tmp/vulcan", src)


def test_gates_a3_checks_records_before_statistics():
    src = open(os.path.join(ROOT, "scripts/a3/gates_a3.py")).read()
    main = src[src.index("def main():"):]
    i_freeze = main.index("require_freeze_commit(")
    i_rec = main.index('require_committed(a.record')
    i_v1 = main.index('require_committed(a.v1_results')
    for call in ("collect(cfg)", "g3_chimeric(cfg)", "collect_ce_s(cfg)"):
        assert i_freeze < main.index(call) and i_rec < main.index(call) and i_v1 < main.index(call)
    assert main.index("require_committed(a.step4b_record") < main.index("finalize(cfg")


def test_step4b_refuses_outside_freeze():
    src = open(os.path.join(ROOT, "scripts/a3/step4b_chimeric.py")).read()
    main = src[src.index("def main():"):]
    assert main.index("require_freeze(cfg)") < main.index("phase_ms1()")
