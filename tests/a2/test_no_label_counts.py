"""A2.9.4: no code under scripts/a2/ may print or store per-label counts before step 5, and admission must not read
gate-DB labels (P1/P2).  Static checks on the admission-phase sources + runtime checks on their outputs.

Admission-phase files (steps 2-4): select_dn2, own_search, admit, build_mixtures, resolve_panels and the library
modules they use.  Step-5 / whitelisted files: gates_a2.py (the gate run itself) and host_threshold.py (reads the
labels of the two human HOST runs only, stores two scalar thresholds)."""
import ast
import os
import re

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
ADMISSION = ["scripts/a2/select_dn2.py", "scripts/a2/own_search.py", "scripts/a2/admit.py",
             "scripts/a2/build_mixtures.py", "scripts/a2/resolve_panels.py", "vulcan/gates/a2_admission.py",
             "vulcan/gates/a2_identity.py", "vulcan/gates/a2_fs.py", "vulcan/gates/a2_panels.py",
             "vulcan/gates/a2_config.py"]
STEP5 = ["scripts/a2/gates_a2.py", "scripts/a2/host_threshold.py"]
# own-species index access (A2.2a / A2.3 step 3): allowed only in these functions
OWN_LABEL_FUNCS = {("vulcan/gates/a2_admission.py", "own_labels"), ("vulcan/gates/a2_admission.py", "certify_run"),
                   ("vulcan/gates/a2_admission.py", "qvalues_psm"), ("scripts/a2/own_search.py", "certify_and_markers")}
# hashing the frozen classifier pickle (A2.6 item 4) is not a label read
HASH_ONLY = {("vulcan/gates/a2_config.py", "check_seqclf"), ("scripts/a2/select_dn2.py", "main")}
GATE_LABEL_TOKENS = re.compile(r"labels\.arrow|proteins\.arrow|read_labels|m3_data|\bMD\b|seqclf|logit")
OWN_LABEL_TOKENS = re.compile(r"\bis_decoy\b|\btdc_pair\b|labels\.arrow|proteins\.arrow")
COUNT_TOKENS = re.compile(r"n_target|n_decoy|frac_target|target_win|\bwin_t\b|\bn_tail\b|tail_target|is_target|"
                          r"n_T\b|n_D\b|target_fraction|decoy_fraction")


def _funcs(tree):
    out = []
    for n in ast.walk(tree):
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)):
            out.append((n.lineno, n.end_lineno, n.name))
    return out


def _enclosing(funcs, line):
    best = None
    for a, b, name in funcs:
        if a <= line <= b and (best is None or a >= best[0]):
            best = (a, name)
    return best[1] if best else None


def _code_lines(path):
    """(line number, text) of code with docstrings/comments removed."""
    src = open(os.path.join(ROOT, path)).read()
    tree = ast.parse(src)
    doc = set()
    for n in ast.walk(tree):
        if isinstance(n, (ast.Module, ast.FunctionDef, ast.ClassDef, ast.AsyncFunctionDef)) and n.body \
                and isinstance(n.body[0], ast.Expr) and isinstance(getattr(n.body[0], "value", None), ast.Constant) \
                and isinstance(n.body[0].value.value, str):
            doc.update(range(n.body[0].lineno, n.body[0].end_lineno + 1))
    lines = []
    for i, line in enumerate(src.splitlines(), 1):
        if i in doc:
            continue
        code = line.split("#", 1)[0]
        lines.append((i, code))
    return tree, lines


@pytest.mark.parametrize("path", ADMISSION)
def test_admission_never_touches_gate_labels(path):
    tree, lines = _code_lines(path)
    funcs = _funcs(tree)
    bad = []
    for i, code in lines:
        m = GATE_LABEL_TOKENS.search(code)
        f = (path, _enclosing(funcs, i))
        if m and f not in OWN_LABEL_FUNCS and not (f in HASH_ONLY and m.group(0) == "seqclf"):
            bad.append((i, code.strip()))
        if OWN_LABEL_TOKENS.search(code) and (path, _enclosing(funcs, i)) not in OWN_LABEL_FUNCS:
            bad.append((i, code.strip()))
    assert not bad, bad


@pytest.mark.parametrize("path", ADMISSION)
def test_admission_has_no_per_label_count_names(path):
    _, lines = _code_lines(path)
    bad = [(i, c.strip()) for i, c in lines if COUNT_TOKENS.search(c)]
    assert not bad, bad


@pytest.mark.parametrize("path", ADMISSION)
def test_admission_does_not_import_gate_statistics(path):
    src = open(os.path.join(ROOT, path)).read()
    assert "a2_stats" not in src and "gates import stats" not in src and "gates_a1" not in src


def test_own_label_access_only_through_guard():
    # own_search reads own labels only via a2_admission.own_labels (which refuses non-own indexes)
    src = open(os.path.join(ROOT, "scripts/a2/own_search.py")).read()
    assert "AD.own_labels(idx)" in src and "labels.arrow" not in src


def test_every_scripts_a2_file_is_classified():
    files = sorted(f"scripts/a2/{f}" for f in os.listdir(os.path.join(ROOT, "scripts/a2")) if f.endswith(".py"))
    unclassified = [f for f in files if f not in ADMISSION + STEP5 + ["scripts/a2/power_checks.py"]]
    assert not unclassified, "new scripts/a2 file must be added to the A2.9.4 label checks: %s" % unclassified


def test_gate_run_checks_committed_membership_first():
    src = open(os.path.join(ROOT, "scripts/a2/gates_a2.py")).read()
    main = src[src.index("def main():"):]
    assert main.index("require_final_membership(") < main.index("collect(cfg)")
    assert main.index("require_final_membership(") < main.index("g3(cfg)")


def test_power_checks_reads_no_data():
    src = open(os.path.join(ROOT, "scripts/a2/power_checks.py")).read()
    assert not re.search(r"labels\.arrow|load_store|parquet|data/raw|/tmp/vulcan", src)
