"""A6.4.6.3 label guard (static part) and 'no per-label output from steps 2-4' (A6.8.4), as
tests/a2/test_no_label_counts.py.

Classification of scripts/a6/*.py:
  LABEL_FREE   never reference labels.arrow / proteins.arrow / is_decoy / tdc_pair / m3_data / seqclf, never compute a
               TDC threshold or competition q-value (double_competition / tail_threshold), never import gate statistics;
  CERT         the certification module (admit_hn1) and its dev driver (dry_run): may read full-index labels, only
               through admit_hn1.LabelView at host full-index rows (runtime guard tested in tests/a6/test_admit.py);
  STEP5        gates_hn1: the only place labels of null rows are read.
"""
import ast
import os
import re

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
A6 = os.path.join(ROOT, "scripts", "a6")
LABEL_FREE = ["select_hn1.py", "mixtures_hn1.py", "reduced_index.py", "status.py", "a6_stats.py",
              "a6_power_dryrun.py", "throughput.py"]
CERT = ["admit_hn1.py", "dry_run.py"]
STEP5 = ["gates_hn1.py"]
LABEL_TOKENS = re.compile(r"labels\.arrow|proteins\.arrow|\bis_decoy\b|\btdc_pair\b|m3_data|\bMD\b|seqclf|"
                          r"read_labels|full_labels|\.labels\b")
TDC_TOKENS = re.compile(r"double_competition|tail_threshold|qvalues_psm|certify_run")
COUNT_TOKENS = re.compile(r"n_target|n_decoy|frac_target|target_win|\bwin_t\b|\bn_tail\b|tail_target|is_target|"
                          r"\bn_T\b|\bn_D\b|target_fraction|decoy_fraction")


def _code_lines(path):
    src = open(path).read()
    tree = ast.parse(src)
    doc = set()
    for n in ast.walk(tree):
        if isinstance(n, (ast.Module, ast.FunctionDef, ast.ClassDef, ast.AsyncFunctionDef)) and n.body \
                and isinstance(n.body[0], ast.Expr) and isinstance(getattr(n.body[0], "value", None), ast.Constant) \
                and isinstance(n.body[0].value.value, str):
            doc.update(range(n.body[0].lineno, n.body[0].end_lineno + 1))
    out = []
    for i, line in enumerate(src.splitlines(), 1):
        if i in doc:
            continue
        out.append((i, line.split("#", 1)[0]))
    return out


def _present(names):
    return [f for f in names if os.path.exists(os.path.join(A6, f))]


@pytest.mark.parametrize("f", _present(LABEL_FREE))
def test_label_free_modules_never_touch_labels(f):
    bad = [(i, c.strip()) for i, c in _code_lines(os.path.join(A6, f)) if LABEL_TOKENS.search(c)]
    assert not bad, bad


@pytest.mark.parametrize("f", [x for x in _present(LABEL_FREE) if x != "a6_stats.py"])
def test_label_free_modules_compute_no_tdc_threshold(f):
    bad = [(i, c.strip()) for i, c in _code_lines(os.path.join(A6, f)) if TDC_TOKENS.search(c)]
    assert not bad, bad


@pytest.mark.parametrize("f", _present(["select_hn1.py", "mixtures_hn1.py", "reduced_index.py", "status.py"]
                                       + CERT))
def test_steps_2_to_4_have_no_per_label_count_names(f):
    bad = [(i, c.strip()) for i, c in _code_lines(os.path.join(A6, f)) if COUNT_TOKENS.search(c)]
    assert not bad, bad


@pytest.mark.parametrize("f", _present(["select_hn1.py", "mixtures_hn1.py", "reduced_index.py"]))
def test_steps_2_to_4_do_not_import_gate_statistics(f):
    src = "\n".join(c for _, c in _code_lines(os.path.join(A6, f)))
    assert not re.search(r"import a6_stats|from a6_stats import (?!cap_factors)|a2_stats|gates import stats as|"
                         r"gates_hn1", src), f


@pytest.mark.parametrize("f", _present(CERT))
def test_certification_module_does_not_import_gate_statistics(f):
    # vulcan.gates.stats.tail_threshold is allowed here: the run's own host 1% threshold (full-index search)
    src = "\n".join(c for _, c in _code_lines(os.path.join(A6, f)))
    assert not re.search(r"import a6_stats|from a6_stats import|a2_stats|gates_hn1|mixtures_hn1", src), f


def test_mixtures_use_label_free_competition_only():
    src = "\n".join(c for _, c in _code_lines(os.path.join(A6, "mixtures_hn1.py")))
    assert "psm_competition" in src and "double_competition" not in src


def test_only_gates_hn1_reads_null_row_labels():
    src = open(os.path.join(A6, "gates_hn1.py")).read()
    # the label read is inside collect(), which sets the completed-test flag before reading
    body = src[src.index("def collect("):src.index("def evaluate(")]
    assert body.index("LabelRead.started = True") < body.index("MD.labels(")
    main = src[src.index("def main():"):]
    assert main.index("guards(") < main.index("collect(")


def test_every_scripts_a6_file_is_classified():
    files = sorted(f for f in os.listdir(A6) if f.endswith(".py"))
    unclassified = [f for f in files if f not in LABEL_FREE + CERT + STEP5]
    assert not unclassified, f"classify new scripts/a6 files for the A6.4.6.3 label guard: {unclassified}"
