"""A3.1.3 freeze semantics (vulcan.gates.a3_config.freeze_status) and the A3.3.4 LOO backstop (bench/t1/a3.py)."""
import copy

import numpy as np

from vulcan.gates import a3_config

CFG = a3_config.load()


def test_freeze_at_head_trees_accepted_and_changed_path_rejected():
    c = copy.deepcopy(CFG)
    head = a3_config._git("rev-parse", "HEAD")
    paths = list(a3_config.FROZEN_TREES)          # A3.1.3 + A5.11.6
    c["freeze"]["commit"] = head
    c["freeze"]["tree_hashes"] = a3_config.tree_hashes(head, paths)
    st = a3_config.freeze_status(c)
    assert st["freeze_is_ancestor"] and st["tree_hashes_match"] and st["head_tree_hashes_match"] and st["code_at_freeze"]
    # an older commit whose frozen code differs from HEAD's: recorded hashes are that commit's, HEAD's differ
    old = "prereg-v0.4.2"
    c["freeze"]["commit"] = old
    c["freeze"]["tree_hashes"] = a3_config.tree_hashes(a3_config._git("rev-parse", old), paths)
    st = a3_config.freeze_status(c)
    assert st["freeze_is_ancestor"] and st["tree_hashes_match"] and not st["head_tree_hashes_match"]
    assert st["code_at_freeze"] is False and st["ok"] is False


def test_freeze_requires_a5_extended_tree_list():
    c = copy.deepcopy(CFG)
    head = a3_config._git("rev-parse", "HEAD")
    c["freeze"]["commit"] = head
    c["freeze"]["tree_hashes"] = a3_config.tree_hashes(head, ["vulcan", "bench", "scripts/a3", "scripts/e7"])
    st = a3_config.freeze_status(c)
    assert st["ok"] is False and "scripts/m3e" in st["reason"] and "scripts/baselines" in st["reason"]


def _curves(rng, R=6, J=5, S=3):
    from bench.t1 import evaluate as EV
    G = len(EV.GRID)
    out = {}
    for m, scale in (("num", 1.15), ("cmp", 1.0)):
        base = np.cumsum(rng.uniform(0, 50, size=G)) * scale
        ids = {lv: base[None, None, :] * rng.uniform(0.9, 1.1, size=(R, J, 1)) for lv in EV.LEVELS}
        cnt = {}
        for lv in EV.LEVELS:
            n_t = ids[lv][:, :, None, :].repeat(S, 2)
            n_e = n_t * np.linspace(0.0, 0.08, G)[None, None, None, :] * rng.uniform(0.8, 1.2, size=(R, J, S, 1))
            cnt[lv] = np.stack([n_t, n_e, 0.1 * n_e, 0.05 * n_e], axis=-1)
        out[m] = {"ids": ids, "cnt": cnt}
    return out


def test_loo_full_ratio_equals_noninferiority_and_loo_is_min():
    from bench.t1 import a3 as A3
    from bench.t1 import evaluate as EV
    curves = _curves(np.random.default_rng(5))
    methods = ["num", "cmp"]
    out = {"methods": {m: {"verdicts": {lv: {} for lv in EV.LEVELS}} for m in methods}}
    ni = EV.noninferiority(out, curves, methods, "num")
    loo = A3.loo_ratios(curves, "num", ["cmp"])
    assert A3.check_loo_full_equals_noninferiority(loo, ni) > 0
    d = loo["cmp"]["peptidoform"]["matched_fdp"]["0.01"]
    assert len(d["loo"]) == 6 and d["loo_min"] == min(d["loo"])
    # LOO of run 0 by hand
    keep = np.arange(1, 6)
    k = int(np.flatnonzero(np.isclose(EV.GRID, 0.01))[0])
    own = lambda m: curves[m]["ids"]["peptidoform"][keep][:, :, k].sum(0).mean()   # noqa: E731
    assert np.isclose(loo["cmp"]["peptidoform"]["own_q"]["0.01"]["loo"][0], own("num") / own("cmp"))


def test_loo_undefined_with_one_run():
    from bench.t1 import a3 as A3
    curves = _curves(np.random.default_rng(7), R=1)
    d = A3.loo_ratios(curves, "num", ["cmp"])["cmp"]["peptidoform"]["own_q"]["0.01"]
    assert d["loo"] == [] and d["loo_min"] is None


def test_a3_method_set_and_r_rule():
    from bench.t1 import a3 as A3
    assert A3.A1_METHODS == ["vulcan", "sage", "sage+mokapot", "sage+pred+mokapot", "comet+mokapot"]
    assert set(A3.A3_COMPARATORS) == {"sage", "sage+mokapot", "sage+pred+mokapot", "comet+mokapot", "sage_chim2",
                                      "sage_chim2+pred+mokapot"}
    assert A3.SAGE_CHIM_R == (2,) and A3.A3_NUMERATORS.split(",")[0] == "vulcan"


def test_closure_guard_follows_python_dash_m_subprocesses():
    import v1_closure_guard as V
    src = 'import subprocess, sys\nsubprocess.run([sys.executable, "-m", "vulcan.cohort", "--x"])\n'
    names = V.imports_of(src, "bench/t1/ours.py")
    assert "vulcan.cohort" in names and "vulcan.cohort.__main__" in names
    cl = V.closure("HEAD", ["bench/t1_driver.py"])
    assert {"vulcan/cohort.py", "vulcan/ms1/__main__.py", "scripts/e7/strata.py"} <= cl
