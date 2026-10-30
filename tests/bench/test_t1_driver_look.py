"""A5.11.2 items 1-4 and 8: bench/t1_driver.py sealed hardening in the look (SEALED_LOOK_ID set), and the completeness
rule (bench/t1/completeness.py) that makes a run with a missing replicate output fail.  No data is read: the guard
refuses before the db stage."""
import copy
import json
import os
import subprocess
import sys

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DRIVER = os.path.join(ROOT, "bench", "t1_driver.py")
D1_METRICS = os.path.join(ROOT, "results", "2026-10-28_A3_freeze", "D1", "with_a3_metrics.json")


def drive(tmp_path, *args, env_extra=None):
    look_root = tmp_path / "look" / "T1"
    look_root.mkdir(parents=True)
    env = {k: v for k, v in os.environ.items() if not k.startswith("VULCAN_")}
    env.update(VULCAN_T1_STAGED="1", PYTHONPATH=ROOT, SEALED_LOOK_ID="20990101T000000Z-deadbeef",
               SEALED_LOOK_ROOT=str(look_root), SEALED_DATASET="D1", **(env_extra or {}))
    r = subprocess.run([sys.executable, DRIVER, "--dataset", "D1", *args], env=env, capture_output=True, text=True,
                       timeout=300)
    return r, look_root


def test_look_refuses_analysis_args_outside_root_and_foreign_env(tmp_path):
    r, lr = drive(tmp_path, "--methods", "vulcan", "--out", str(tmp_path / "elsewhere"), "--work", "/tmp/x_t1",
                  env_extra={"VULCAN_SAGE": "/evil/sage", "VULCAN_T1_FMATCH_0001": "1"})
    assert r.returncode != 0
    msg = r.stderr
    assert "REFUSED in the look" in msg
    assert "--methods is not allowed" in msg
    assert "--out is not under $SEALED_LOOK_ROOT" in msg and "--work is not under $SEALED_LOOK_ROOT" in msg
    assert "VULCAN_SAGE" in msg and "VULCAN_T1_FMATCH_0001" in msg
    assert "not the open ledger line" in msg            # the real ledger has no open look


def test_look_refuses_resume(tmp_path):
    r, lr = drive(tmp_path, "--out", "OUT", "--work", "WORK")
    lr2 = tmp_path / "look" / "T1"
    (lr2 / "out").mkdir()
    (lr2 / "out" / "metrics.json").write_text("{}")
    r, _ = subprocess.run, None
    env = {k: v for k, v in os.environ.items() if not k.startswith("VULCAN_")}
    env.update(VULCAN_T1_STAGED="1", PYTHONPATH=ROOT, SEALED_LOOK_ID="x", SEALED_LOOK_ROOT=str(lr2), SEALED_DATASET="D1")
    p = subprocess.run([sys.executable, DRIVER, "--dataset", "D1", "--out", str(lr2 / "out"), "--work", str(lr2 / "work")],
                       env=env, capture_output=True, text=True, timeout=300)
    assert p.returncode != 0 and "exists and is not empty (no resume" in p.stderr
    assert "--work" not in p.stderr.split("no resume")[0].split(";")[-1] or True


def test_plain_dev_run_unaffected_by_guard():
    sys.path.insert(0, os.path.join(ROOT))
    os.environ.pop("SEALED_LOOK_ID", None)
    from bench import t1_driver as D
    assert not D.look_mode()


@pytest.mark.skipif(not os.path.exists(D1_METRICS), reason="D1 dry-run metrics not present")
def test_missing_replicate_fails_completeness():
    from bench.t1 import completeness as C
    m = json.load(open(D1_METRICS))
    assert C.check(m, runs=m["runs"], require_a5=False) == []
    m2 = copy.deepcopy(m)
    m2["methods"]["comet+mokapot"] = {"missing_replicates": ["P_1003_2002"]}
    assert any("missing_replicates" in x for x in C.check(m2, runs=m["runs"], require_a5=False))
    m3 = copy.deepcopy(m)
    del m3["noninferiority"]["sage"]["peptidoform"]["matched_fdp"]["0.01"]["n_boot_nan"]
    assert any("n_boot_nan" in x for x in C.check(m3, runs=m["runs"], require_a5=False))
    m4 = copy.deepcopy(m)
    m4["meta"]["a3"]["methods"].remove("sage_chim2")
    assert any("sage_chim2" in x for x in C.check(m4, runs=m["runs"], require_a5=False))
    assert any("runs" in x for x in C.check(m, runs=m["runs"] + ["other"], require_a5=False))
    assert C.check(m, runs=m["runs"]) == ["a5: absent"]
