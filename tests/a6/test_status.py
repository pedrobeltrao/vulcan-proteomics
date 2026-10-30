"""A6.2 status file: closed enum; absent / unverifiable -> 'stopped: node time limit'."""
import json
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts", "a6"))
import status as S  # noqa: E402


def test_enum_is_the_closed_a6_set():
    assert S.ENUM == ("PASS", "FAILED", "not evaluable", "invalid pool", "not completed", "label-free stop",
                      "not tagged", "stopped: node time limit")
    assert S.STATUS_PATH.endswith("results/2026-10-28_A6_hn1/a6_status.json")


@pytest.mark.parametrize("bad", ["pass", "Passed", "OK", "", "not evaluable (invalid pool)", None])
def test_writer_refuses_values_outside_the_enum(tmp_path, bad):
    with pytest.raises(ValueError):
        S.write_status(bad, path=str(tmp_path / "s.json"))
    assert not (tmp_path / "s.json").exists()


def test_pass_needs_commit_and_output(tmp_path):
    with pytest.raises(ValueError):
        S.write_status("PASS", path=str(tmp_path / "s.json"))
    with pytest.raises(ValueError):
        S.write_status("FAILED", step5_commit="abc", path=str(tmp_path / "s.json"))


def test_absent_file_is_node_time_limit(tmp_path):
    r = S.read_status(str(tmp_path / "missing.json"), require_committed=False)
    assert r["status"] == "stopped: node time limit" and not r["verified"]


def test_unparseable_or_outside_enum_is_node_time_limit(tmp_path):
    p = tmp_path / "s.json"
    p.write_text("{not json")
    assert S.read_status(str(p), require_committed=False)["status"] == "stopped: node time limit"
    p.write_text(json.dumps({"status": "pass"}))
    assert S.read_status(str(p), require_committed=False)["status"] == "stopped: node time limit"


def test_roundtrip_label_free_stop(tmp_path):
    p = str(tmp_path / "s.json")
    S.write_status("not tagged", note="infeasible", path=p)
    r = S.read_status(p, require_committed=False)
    assert r["status"] == "not tagged" and r["verified"]


def test_uncommitted_file_is_node_time_limit(tmp_path):
    p = str(tmp_path / "s.json")
    S.write_status("label-free stop", path=p)
    assert S.read_status(p, require_committed=True)["status"] == "stopped: node time limit"


def test_pass_with_tampered_gate_output_is_rejected(tmp_path, monkeypatch):
    out = os.path.join(S.ROOT, "results", "2026-10-28_A6_hn1", "_test_gate_output.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    try:
        open(out, "w").write('{"x": 1}')
        p = str(tmp_path / "s.json")
        S.write_status("PASS", step5_commit="a" * 40, gate_output=out, path=p)
        assert S.read_status(p, require_committed=False)["status"] == "PASS"
        open(out, "w").write('{"x": 2}')
        assert S.read_status(p, require_committed=False)["status"] == "stopped: node time limit"
    finally:
        os.remove(out)
