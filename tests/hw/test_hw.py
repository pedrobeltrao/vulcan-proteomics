import json
import os
import subprocess
import sys

from vulcan import hw

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def test_n_threads_affinity(monkeypatch):
    monkeypatch.delenv("VULCAN_THREADS", raising=False)
    assert hw.n_threads() == len(os.sched_getaffinity(0))


def test_n_threads_override(monkeypatch):
    monkeypatch.setenv("VULCAN_THREADS", "4")
    assert hw.n_threads() == 4
    monkeypatch.setenv("VULCAN_THREADS", "garbage")
    assert hw.n_threads() == len(os.sched_getaffinity(0))
    monkeypatch.setenv("VULCAN_THREADS", "0")
    assert hw.n_threads() >= 1


def test_pin_thread_pools(monkeypatch):
    for k in hw.THREAD_VARS:
        monkeypatch.setenv(k, "99")
    assert hw.pin_thread_pools(4) == 4
    for k in hw.THREAD_VARS:
        assert os.environ[k] == "4"
    assert hw.thread_env(2)["RAYON_NUM_THREADS"] == "2"


def test_fingerprint_keys():
    fp = hw.fingerprint()
    for k in ("hostname", "cpu_model", "affinity_count", "physical_cores_in_affinity", "avx512", "ram_gb", "gpus", "driver"):
        assert k in fp
    assert fp["affinity_count"] == len(os.sched_getaffinity(0))
    assert fp["physical_cores_in_affinity"] >= 1
    json.dumps(fp)


def test_hwinfo_cli(tmp_path):
    out = tmp_path / "hw.json"
    subprocess.check_call([sys.executable, os.path.join(ROOT, "bench", "hwinfo.py"), str(out)])
    assert json.load(open(out))["affinity_count"] >= 1
    s = subprocess.check_output([sys.executable, os.path.join(ROOT, "bench", "hwinfo.py")], text=True)
    assert json.loads(s)["hostname"]
