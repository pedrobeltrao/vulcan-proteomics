"""A5.11.1 item 10: bench/sealed_run.py one look in code (DESIGN §6.0d).  Every test runs in a temporary git repo via
--root, with dummy "sealed" files and dummy commands (no real data, no real ledger).  --dev is used so that the A3
freeze check (which needs the real vulcan/ tree) is skipped; every other check is the production code path."""
import datetime
import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import time

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SR = os.path.join(ROOT, "bench", "sealed_run.py")
sys.path.insert(0, os.path.join(ROOT, "bench"))
import sealed_run as SRM  # noqa: E402

DUMMY = r'''
import json, os, subprocess, sys, time
root, ctl = sys.argv[1], sys.argv[2]
c = json.load(open(ctl)) if os.path.exists(ctl) else {}
ds = os.environ["SEALED_DATASET"]
b = c.get(ds, {})
repo = os.getcwd()
rec = {"env": sorted(os.environ), "ledger": open(os.path.join(repo, "results/SEALED_LEDGER.md")).read(),
       "files_readable": all(os.access(p, os.R_OK) for p in os.environ["SEALED_FILES"].split()),
       "look_id": os.environ.get("SEALED_LOOK_ID"), "root": os.environ.get("SEALED_LOOK_ROOT")}
json.dump(rec, open(os.path.join(os.path.dirname(ctl), f"seen_{ds}.json"), "w"))
print("STAGE start", flush=True)
print("SECRET 12345 identifications", flush=True)
if b.get("sleep"):
    time.sleep(b["sleep"])
if b.get("leak"):
    open(os.path.join(repo, "results", "leaked.txt"), "w").write("x")
os.makedirs(os.path.join(root, "out"), exist_ok=True)
open(os.path.join(root, "out", "metrics.json"), "w").write("{}")
open(os.path.join(root, "out", "manifest.sha256"), "w").write("abc  metrics.json\n")
if b.get("crash"):
    raise ValueError("secret message 999")
print("STAGE done", flush=True)
sys.exit(b.get("exit", 0))
'''


def sh(*a, cwd=None, check=True):
    r = subprocess.run(list(a), cwd=cwd, capture_output=True, text=True)
    if check and r.returncode:
        raise RuntimeError(r.stderr)
    return r.stdout.strip()


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def write_cfg(repo, tmp, c9_commit, c9_sha, verdict="NO-GO", hashes=None, deviations=None):
    import yaml
    cfg = {
        "c9": {"verdict": verdict, "commit": c9_commit, "verdict_file": "results/c9/c9_power.json",
               "verdict_file_sha256": c9_sha},
        "look": {
            "commands": {ds: [sys.executable, "dummy/cmd.py", "{root}", str(tmp / "ctl" / "ctl.json")]
                         for ds in ("T5", "T1")},
            "purpose": {"T5": "arm M", "T1": "v1 primary + v2_sep secondary"},
            "env": {"PATH": os.environ["PATH"], "HOME": str(tmp)},
            "root": str(tmp / "lookroot"), "time_limit": "1h",
            "console_whitelist": [r"^STAGE \w+$"],
            "leak_paths": ["results"],
            "dev_files": {"T5": ["data/sealed/X5/a.raw"], "T1": ["data/sealed/X1/b.raw", "data/sealed/X1/c.raw"]},
        },
        "hashes": hashes if hashes is not None else {"dummy/cmd.py": sha(repo / "dummy" / "cmd.py")},
        "env_manifests": {},
        "deviations": deviations or [],
    }
    (repo / "config").mkdir(exist_ok=True)
    (repo / "config" / "a5.yaml").write_text(yaml.safe_dump(cfg))


def commit(repo, msg, *paths):
    sh("git", "add", *(paths or ["-A"]), cwd=repo)
    sh("git", "commit", "-q", "-m", msg, cwd=repo)


class RepoPath(type(__import__("pathlib").Path())):
    pass


@pytest.fixture
def repo(tmp_path):
    r = RepoPath(tmp_path / "repo")
    r.mkdir()
    sh("git", "init", "-q", cwd=r)
    sh("git", "config", "user.email", "t@t", cwd=r)
    sh("git", "config", "user.name", "t", cwd=r)
    (r / ".gitignore").write_text("data/sealed/\n")
    (r / "bench").mkdir()
    shutil.copy(os.path.join(ROOT, "bench", "env_manifest.py"), r / "bench" / "env_manifest.py")
    shutil.copy(SR, r / "bench" / "sealed_run.py")
    (r / "dummy").mkdir()
    (r / "dummy" / "cmd.py").write_text(DUMMY)
    (r / "results" / "c9").mkdir(parents=True)
    (r / "results" / "c9" / "c9_power.json").write_text(json.dumps({
        "summary": {SRM.C9_QUANTITY: {"P_v1_minimum": 0.5}},
        "walk_cells": [{"P_v1_minimum_by_run_budget": {"200": 0.5}}]}))
    sh_lines = []
    for rel in ("X5/a.raw", "X1/b.raw", "X1/c.raw"):
        p = r / "data" / "sealed" / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(rel.encode())
        sh_lines.append(f"{sha(p)}  sealed/{rel}")
    (r / "data" / "manifests").mkdir(parents=True)
    (r / "data" / "manifests" / "sealed.sha256").write_text("\n".join(sh_lines) + "\n")
    commit(r, "base")
    c9c = sh("git", "rev-parse", "HEAD", cwd=r)
    write_cfg(r, tmp_path, c9c, sha(r / "results" / "c9" / "c9_power.json"))
    commit(r, "config")
    for t in ("prereg-v0.4.3", "prereg-v0.4.4", "milestone-t1"):
        sh("git", "tag", t, cwd=r)
    (tmp_path / "ctl").mkdir()
    r.tmp = tmp_path
    r.c9c = c9c
    return r


def ctl(repo, **kw):
    (repo.tmp / "ctl" / "ctl.json").write_text(json.dumps(kw))


def run(repo, *args, env=None, timeout=120):
    e = dict(os.environ, **(env or {}))
    return subprocess.run([sys.executable, SR, "--root", str(repo), *args], capture_output=True, text=True, env=e,
                          timeout=timeout)


def look(repo, *extra, **kw):
    return run(repo, "--look", "a5", "--dev", *extra, **kw)


def ledger(repo):
    return (repo / "results" / "SEALED_LEDGER.md").read_text()


def test_plain_t1_and_t5_refused(repo):
    for ds in ("T1", "T5"):
        r = run(repo, "--dataset", ds, "--purpose", "x", "--", "true")
        assert r.returncode and "only by the one look" in r.stderr


def test_dev_refused_on_project_root():
    r = subprocess.run([sys.executable, SR, "--look", "a5", "--dev"], capture_output=True, text=True)
    assert r.returncode and "temporary clone" in r.stderr


def test_look_success_then_second_refused_even_after_rewrite(repo):
    ctl(repo)
    r = look(repo, env={"VULCAN_SAGE": "/evil/sage", "SECRET_VAR": "1"})
    assert r.returncode == 0, r.stderr + r.stdout
    # console whitelist: stage lines shown, data lines not
    assert "STAGE start" in r.stdout and "SECRET" not in r.stdout and "12345" not in r.stdout
    seen5 = json.load(open(repo.tmp / "ctl" / "seen_T5.json"))
    assert not os.path.exists(f"/tmp/vulcan_look_{seen5['look_id']}")          # scratch removed
    # env allow-list: caller variables dropped
    assert "VULCAN_SAGE" not in seen5["env"] and "SECRET_VAR" not in seen5["env"]
    assert {"SEALED_LOOK_ID", "SEALED_LOOK_ROOT", "SEALED_FILES"} <= set(seen5["env"])
    # both ledger lines written (one commit) before any file was unlocked
    lid = seen5["look_id"]
    assert f"| T5 |" in seen5["ledger"] and f"| T1 |" in seen5["ledger"] and seen5["ledger"].count(f"look={lid}") == 2
    assert seen5["files_readable"]
    opened = sh("git", "log", "--format=%s", cwd=repo).splitlines()
    assert any("opened (T5, T1)" in s for s in opened)
    led = ledger(repo)
    assert f"result: look={lid} exit=0" in led and "leak_test=pass" in led
    assert led.index("| T5 | result") < led.index("| T1 | result")
    # the files are relocked
    assert not os.access(repo / "data" / "sealed" / "X1" / "b.raw", os.R_OK)
    # second look refused (ledger-only commits are skipped: the refusal is the one-look one, not "no milestone tag")
    r2 = look(repo)
    assert r2.returncode and "one look" in r2.stderr
    # rewrite the ledger file (remove the rows): the git-history check still refuses
    p = repo / "results" / "SEALED_LEDGER.md"
    p.write_text("\n".join(l for l in p.read_text().splitlines() if "look=" not in l) + "\n")
    commit(repo, "rewrite ledger")
    sh("git", "tag", "milestone-t2", cwd=repo)
    r3 = look(repo)
    assert r3.returncode and "one look" in r3.stderr
    # retry after a success is refused
    r4 = look(repo, "--retry", lid)
    assert r4.returncode
    # the claim script reads the manifest hashes and the base commit from the ledger lines of the look
    sys.path.insert(0, ROOT)
    from bench.t1 import claims_a5 as CL
    lr, base = CL.ledger_look(str(repo), lid)
    assert set(lr) == {"T1", "T5"} and all(len(v) == 64 for v in lr.values())
    assert base == sh("git", "rev-parse", "milestone-t1", cwd=repo)


def test_missing_prereg_tag_refused(repo):
    sh("git", "tag", "-d", "prereg-v0.4.4", cwd=repo)
    r = look(repo)
    assert r.returncode and "prereg-v0.4.4" in r.stderr and "SEALED_LEDGER" not in sh("git", "ls-files", cwd=repo)


def test_no_milestone_refused(repo):
    sh("git", "tag", "-d", "milestone-t1", cwd=repo)
    r = look(repo)
    assert r.returncode and "milestone" in r.stderr


def test_hash_mismatch_refused(repo):
    write_cfg(repo, repo.tmp, repo.c9c, sha(repo / "results" / "c9" / "c9_power.json"), hashes={"dummy/cmd.py": "0" * 64})
    commit(repo, "bad hash")
    sh("git", "tag", "-f", "milestone-t1", cwd=repo)
    r = look(repo)
    assert r.returncode and "hash mismatch" in r.stderr


def test_c9_mismatch_refused(repo):
    write_cfg(repo, repo.tmp, repo.c9c, sha(repo / "results" / "c9" / "c9_power.json"), verdict="GO")
    commit(repo, "c9 says GO but file gives NO-GO")
    sh("git", "tag", "-f", "milestone-t1", cwd=repo)
    r = look(repo)
    assert r.returncode and "C9 record says GO" in r.stderr
    write_cfg(repo, repo.tmp, repo.c9c, "f" * 64)
    commit(repo, "c9 sha wrong")
    sh("git", "tag", "-f", "milestone-t1", cwd=repo)
    r = look(repo)
    assert r.returncode and "C9 record sha256" in r.stderr


def test_go_branch_needs_a6_before_time_box(repo):
    p = repo / "results" / "c9" / "c9_power.json"
    p.write_text(json.dumps({"summary": {SRM.C9_QUANTITY: {"P_v1_minimum": 1.0}},
                             "walk_cells": [{"P_v1_minimum_by_run_budget": {"200": 0.97}}]}))
    commit(repo, "c9 GO")
    c9c = sh("git", "rev-parse", "HEAD", cwd=repo)
    write_cfg(repo, repo.tmp, c9c, sha(p), verdict="GO")
    commit(repo, "cfg GO")
    sh("git", "tag", "-f", "milestone-t1", cwd=repo)
    ctl(repo)
    r = look(repo)
    if datetime.datetime.now(datetime.timezone.utc) < SRM.A6_TIME_BOX_END:
        assert r.returncode and "A6 status file" in r.stderr
        st = repo / SRM.A6_STATUS_REL
        st.parent.mkdir(parents=True, exist_ok=True)
        st.write_text(json.dumps({"status": "FAILED", "step5_commit": None, "gate_output_sha256": None}))
        commit(repo, "a6 status FAILED")
        sh("git", "tag", "-f", "milestone-t1", cwd=repo)
        r = look(repo)
        assert r.returncode and "prereg-v0.4.5" in r.stderr      # a gate-run status needs the A6 tag
        st.write_text(json.dumps({"status": "not tagged", "step5_commit": None, "gate_output_sha256": None}))
        commit(repo, "a6 status not tagged")
        sh("git", "tag", "-f", "milestone-t1", cwd=repo)
        r = look(repo)
    assert r.returncode == 0, r.stderr
    assert "c9=GO" in ledger(repo)


def test_leak_detected_fails_look_and_deletes(repo):
    ctl(repo, T1={"leak": True})
    r = look(repo)
    assert r.returncode == 1
    led = ledger(repo)
    assert "leak_test=fail (1 files)" in led and "FAILED" in led
    assert not (repo / "results" / "leaked.txt").exists()
    lid = json.load(open(repo.tmp / "ctl" / "seen_T1.json"))["look_id"]
    assert not (repo.tmp / "lookroot" / lid).exists()


def test_failure_summary_frames_only_and_retry_rules(repo):
    ctl(repo, T5={"crash": True})
    r = look(repo)
    assert r.returncode == 1
    lid = json.load(open(repo.tmp / "ctl" / "seen_T5.json"))["look_id"]
    led = ledger(repo)
    assert f"| T5 | result: look={lid} exit=1" in led and f"| T1 | result: look={lid} exit=0" in led  # both run
    summ = [f for f in os.listdir(repo / "results") if f.endswith("_failure.md")]
    assert len(summ) == 1
    txt = (repo / "results" / summ[0]).read_text()
    assert "ValueError" in txt and "cmd.py" in txt and "secret message" not in txt and "999" not in txt
    assert "12345" not in txt
    assert not (repo.tmp / "lookroot" / lid).exists()
    # retry without a deviation entry: refused
    sh("git", "tag", "milestone-fix1", cwd=repo)
    r = look(repo, "--retry", lid)
    assert r.returncode and "deviations" in r.stderr
    # deviation entry but the same milestone tag: refused
    write_cfg(repo, repo.tmp, repo.c9c, sha(repo / "results" / "c9" / "c9_power.json"),
              deviations=[{"date": "2026-10-28", "retry_of": lid, "fix": "x"}])
    commit(repo, "deviation")
    sh("git", "tag", "-f", "milestone-t1", cwd=repo)
    r = look(repo, "--retry", lid)
    assert r.returncode and "NEW milestone" in r.stderr
    sh("git", "tag", "-d", "milestone-t1", cwd=repo)
    sh("git", "tag", "-d", "milestone-fix1", cwd=repo)
    sh("git", "tag", "milestone-fix2", cwd=repo)
    ctl(repo, T1={"exit": 2})
    r = look(repo, "--retry", lid)                     # retry 1 (fails again)
    assert r.returncode == 1, r.stderr
    lid2 = json.load(open(repo.tmp / "ctl" / "seen_T1.json"))["look_id"]
    assert lid2 != lid and "attempt=2" in ledger(repo)
    # retrying a look that is not the latest: refused
    r = look(repo, "--retry", lid)
    assert r.returncode and "latest" in r.stderr
    write_cfg(repo, repo.tmp, repo.c9c, sha(repo / "results" / "c9" / "c9_power.json"),
              deviations=[{"date": "2026-10-28", "retry_of": lid, "fix": "x"},
                          {"date": "2026-10-28", "retry_of": lid2, "fix": "y"}])
    commit(repo, "deviation 2")
    sh("git", "tag", "-d", "milestone-fix2", cwd=repo)
    sh("git", "tag", "milestone-fix3", cwd=repo)
    ctl(repo, T5={"exit": 1})
    r = look(repo, "--retry", lid2)                    # retry 2 (fails)
    assert r.returncode == 1, r.stderr
    lid3 = json.load(open(repo.tmp / "ctl" / "seen_T5.json"))["look_id"]
    write_cfg(repo, repo.tmp, repo.c9c, sha(repo / "results" / "c9" / "c9_power.json"),
              deviations=[{"date": "2026-10-28", "retry_of": x, "fix": "z"} for x in (lid, lid2, lid3)])
    commit(repo, "deviation 3")
    sh("git", "tag", "-d", "milestone-fix3", cwd=repo)
    sh("git", "tag", "milestone-fix4", cwd=repo)
    r = look(repo, "--retry", lid3)                    # a third retry: refused
    assert r.returncode and "retries already used" in r.stderr


def _wait_for(p, t=60):
    t0 = time.time()
    while not p.exists():
        if time.time() - t0 > t:
            raise TimeoutError(p)
        time.sleep(0.2)


def test_sigterm_writes_failure_line(repo):
    ctl(repo, T5={"sleep": 60})
    e = dict(os.environ)
    pr = subprocess.Popen([sys.executable, SR, "--root", str(repo), "--look", "a5", "--dev"], env=e,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    _wait_for(repo.tmp / "ctl" / "seen_T5.json")
    time.sleep(0.5)
    pr.send_signal(signal.SIGTERM)
    pr.wait(timeout=60)
    assert pr.returncode == 143
    led = ledger(repo)
    lid = json.load(open(repo.tmp / "ctl" / "seen_T5.json"))["look_id"]
    assert f"| T5 | result: look={lid} failed (terminated)" in led and f"| T1 | result: look={lid} failed (terminated)" in led
    assert not os.access(repo / "data" / "sealed" / "X5" / "a.raw", os.R_OK)
    assert not (repo.tmp / "lookroot" / lid).exists()
    _, looks = SRM.look_state(str(repo))
    assert lid not in SRM.unterminated(looks) and not SRM.look_succeeded(looks[lid])


def test_unterminated_blocks_until_recover(repo):
    ctl(repo, T5={"sleep": 60})
    pr = subprocess.Popen([sys.executable, SR, "--root", str(repo), "--look", "a5", "--dev"],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, start_new_session=True)
    _wait_for(repo.tmp / "ctl" / "seen_T5.json")
    os.killpg(pr.pid, signal.SIGKILL)                  # node loss: no result line
    pr.wait(timeout=30)
    lid = json.load(open(repo.tmp / "ctl" / "seen_T5.json"))["look_id"]
    assert (repo.tmp / "lookroot" / lid).exists()
    for p in (repo / "data" / "sealed").rglob("*.raw"):
        os.chmod(p, 0o000)                             # (the epilogue would relock; the test does it by hand)
    sh("git", "tag", "milestone-fix", cwd=repo)
    r = look(repo)
    assert r.returncode and f"--recover {lid}" in r.stderr
    r = look(repo, "--retry", lid)
    assert r.returncode and "--recover" in r.stderr
    r = run(repo, "--recover", lid)
    assert r.returncode == 0, r.stderr
    assert f"| T5 | result: look={lid} failed (unterminated)" in ledger(repo)
    assert not (repo.tmp / "lookroot" / lid).exists()
    r = run(repo, "--recover", lid)                    # nothing left to recover
    assert r.returncode
    r = look(repo)                                     # a fresh look is still refused (one look)
    assert r.returncode and "one look" in r.stderr


def test_look_takes_no_command(repo):
    r = run(repo, "--look", "a5", "--dev", "--", "echo", "hi")
    assert r.returncode and "takes no command" in r.stderr
