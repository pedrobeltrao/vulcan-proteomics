#!/usr/bin/env python3
"""A5.11.8 dev `--look` dry run (DESIGN §6.0d): the whole look mechanics on dev data, in a temporary CLONE (never the
project ledger, never sealed data).

    python scripts/a5/dev_look_dryrun.py --commit <sha> --base <dir outside the repo> [--go]

What it does:
  1. `git clone --shared` of the repository at <commit> into <base>/clone (tags created there stay there).  Ignored
     resources (conda_envs, tools, data/raw, data/fasta, data/cache) are symlinked from the project.
  2. Stand-in "sealed" sets: small placeholder files under clone data/sealed/DEV/ (one per dev run, named after it),
     listed in the clone's data/manifests/sealed.sha256 and in config look.dev_files.  sealed_run locks/unlocks and
     sha256-verifies them exactly as it would T1/T5 files.  The commands read the real dev raws: T5 = D7 (LM2486/2487/
     2490, entrapment seeds 2001-2003, decoy seed 1001, as the dev arm-M design) via bench/t5_armM.py --mode dev; T1 =
     D1 via bench/t1_driver.py --dataset D1 --a3 (all look-mode hardening active: SEALED_LOOK_ID is set).
  3. The clone's config/a5.yaml look section is the production one except: commands (dev data), root (<base>/lookroot),
     dev_files, leak_paths (the clone's results/ and work/ plus the production /tmp globs), VULCAN_LOCAL_ENV (the
     existing staged env, to save the copy time; the production look rebuilds it).  `hashes` are recomputed in the
     clone for every listed code file; env manifests are the real ones.  Fake tags prereg-v0.4.3/4/5 and
     milestone-a5dev, and (GO branch before the A6 time box) a fake A6 status file "not completed", exist in the clone
     only.
  4. Runs `bench/sealed_run.py --root <clone> --look a5 --dev`, times it, then runs the claim script on both outputs
     (dev mode: no ledger sha check is skipped; the ledger result lines of the clone are used) and the wording lint.
Writes <out>/dev_look.json (wall times, exit codes, ledger lines, console transcript, leak test, claim-script exit,
lint result).  The v1 bit-identity check (v1 columns with and without the preceding arm-M step) is done by
compare_v1 against a plain `--a3` D1 run of the same commit (the A3.1.6 dry run's "with" run).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PY = os.path.join(ROOT, "conda_envs", "pe", "bin", "python")
D7_RUNS = ["LM2486_DDA", "LM2487_DDA", "LM2490_DDA"]
D1_RUNS = ["Velos005137"]
LINKS = ["conda_envs", "tools", "data/raw", "data/fasta", "data/cache"]
HASHED = ["bench/sealed_run.py", "bench/t1_driver.py", "bench/t5_armM.py", "bench/env_manifest.py"]


def sh(*a, cwd=None, check=True, env=None):
    r = subprocess.run(list(a), cwd=cwd, capture_output=True, text=True, env=env)
    if check and r.returncode:
        raise RuntimeError(f"{a}: {r.stderr[-2000:]}")
    return r.stdout.strip()


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def compare_v1(m_look: dict, m_ref: dict) -> dict:
    """v1 columns: methods.vulcan, noninferiority (numerator vulcan), matched_fdp_curves.vulcan -- exact equality."""
    out = {}
    for k, f in (("methods.vulcan", lambda m: m["methods"]["vulcan"]),
                 ("noninferiority", lambda m: m["noninferiority"]),
                 ("matched_fdp_curves.vulcan", lambda m: m["matched_fdp_curves"]["vulcan"])):
        a, b = f(m_look), f(m_ref)
        out[k] = json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--commit", default="HEAD")
    ap.add_argument("--base", required=True)
    ap.add_argument("--out", required=True, help="where dev_look.json goes")
    ap.add_argument("--ref-metrics", default=None, help="metrics.json of a plain --a3 D1 run of the same commit")
    ap.add_argument("--t1-extra", nargs="*", default=[])
    a = ap.parse_args()
    import yaml
    commit = sh("git", "-C", ROOT, "rev-parse", a.commit)
    base = os.path.abspath(a.base)
    assert not base.startswith(ROOT + os.sep), "base must be outside the repository"
    clone = os.path.join(base, "clone")
    if os.path.exists(base):
        shutil.rmtree(base)
    os.makedirs(base)
    sh("git", "clone", "-q", "--shared", "--no-checkout", ROOT, clone)
    sh("git", "checkout", "-q", commit, cwd=clone)
    sh("git", "checkout", "-q", "-b", "a5dev", cwd=clone)
    sh("git", "config", "user.email", "a5dev@local", cwd=clone)
    sh("git", "config", "user.name", "a5 dev dry run", cwd=clone)
    for t in sh("git", "tag", cwd=clone).split():
        sh("git", "tag", "-d", t, cwd=clone)
    with open(os.path.join(clone, ".git", "info", "exclude"), "a") as fh:
        fh.write("".join(f"/{l}\n" for l in LINKS))
    for l in LINKS:
        d = os.path.join(clone, l)
        os.makedirs(os.path.dirname(d), exist_ok=True)
        if not os.path.lexists(d):
            os.symlink(os.path.join(ROOT, l), d)
    # stand-in sealed files
    dev_files, lines = {}, []
    for ds, runs in (("T5", D7_RUNS), ("T1", D1_RUNS)):
        dev_files[ds] = []
        for r in runs:
            rel = f"data/sealed/DEV/{ds}_{r}.raw"
            p = os.path.join(clone, rel)
            os.makedirs(os.path.dirname(p), exist_ok=True)
            open(p, "w").write(f"stand-in for dev run {r}\n")
            dev_files[ds].append(rel)
            lines.append(f"{sha(p)}  {rel[len('data/'):]}")
    open(os.path.join(clone, "data/manifests/sealed.sha256"), "w").write("\n".join(lines) + "\n")
    cfgp = os.path.join(clone, "config/a5.yaml")
    cfg = yaml.safe_load(open(cfgp))
    lk = cfg["look"]
    d7 = [os.path.join(ROOT, "data/raw/PXD031616", r + ".raw") for r in D7_RUNS]
    lk["commands"] = {
        "T5": [PY, "bench/t5_armM.py", "--mode", "dev", "--runs", *d7, "--entrap-seeds", "2001", "2002", "2003",
               "--decoy-seeds", "1001", "--out", "{root}/out", "--work", "{scratch}/T5/work"],
        "T1": [PY, "bench/t1_driver.py", "--dataset", "D1", "--a3", "--out", "{root}/out", "--work",
               "{scratch}/T1/work", "--rt-off-control", *a.t1_extra],
    }
    lk["root"] = os.path.join(base, "lookroot")
    lk["dev_files"] = dev_files
    lk["leak_paths"] = ["results", "work", "/tmp/vulcan_t1_*", "/tmp/vulcan_m3e*"]
    staged = f"/tmp/vulcan_env_{os.getuid()}"
    lk["env"]["VULCAN_LOCAL_ENV"] = staged
    lk["env"]["VULCAN_TRFP_BIN"] = os.path.join(staged, "trfp")
    lk["env"]["VULCAN_RAW_DLL"] = os.path.join(staged, "vulcanraw", "VulcanRaw.dll")
    lk["env"]["VULCAN_DOTNET_ROOT"] = os.path.join(staged, "dotnet")
    lk["env"]["PROJECT_ROOT"] = clone
    if lk.get("time_limit") in (None, "TODO"):
        lk["time_limit"] = "dev"
    hashed = list(HASHED) + sorted(os.path.relpath(p, clone) for p in __import__("glob").glob(
        os.path.join(clone, "bench/t1/*.py")))
    cfg["hashes"] = {h: sha(os.path.join(clone, h)) for h in hashed if os.path.exists(os.path.join(clone, h))}
    if not isinstance(cfg.get("env_manifests"), dict):
        cfg["env_manifests"] = {}
    yaml.safe_dump(cfg, open(cfgp, "w"), sort_keys=False, width=200)
    st = os.path.join(clone, "results/2026-10-28_A6_hn1/a6_status.json")
    os.makedirs(os.path.dirname(st), exist_ok=True)
    if not os.path.exists(st):
        json.dump({"status": "not completed", "step5_commit": None, "gate_output_sha256": None,
                   "note": "FAKE (dev look dry run clone only)"}, open(st, "w"))
    sh("git", "add", "-A", "--", "config", "data/manifests/sealed.sha256", "results/2026-10-28_A6_hn1", cwd=clone)
    sh("git", "commit", "-q", "-m", "a5 dev look dry run: dev config, stand-in sealed manifest, fake A6 status", cwd=clone)
    for t in ("prereg-v0.4.3", "prereg-v0.4.4", "prereg-v0.4.5", "milestone-a5dev"):
        sh("git", "tag", t, cwd=clone)
    assert sh("git", "status", "--porcelain", cwd=clone) == "", sh("git", "status", "--porcelain", cwd=clone)
    env = {k: v for k, v in os.environ.items() if not k.startswith("VULCAN_")}
    env["VULCAN_SAGE"] = "/caller/set/must/be/dropped"
    t0 = time.time()
    p = subprocess.run([PY, os.path.join(clone, "bench/sealed_run.py"), "--root", clone, "--look", "a5", "--dev"],
                       capture_output=True, text=True, env=env)
    wall = time.time() - t0
    led = open(os.path.join(clone, "results/SEALED_LEDGER.md")).read()
    rec = {"commit": commit, "wall_s": wall, "sealed_run_exit": p.returncode, "console_stdout": p.stdout,
           "console_stderr": p.stderr, "ledger": led}
    lid = None
    for l in led.splitlines():
        if "look=" in l and "| T5 |" in l and "result" not in l:
            lid = l.split("look=")[1].split()[0]
    rec["look_id"] = lid
    lroot = os.path.join(base, "lookroot", lid) if lid else None
    rec["look_root_exists"] = bool(lroot and os.path.isdir(lroot))
    if rec["look_root_exists"]:
        lt = os.path.join(lroot, "leak_test.txt")
        rec["leak_test_files"] = open(lt).read().splitlines() if os.path.exists(lt) else None
        # ledger result manifest hashes
        res = {}
        for l in led.splitlines():
            if f"result: look={lid}" in l and "manifest_sha256=" in l:
                ds = l.split("|")[4].strip()
                res[ds] = l.split("manifest_sha256=")[1].split()[0]
        rec["ledger_manifest_sha256"] = res
        for ds in ("T5", "T1"):
            man = os.path.join(lroot, ds, "out", "manifest.sha256")
            rec[f"{ds}_manifest_matches_ledger"] = os.path.exists(man) and sha(man) == res.get(ds)
            rec[f"{ds}_work_left"] = sorted(os.listdir(os.path.join(f"/tmp/vulcan_look_{lid}", ds, "work"))) \
                if os.path.isdir(os.path.join(f"/tmp/vulcan_look_{lid}", ds, "work")) else None
        # claim script, once, on the dev outputs
        cl_out = os.path.join(base, "claims")
        cp = subprocess.run([PY, "bench/t1/claims_a5.py", "--t1", os.path.join(lroot, "T1", "out", "metrics.json"),
                             "--t5", os.path.join(lroot, "T5", "out", "armM.json"), "--config", "config/a5.yaml",
                             "--out", cl_out, "--root", clone, "--look-id", lid],
                            cwd=clone, capture_output=True, text=True, env={**env, "PYTHONPATH": clone})
        rec["claims_exit"] = cp.returncode
        rec["claims_stderr_tail"] = cp.stderr[-3000:]
        rec["claims_files"] = sorted(os.listdir(cl_out)) if os.path.isdir(cl_out) else []
        if a.ref_metrics and os.path.exists(os.path.join(lroot, "T1", "out", "metrics.json")):
            rec["v1_bit_identical_with_and_without_armM_step"] = compare_v1(
                json.load(open(os.path.join(lroot, "T1", "out", "metrics.json"))), json.load(open(a.ref_metrics)))
    os.makedirs(a.out, exist_ok=True)
    json.dump(rec, open(os.path.join(a.out, "dev_look.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in rec.items() if k not in ("ledger", "console_stdout")}, indent=1)[:4000])


if __name__ == "__main__":
    main()
