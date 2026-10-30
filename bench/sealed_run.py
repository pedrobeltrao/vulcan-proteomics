#!/usr/bin/env python3
"""Run a command against SEALED data (DESIGN.md 6.6).  The only sanctioned way to read data/sealed/.

    python3 bench/sealed_run.py --dataset T2 --purpose "M3 headline ID/quant run" -- <command ...>
    python3 bench/sealed_run.py --look a5                      # the ONE joint T5 + T1 look (DESIGN §6.0d, A5.5)
    python3 bench/sealed_run.py --recover <look-id>            # close an unterminated look (A5.5.6)
    python3 bench/sealed_run.py --look a5 --retry <look-id>    # one of at most two retries (A5.5.6)

Plain runs (any dataset NOT in ONE_LOOK = {T1, T5}); rules, in order (any failure aborts before a ledger line is written
or a file is unlocked):
  1. the git working tree is clean (`git status --porcelain` empty; ignored files do not count);
  2. HEAD carries a tag matching `milestone-*`.  Ledger-only commits made by this script (subject starts with
     "sealed-ledger:", touching only results/SEALED_LEDGER.md) are skipped when looking for the tag, so several sealed
     runs per milestone are possible, each leaving its own ledger line;
  3. every sealed file of --dataset exists in data/sealed/ (listed from data/manifests/dev_sealed.tsv, role=sealed,
     set_id == dataset, non-truth rows) and, with --verify, matches data/manifests/sealed.sha256 (this reads the files);
  4. a line (UTC date, sha, tag, dataset, command, purpose) is appended to results/SEALED_LEDGER.md and committed;
  5. the files are chmod 0400, the command runs (no shell) with SEALED_FILES (space separated) and SEALED_DATASET in
     its environment and `{files}` args expanded, then the files are chmod 000 again (also on error, SIGINT, SIGTERM);
  6. a result line (exit code, wall time) is appended to the ledger and committed (best effort).
--dry-run checks 1-3 and prints what would happen, touching nothing.  --root overrides the project root (tests).

ONE LOOK for T1 and T5 (Amendment A5, DESIGN §6.0d; reviews/amendment_A5_draft.md A5.5 and A5.11.1, normative):
  * `--dataset T1` / `--dataset T5` are REFUSED; the only way to open them is `--look a5`, which takes NO command: the
    commands, purposes, child environment (an allow-list with explicit values; every other variable is dropped), look
    root, console whitelist and leak paths come from config/a5.yaml `look`.
  * Preconditions (A5.5.1), all checked before any ledger line: clean tree; a `milestone-*` tag on the base commit
    (after skipping ledger-only commits); `prereg-v0.4.3` and `prereg-v0.4.4` are ancestors of the base commit; the
    C9 record (config c9) matches the pilot file at its commit (sha256) and the recomputed GO/NO-GO; in the GO branch
    before the A6 time box ends (D31: 18:00 CET 2026-10-28) also `prereg-v0.4.5` and a committed A6 status file; the A3
    freeze (vulcan.gates.a3_config.freeze_status, tree list extended by A5.11.6) holds; every code-file sha256 and
    environment manifest digest of config `hashes` / `env_manifests` matches; the ledger -- the current file AND every
    line ever added in its git history (`git log -p`) -- has no T1 or T5 row (or the retry conditions hold); no look is
    unterminated; `--verify` is forced.
  * One commit writes both ledger lines (T5, then T1) with a shared look id BEFORE any file is unlocked.  T5 then T1
    run, each with its own fresh `<look.root>/<look-id>/<ds>/` (mode 0700) and with SEALED_LOOK_ID / SEALED_LOOK_ROOT
    exported.  Both always run.  Child stdout+stderr go to `<ds>/logs/console.log`; the terminal shows only lines that
    match the whitelist.  Each dataset's result line (exit, wall, sha256 of `<ds>/out/manifest.sha256`) is committed
    as soon as it ends.  Then the leak test (no file newer than the look start under the configured paths outside the
    look root and node-local scratch; the ledger excepted) writes its own result line.
  * A look FAILS if a command exits non-zero, the leak test fails, or a result line is missing.  A failed look writes a
    failure summary (dataset, last whitelisted stage line, exception type, traceback file:line:function frames; no
    messages, no locals) to results/<date>_look_<id>_failure.md, deletes the look root and node-local scratch unread,
    and commits.  SIGTERM / SIGINT write "failed (terminated)" result lines first.
  * `--recover <id>` closes an unterminated look ("failed (unterminated)"), deletes its root and scratch (reads only
    names, prints none) and commits.  `--retry <id>`: only for the latest look, only if it failed, never after a
    successful look, at most two retries, only with a dated `config/a5.yaml deviations` entry `retry_of: <id>` and a
    new milestone tag.
  * `--dev` (A5.11.8 dev dry run) runs the same code in a temporary clone (refused on the project root itself), where
    the harness provides fake sealed sets (dev raws) and fake tags.  It skips only the A3 freeze check.
"""
import argparse, csv, datetime, fnmatch, glob, hashlib, os, re, shlex, shutil, signal, subprocess, sys, time

LEDGER_REL = "results/SEALED_LEDGER.md"
HEADER = ("# SEALED_LEDGER\n\nEvery use of a sealed dataset (data/sealed/, DESIGN.md 6.6) is appended here by `bench/sealed_run.py` "
          "*before* the files are unlocked.\nDo not edit by hand.  A line without a following `result` line means the run was "
          "interrupted or is in progress.\n\n| date (UTC) | sha | tag | dataset | command | purpose |\n|---|---|---|---|---|---|\n")
ONE_LOOK = ("T5", "T1")                 # A5.11.1 item 1; also the run order (A5.5.2: T5 first, then T1)
REQUIRED_TAGS = ("prereg-v0.4.3", "prereg-v0.4.4")
A6_TAG = "prereg-v0.4.5"
A6_STATUS_REL = "results/2026-10-28_A6_hn1/a6_status.json"
A6_STATUS_ENUM = ("PASS", "FAILED", "not evaluable", "invalid pool", "not completed", "label-free stop", "not tagged",
                  "stopped: node time limit")
A6_TIME_BOX_END = datetime.datetime(2026, 10, 28, 17, 0, 0, tzinfo=datetime.timezone.utc)   # 18:00 CET (D31.3)
C9_GO_THRESHOLD = 0.8                   # A5.1.1
C9_QUANTITY = "f0.20/own/R3/kappa1.0/budget200"   # summary[...]["P_v1_minimum"] (config c9.quantity)
C9_BUDGET = "200"                       # and (stricter) the minimum over every walk cell at 200 runs
MAX_RETRIES = 2
CONFIG_REL = "config/a5.yaml"
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def die(msg):
    sys.exit(f"sealed_run: REFUSED, {msg}")


def git(root, *a, check=True):
    r = subprocess.run(["git", "-C", root, *a], capture_output=True, text=True)
    if check and r.returncode:
        sys.exit(f"sealed_run: git {' '.join(a)} failed: {r.stderr.strip()}")
    return r.stdout.strip()


def base_commit_and_tag(root):
    """Skip ledger-only commits from HEAD; return (sha, milestone tag) of the first other commit, tag None if absent."""
    rev = "HEAD"
    for _ in range(1000):
        subj = git(root, "log", "-1", "--format=%s", rev)
        files = git(root, "show", "--name-only", "--format=", rev).split()
        if subj.startswith("sealed-ledger:") and LEDGER_REL in files and all(
                f == LEDGER_REL or re.fullmatch(r"results/\d{4}-\d\d-\d\d_look_[0-9A-Za-z_-]+_failure\.md", f)
                for f in files):                # ledger lines (+ a failed look's frames-only summary, A5.5.6)
            rev += "~1"; continue
        tags = [t for t in git(root, "tag", "--points-at", rev).split() if fnmatch.fnmatch(t, "milestone-*")]
        return git(root, "rev-parse", rev), (sorted(tags)[-1] if tags else None)
    sys.exit("sealed_run: too many ledger-only commits")


def sealed_files(root, dataset):
    out = []
    with open(os.path.join(root, "data/manifests/dev_sealed.tsv"), newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["role"] == "sealed" and r["set_id"] == dataset and not r["id"].endswith("-truth"):
                out.append(os.path.join(root, "data/sealed", r["accession"], r["filename"]))
    if not out:
        sys.exit(f"sealed_run: no sealed files for dataset '{dataset}' in dev_sealed.tsv")
    missing = [p for p in out if not os.path.exists(p)]
    if missing:
        sys.exit("sealed_run: missing sealed files:\n  " + "\n  ".join(missing))
    return out


def verify(root, files):
    want = {}
    for l in open(os.path.join(root, "data/manifests/sealed.sha256")):
        if l.strip():
            h, p = l.rstrip("\n").split("  ", 1); want[os.path.join(root, "data", p)] = h
    for p in files:
        os.chmod(p, 0o400)
        try:
            h = hashlib.sha256()
            with open(p, "rb") as fh:
                for b in iter(lambda: fh.read(1 << 22), b""):
                    h.update(b)
        finally:
            os.chmod(p, 0o000)
        if want.get(p) != h.hexdigest():
            sys.exit(f"sealed_run: sha256 mismatch for {p}")
    print(f"sealed_run: sha256 verified for {len(files)} file(s)")


def ledger_commit(root, text, msg):
    path = os.path.join(root, LEDGER_REL)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if not os.path.exists(path):
        open(path, "w").write(HEADER)
    with open(path, "a") as fh:
        fh.write(text)
    git(root, "add", LEDGER_REL)
    git(root, "commit", "-q", "-m", "sealed-ledger: " + msg, "--", LEDGER_REL)


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def utcnow():
    return datetime.datetime.now(datetime.timezone.utc)


# ------------------------------------------------------------------------------------------------ ledger history
def _cells(line):
    s = line.strip()
    if not (s.startswith("|") and s.endswith("|")) or s.startswith("|---") or s.startswith("| date (UTC)"):
        return None
    return [c.strip() for c in re.split(r"(?<!\\)\|", s[1:-1])]


def ledger_lines(root):
    """Every ledger table row ever written: the current file plus every '+' line of `git log -p` on the ledger."""
    rows = []
    p = os.path.join(root, LEDGER_REL)
    if os.path.exists(p):
        rows += [l.rstrip("\n") for l in open(p)]
    hist = git(root, "log", "-p", "--format=", "--", LEDGER_REL, check=False)
    rows += [l[1:] for l in hist.splitlines() if l.startswith("+|")]
    seen, out = set(), []
    for l in rows:
        c = _cells(l)
        if c and len(c) >= 6 and l not in seen:
            seen.add(l); out.append(c)
    return out


_LID = re.compile(r"look=([0-9A-Za-z_-]+)")


def look_state(root):
    """-> (one_look_rows, looks).  looks: {id: {"datasets": [...], "results": {ds: text}, "leak": text|None,
    "tag": tag, "order": i}}.  one_look_rows: every opening row ever written for a ONE_LOOK dataset."""
    rows = ledger_lines(root)
    looks, one = {}, []
    for c in rows:
        ds, cmd, purpose = c[3], c[4], c[5]
        if cmd.startswith("result:"):
            m = _LID.search(cmd)
            if not m:
                continue
            L = looks.setdefault(m.group(1), {"datasets": [], "results": {}, "leak": None, "tag": None,
                                              "order": len(looks)})
            if ds == "look":
                L["leak"] = cmd
            else:
                L["results"][ds] = cmd
            continue
        if ds in ONE_LOOK:
            one.append(c)
            m = _LID.search(purpose)
            lid = m.group(1) if m else f"legacy-{len(one)}"
            L = looks.setdefault(lid, {"datasets": [], "results": {}, "leak": None, "tag": None, "order": len(looks)})
            L["datasets"].append(ds)
            L["tag"] = c[2]
    return one, looks


def look_succeeded(L):
    return (all(ds in L["results"] and "exit=0," in L["results"][ds] for ds in ONE_LOOK)
            and L["leak"] is not None and "leak_test=pass" in L["leak"])


def unterminated(looks):
    return [lid for lid, L in looks.items()
            if any(ds not in L["results"] for ds in L["datasets"]) or (L["datasets"] and L["leak"] is None)]


# ------------------------------------------------------------------------------------------------ config + checks
def load_config(root):
    try:
        import yaml
    except ImportError:
        die("PyYAML is needed to read config/a5.yaml (run under the project env)")
    p = os.path.join(root, CONFIG_REL)
    if not os.path.exists(p):
        die(f"{CONFIG_REL} missing")
    return yaml.safe_load(open(p))


def is_ancestor(root, tag, sha):
    if subprocess.run(["git", "-C", root, "rev-parse", "-q", "--verify", f"refs/tags/{tag}"],
                      capture_output=True).returncode:
        return False
    return subprocess.run(["git", "-C", root, "merge-base", "--is-ancestor", tag, sha]).returncode == 0


def git_show(root, rev, path):
    r = subprocess.run(["git", "-C", root, "show", f"{rev}:{path}"], capture_output=True)
    return r.stdout if r.returncode == 0 else None


def c9_branch(root, cfg, base):
    """A5.1.1: recompute GO / NO-GO from the committed pilot verdict file; refuse on any mismatch."""
    import json
    c9 = cfg.get("c9") or {}
    for k in ("verdict_file", "commit", "verdict_file_sha256", "verdict"):
        if not c9.get(k) or c9.get(k) == "TODO":
            die(f"config c9.{k} not recorded")
    if subprocess.run(["git", "-C", root, "merge-base", "--is-ancestor", str(c9["commit"]), base]).returncode:
        die(f"C9 commit {c9['commit']} is not an ancestor of the base commit")
    blob = git_show(root, c9["commit"], c9["verdict_file"])
    if blob is None:
        die(f"C9 file {c9['verdict_file']} not at commit {c9['commit']}")
    if hashlib.sha256(blob).hexdigest() != c9["verdict_file_sha256"]:
        die("C9 record sha256 does not match the pilot file")
    d = json.loads(blob)
    try:
        q = float(d["summary"][C9_QUANTITY]["P_v1_minimum"])
        cells = [float(cell["P_v1_minimum_by_run_budget"][C9_BUDGET]) for cell in d["walk_cells"]]
    except (KeyError, TypeError, ValueError):
        die("C9 pilot file lacks the recorded quantity")
    verdict = "GO" if q >= C9_GO_THRESHOLD and cells and min(cells) >= C9_GO_THRESHOLD else "NO-GO"
    if verdict != c9["verdict"]:
        die(f"C9 record says {c9['verdict']} but the pilot file gives {verdict}")
    return verdict


A6_TAGLESS = ("not tagged", "label-free stop", "stopped: node time limit", "not completed")   # no gate run needed a tag


def a6_ready(root, base, now=None):
    """GO branch: before the A6 time box ends the look needs a committed A6 status file in the closed enum, and, for a
    status produced by the gate run (PASS / FAILED / not evaluable / invalid pool), `prereg-v0.4.5` as an ancestor."""
    import json
    now = now or utcnow()
    if now >= A6_TIME_BOX_END:
        return "time box expired"
    blob = git_show(root, base, A6_STATUS_REL)
    try:
        st = json.loads(blob)["status"] if blob else None
    except Exception:
        st = None
    if st not in A6_STATUS_ENUM:
        die(f"GO branch: no committed A6 status file at {A6_STATUS_REL} before the A6 time box ends "
            f"({A6_TIME_BOX_END.isoformat()})")
    if st not in A6_TAGLESS and not is_ancestor(root, A6_TAG, base):
        die(f"GO branch: A6 status {st!r} but {A6_TAG} is not an ancestor of the base commit")
    return st


def check_hashes(root, cfg, dev):
    sys.dont_write_bytecode = True
    sys.path.insert(0, os.path.join(root, "bench"))
    import env_manifest as EM
    bad = []
    hashes = cfg.get("hashes") or {}
    if not isinstance(hashes, dict) or not hashes:
        die("config hashes missing")
    for rel, h in sorted(hashes.items()):
        p = os.path.join(root, rel)
        if not os.path.isfile(p) or sha256_file(p) != h:
            bad.append(rel)
    for name, e in sorted((cfg.get("env_manifests") or {}).items()):
        if EM.digest(e["path"]) != e["digest"]:
            bad.append(f"env:{name}")
    if not dev and not cfg.get("env_manifests"):
        bad.append("env_manifests missing")
    if bad:
        die("hash mismatch: " + ", ".join(bad))


def check_freeze(root):
    code = ("import json,sys; from vulcan.gates import a3_config as A; s=A.freeze_status(); "
            "print(json.dumps({'ok': bool(s.get('ok')), 'reason': s.get('reason')}))")
    r = subprocess.run([sys.executable, "-c", code], cwd=root, capture_output=True, text=True,
                       env={**os.environ, "PYTHONPATH": root})
    import json
    try:
        st = json.loads(r.stdout.strip().splitlines()[-1])
    except Exception:
        die(f"A3 freeze check failed to run: {r.stderr.strip()[-300:]}")
    if not st["ok"]:
        die(f"A3 freeze does not hold: {st.get('reason')}")


# ------------------------------------------------------------------------------------------------ look mechanics
def leak_test(root, cfg, t0, look_dir, scratch):
    """Files newer than t0 under the configured paths, outside the look root and node-local scratch; ledger excepted."""
    excl = [os.path.realpath(look_dir), os.path.realpath(scratch), os.path.realpath(os.path.join(root, LEDGER_REL))]
    found = []
    for pat in cfg["look"].get("leak_paths", []):
        pat = pat if os.path.isabs(pat) else os.path.join(root, pat)
        for top in glob.glob(pat):
            cands = [top] if not os.path.isdir(top) else (os.path.join(dp, f) for dp, _, fns in os.walk(top) for f in fns)
            for p in cands:
                rp = os.path.realpath(p)
                if any(rp == e or rp.startswith(e + os.sep) for e in excl):
                    continue
                try:
                    if os.lstat(p).st_mtime > t0:
                        found.append(p)
                except OSError:
                    pass
    return found


_FRAME = re.compile(r'^\s*File "([^"]+)", line (\d+), in (\S+)')
_EXC = re.compile(r"^([A-Za-z_][\w.]*(?:Error|Exception|Exit|Interrupt|Warning))\b")


def failure_summary(log_path, whitelist):
    """Dataset log -> frames only: last whitelisted stage line, exception type(s), traceback file:line:function."""
    out = {"stage": None, "frames": [], "exception_types": []}
    if not os.path.exists(log_path):
        return out
    for line in open(log_path, errors="replace"):
        line = line.rstrip("\n")
        if any(w.search(line) for w in whitelist):
            out["stage"] = line
        m = _FRAME.match(line)
        if m:
            f = m.group(1)
            f = os.path.relpath(f, PROJECT_ROOT) if f.startswith(PROJECT_ROOT) else os.path.basename(f)
            out["frames"].append(f"{f}:{m.group(2)}:{m.group(3)}")
        m = _EXC.match(line)
        if m:
            out["exception_types"].append(m.group(1))
    out["frames"] = out["frames"][-40:]
    out["exception_types"] = out["exception_types"][-5:]
    return out


def expand(cmd, ds_root, files, scratch=""):
    out = []
    for c in cmd:
        if c == "{files}":
            out += files
        else:
            out.append(c.replace("{root}", ds_root).replace("{scratch}", scratch))
    return out


def rm_tree(p):
    if os.path.lexists(p):
        for dp, dns, fns in os.walk(p):
            for d in dns:
                try: os.chmod(os.path.join(dp, d), 0o700)
                except OSError: pass
        shutil.rmtree(p, ignore_errors=True)


def real_sealed_present(root):
    """True if any sealed file listed in data/manifests/dev_sealed.tsv exists under root/data/sealed (the project)."""
    p = os.path.join(root, "data/manifests/dev_sealed.tsv")
    if not os.path.exists(p):
        return False
    with open(p, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("role") == "sealed" and os.path.lexists(os.path.join(root, "data/sealed", r["accession"],
                                                                           r["filename"])):
                return True
    return False


def run_look(a, root):
    dev = a.dev
    if dev and real_sealed_present(root):
        die("--dev runs only in a temporary clone (--root) without sealed data, never on the project root")
    if a.command:
        die("--look takes no command; the commands come from config/a5.yaml look.commands")
    if git(root, "status", "--porcelain"):
        die("git tree is not clean (commit or stash first)")
    cfg = load_config(root)
    look = cfg.get("look") or {}
    for k in ("commands", "purpose", "env", "root", "console_whitelist", "leak_paths", "time_limit"):
        if k not in look or look[k] in (None, "TODO"):
            die(f"config look.{k} missing")
    sha, tag = base_commit_and_tag(root)
    if not tag:
        die("the base commit has no milestone-* tag")
    for t in REQUIRED_TAGS:
        if not is_ancestor(root, t, sha):
            die(f"tag {t} is missing or not an ancestor of the base commit")
    branch = c9_branch(root, cfg, sha)
    a6 = a6_ready(root, sha) if branch == "GO" else "n/a (NO-GO)"
    if not dev:
        check_freeze(root)
    check_hashes(root, cfg, dev)
    one, looks = look_state(root)
    open_ = unterminated(looks)
    if open_:
        die(f"unterminated look {open_[0]}: run `sealed_run.py --recover {open_[0]}` first")
    attempt = 1
    if a.retry:
        if a.retry not in looks:
            die(f"--retry {a.retry}: no such look in the ledger history")
        latest = max(looks, key=lambda k: looks[k]["order"])
        if a.retry != latest:
            die(f"--retry {a.retry}: only the latest look ({latest}) can be retried")
        if any(look_succeeded(L) for L in looks.values()):
            die("--retry: an earlier look succeeded (two exit-0 result lines and a passed leak test)")
        if len(looks) - 1 >= MAX_RETRIES:
            die(f"--retry: {MAX_RETRIES} retries already used")
        devs = [d for d in (cfg.get("deviations") or []) if isinstance(d, dict) and d.get("retry_of") == a.retry
                and d.get("date")]
        if not devs:
            die(f"--retry: config/a5.yaml deviations has no dated entry with retry_of: {a.retry}")
        if tag in {L["tag"] for L in looks.values()}:
            die("--retry: needs a NEW milestone tag on the fix commit")
        attempt = len(looks) + 1
    elif one:
        die("one look (§6.0d A5.5.2): the ledger history already has a T1/T5 line; only --retry is possible")
    need = float(look.get("min_scratch_free_gb") or 0)
    free = shutil.disk_usage("/tmp").free / 2**30
    if free < need:
        die(f"node-local /tmp has {free:.0f} GB free < look.min_scratch_free_gb {need:.0f} (free space first)")
    files = {ds: (dev_files(root, cfg, ds) if dev else sealed_files(root, ds)) for ds in ONE_LOOK}
    for ds in ONE_LOOK:
        verify(root, files[ds])                 # forced (A5.5.1.6)
    lid = utcnow().strftime("%Y%m%dT%H%M%SZ") + "-" + sha[:8]
    look_dir = os.path.join(look["root"], lid)
    scratch = f"/tmp/vulcan_look_{lid}"
    if os.path.exists(look_dir):
        die(f"look root {look_dir} exists")
    os.makedirs(look["root"], exist_ok=True)
    os.mkdir(look_dir, 0o700)
    os.chmod(look_dir, 0o700)
    for ds in ONE_LOOK:
        os.makedirs(os.path.join(look_dir, ds, "logs"), mode=0o700)
    rm_tree(scratch)
    os.makedirs(scratch, mode=0o700)
    wl = [re.compile(w) for w in look["console_whitelist"]]
    esc = lambda s: str(s).replace("|", "\\|").replace("\n", " ")
    now = utcnow().strftime("%Y-%m-%d %H:%M:%SZ")
    text = "".join(f"| {now} | {sha[:12]} | {tag} | {ds} | `{esc(shlex.join(look['commands'][ds]))}` | "
                   f"look={lid} attempt={attempt} c9={branch} a6={esc(a6)}; {esc(look['purpose'][ds])} |\n"
                   for ds in ONE_LOOK)
    t0 = time.time()
    ps_log(look_dir, "start")
    ledger_commit(root, text, f"look a5 {lid} opened (T5, T1) at {tag}")
    print(f"sealed_run: look {lid} opened at {tag} @ {sha[:12]} (c9 {branch}); T5 then T1", flush=True)
    state = {"done": set(), "rc": {}, "proc": None}

    def relock_all():
        for fs in files.values():
            for p in fs:
                try: os.chmod(p, 0o000)
                except OSError: pass

    def on_signal(signum, _frame):
        if state["proc"] is not None and state["proc"].poll() is None:
            try: state["proc"].terminate()
            except OSError: pass
        relock_all()
        why = "terminated" if signum == signal.SIGTERM else "interrupted"
        txt = "".join(f"|  |  |  | {ds} | result: look={lid} failed ({why}) |  |\n" for ds in ONE_LOOK
                      if ds not in state["done"])
        txt += f"|  |  |  | look | result: look={lid} leak_test=not run ({why}) |  |\n"
        try:
            if state["proc"] is not None:
                state["proc"].wait(timeout=20)
        except Exception:
            pass
        try:
            if state.get("log") is not None:
                state["log"].close()
        except Exception:
            pass
        for _ in range(3):                              # NFS: open handles leave .nfs* files for a moment
            rm_tree(look_dir); rm_tree(scratch)         # a failed look's outputs are deleted unread (A5.5.6)
            if not os.path.exists(look_dir):
                break
            time.sleep(1)
        txt += f"|  |  |  | look | result: look={lid} FAILED ({why}), outputs deleted unread |  |\n"
        try:
            ledger_commit(root, txt, f"look {lid} failed ({why})")
        except SystemExit:
            pass
        os._exit(143 if signum == signal.SIGTERM else 130)
    signal.signal(signal.SIGTERM, on_signal)
    signal.signal(signal.SIGINT, on_signal)

    for ds in ONE_LOOK:
        ds_root = os.path.join(look_dir, ds)
        env = {str(k): str(v) for k, v in look["env"].items()}
        env.update(SEALED_FILES=" ".join(files[ds]), SEALED_DATASET=ds, SEALED_LOOK_ID=lid, SEALED_LOOK_ROOT=ds_root,
                   SEALED_LOOK_SCRATCH=scratch)
        cmd = expand(look["commands"][ds], ds_root, files[ds], scratch)
        rc, t1 = 1, time.time()
        print(f"sealed_run: {ds} start", flush=True)
        try:
            for p in files[ds]:
                os.chmod(p, 0o400)
            with open(os.path.join(ds_root, "logs", "console.log"), "w") as lf:
                state["log"] = lf
                pr = subprocess.Popen(cmd, env=env, cwd=root, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                      text=True, errors="replace", bufsize=1)
                state["proc"] = pr
                for line in pr.stdout:
                    lf.write(line); lf.flush()
                    if any(w.search(line.rstrip("\n")) for w in wl):
                        print(f"  {ds}| {line.rstrip()}", flush=True)
                rc = pr.wait()
        except OSError as e:
            print(f"sealed_run: {ds} could not start ({type(e).__name__})", flush=True)
            rc = 127
        finally:
            for p in files[ds]:
                try: os.chmod(p, 0o000)
                except OSError: pass
        man = os.path.join(ds_root, "out", "manifest.sha256")
        mh = sha256_file(man) if os.path.exists(man) else "none"
        state["rc"][ds] = rc
        ledger_commit(root, f"|  |  |  | {ds} | result: look={lid} exit={rc}, wall={time.time() - t1:.0f}s, "
                            f"manifest_sha256={mh} |  |\n", f"look {lid} {ds} result exit={rc}")
        state["done"].add(ds)
        print(f"sealed_run: {ds} exit={rc}", flush=True)
    leaks = leak_test(root, cfg, t0, look_dir, scratch)
    with open(os.path.join(look_dir, "leak_test.txt"), "w") as fh:
        fh.write("\n".join(leaks) + ("\n" if leaks else ""))
    ok_leak = not leaks
    ps_log(look_dir, "end")
    failed = (not ok_leak) or any(state["rc"][ds] != 0 for ds in ONE_LOOK)
    summ_rel = None
    if failed:
        import json
        summ = {ds: dict(exit=state["rc"][ds], **failure_summary(os.path.join(look_dir, ds, "logs", "console.log"), wl))
                for ds in ONE_LOOK}
        summ["leak_test"] = {"pass": ok_leak, "n_files": len(leaks)}
        summ_rel = f"results/{datetime.date.today().isoformat()}_look_{lid}_failure.md"
        os.makedirs(os.path.join(root, os.path.dirname(summ_rel)), exist_ok=True)
        with open(os.path.join(root, summ_rel), "w") as fh:
            fh.write(f"# Failed A5 look {lid} (frames only; A5.5.6)\n\n```json\n{json.dumps(summ, indent=1)}\n```\n")
        rm_tree(look_dir); rm_tree(scratch)
        for p in leaks:
            try: os.remove(p)
            except OSError: pass
        git(root, "add", summ_rel)
    path = os.path.join(root, LEDGER_REL)
    with open(path, "a") as fh:
        fh.write(f"|  |  |  | look | result: look={lid} leak_test={'pass' if ok_leak else f'fail ({len(leaks)} files)'}"
                 f"{'; FAILED, outputs deleted unread, summary ' + summ_rel if failed else ''} |  |\n")
    git(root, "add", LEDGER_REL)
    git(root, "commit", "-q", "-m", f"sealed-ledger: look {lid} {'FAILED' if failed else 'complete'}", "--",
        LEDGER_REL, *([summ_rel] if summ_rel else []))
    rm_tree(scratch)                                    # node-local scratch (replicate indexes are not kept)
    print(f"sealed_run: look {lid} {'FAILED' if failed else 'complete'}; leak test {'pass' if ok_leak else 'FAIL'}"
          f"; outputs {'deleted' if failed else 'in ' + look_dir}", flush=True)
    return 1 if failed else 0


def ps_log(look_dir, when):
    try:
        r = subprocess.run(["ps", "-eo", "user,pid,etime,pcpu,args"], capture_output=True, text=True)
        open(os.path.join(look_dir, f"ps_{when}.txt"), "w").write(r.stdout)
    except OSError:
        pass


def dev_files(root, cfg, ds):
    """--dev: the stand-in files of the A5.11.8 dry run (config look.dev_files, paths relative to the clone)."""
    fs = [os.path.join(root, p) for p in ((cfg["look"].get("dev_files") or {}).get(ds) or [])]
    if not fs or not all(os.path.exists(p) for p in fs):
        die(f"--dev: look.dev_files.{ds} missing")
    return fs


def recover(a, root):
    if git(root, "status", "--porcelain"):
        die("git tree is not clean")
    cfg = load_config(root)
    _, looks = look_state(root)
    lid = a.recover
    if lid not in looks or lid not in unterminated(looks):
        die(f"--recover {lid}: not an unterminated look")
    L = looks[lid]
    txt = "".join(f"|  |  |  | {ds} | result: look={lid} failed (unterminated) |  |\n" for ds in ONE_LOOK
                  if ds in L["datasets"] and ds not in L["results"])
    if L["leak"] is None:
        txt += f"|  |  |  | look | result: look={lid} leak_test=not run (unterminated); outputs deleted unread |  |\n"
    rm_tree(os.path.join((cfg.get("look") or {}).get("root", "/nonexistent"), lid))
    rm_tree(f"/tmp/vulcan_look_{lid}")
    ledger_commit(root, txt, f"look {lid} recovered (failed, unterminated)")
    print(f"sealed_run: look {lid} recorded as failed (unterminated); root deleted")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", help="sealed set id from DESIGN 6.6 (T2, T3, T4); T1/T5 only via --look a5")
    ap.add_argument("--purpose")
    ap.add_argument("--verify", action="store_true", help="check sha256 of the sealed files against sealed.sha256 first")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--look", choices=["a5"], help="the one joint T5 + T1 look (DESIGN §6.0d)")
    ap.add_argument("--retry", metavar="LOOK_ID")
    ap.add_argument("--recover", metavar="LOOK_ID")
    ap.add_argument("--dev", action="store_true", help="A5.11.8 dev dry run in a temporary clone (--root)")
    ap.add_argument("--root", default=PROJECT_ROOT)
    ap.add_argument("command", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    a.command = a.command[1:] if a.command[:1] == ["--"] else a.command
    root = a.root
    if a.recover:
        sys.exit(recover(a, root))
    if a.look:
        if a.dataset or a.dry_run:
            die("--look takes neither --dataset nor --dry-run")
        sys.exit(run_look(a, root))
    if a.retry or a.dev:
        die("--retry / --dev only with --look a5")
    if not a.dataset or not a.purpose:
        die("give --dataset and --purpose (or --look a5)")
    if a.dataset in ONE_LOOK:
        die(f"{a.dataset} is opened only by the one look: `sealed_run.py --look a5` (DESIGN §6.0d A5.11.1)")
    cmd = a.command
    if not cmd:
        sys.exit("sealed_run: give the command after `--`")
    if git(root, "status", "--porcelain"):
        sys.exit("sealed_run: REFUSED, git tree is not clean (commit or stash first)")
    sha, tag = base_commit_and_tag(root)
    if not tag:
        sys.exit("sealed_run: REFUSED, HEAD has no milestone-* tag (tag the commit you want evaluated)")
    files = sealed_files(root, a.dataset)
    if a.verify:
        verify(root, files)
    cmd = [x for c in cmd for x in (files if c == "{files}" else [c])]
    if a.dry_run:
        print(f"sealed_run (dry run): tree clean, tag {tag} @ {sha[:12]}, {len(files)} file(s) of {a.dataset}; would run: {shlex.join(cmd)}")
        return
    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ")
    esc = lambda s: s.replace("|", "\\|").replace("\n", " ")
    ledger_commit(root, f"| {now} | {sha[:12]} | {tag} | {a.dataset} | `{esc(shlex.join(cmd))}` | {esc(a.purpose)} |\n", f"{a.dataset} at {tag}")

    def relock(*_):
        for p in files:
            try: os.chmod(p, 0o000)
            except OSError: pass
    signal.signal(signal.SIGTERM, lambda *_: (relock(), sys.exit(143)))
    rc, t0 = 1, time.time()
    try:
        for p in files:
            os.chmod(p, 0o400)
        env = dict(os.environ, SEALED_FILES=" ".join(files), SEALED_DATASET=a.dataset)
        rc = subprocess.run(cmd, env=env, cwd=root).returncode
    finally:
        relock()
    try:
        ledger_commit(root, f"|  |  |  | {a.dataset} | result: exit={rc}, wall={time.time() - t0:.0f}s |  |\n", f"{a.dataset} result exit={rc}")
    except SystemExit:
        print("sealed_run: warning: could not commit the result line", file=sys.stderr)
    sys.exit(rc)


if __name__ == "__main__":
    main()
