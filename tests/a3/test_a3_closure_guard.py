"""scripts/a3/v1_closure_guard.py: the transitive import closure of the v1 path and the allow-list rule."""
import v1_closure_guard as V


def test_imports_of_handles_function_level_and_relative():
    src = "import a.b\nfrom .x import y\ndef f():\n    from c import d\n    import e\n"
    names = V.imports_of(src, "pkg/sub/mod.py")
    assert {"a.b", "pkg.sub.x", "pkg.sub.x.y", "c", "c.d", "e"} <= set(names)


def test_closure_at_tag_contains_the_v1_path():
    cl = V.closure(V.TAG)
    for f in ("scripts/a2/gates_a2.py", "vulcan/gates/a2_stats.py", "vulcan/gates/stats.py", "vulcan/stats/g2.py",
              "vulcan/stats/tdc.py", "vulcan/gates/homolog.py", "vulcan/rescore/static.py", "vulcan/gpu/engine.py",
              "vulcan/index/build.py", "scripts/m3e/gates_a1.py", "scripts/m3e/seqclf.py", "scripts/m2/e2_correctness.py",
              "vulcan/gpu/kernels.cu", "config/a2.yaml", "bench/t1_driver.py", "bench/t1/evaluate.py"):
        assert f in cl, f
    assert not any(f.startswith("scripts/a3/") for f in cl)


def test_gates_a2_path_passes_and_unallowed_change_fails(monkeypatch):
    rep = V.check("HEAD")
    assert rep["by_entry"]["scripts/a2/gates_a2.py"]["pass"], rep["by_entry"]["scripts/a2/gates_a2.py"]
    real = V.diff_sha

    def fake(commit, path):
        sha, n = real(commit, path)
        return ("f" * 64, 1) if path == "vulcan/gpu/engine.py" else (sha, n)
    monkeypatch.setattr(V, "diff_sha", fake)
    rep = V.check("HEAD")
    assert "vulcan/gpu/engine.py" in rep["by_entry"]["scripts/a2/gates_a2.py"]["unallowed"]
    assert rep["pass"] is False


def test_allow_entries_need_committed_evidence(monkeypatch):
    entries = V.load_allow()
    bad = [dict(e, evidence=["results/does_not_exist.json"]) for e in entries]
    monkeypatch.setattr(V, "load_allow", lambda: bad)
    rep = V.check("HEAD")
    assert "vulcan/gpu/kernels.cu" in rep["unallowed"]
