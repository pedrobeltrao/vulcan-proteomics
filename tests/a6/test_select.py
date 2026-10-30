"""A6.3 selection logic on synthetic inputs (no PRIDE listing, no download)."""
import datetime as dt
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts", "a6"))
import select_hn1 as SH  # noqa: E402

SHA = "1a2b3c4d" + "0" * 32
GUARD = {"sha1": set(), "name_size": set(), "names": set(), "sha256": set()}


def test_seed_from_tag_sha():
    assert SH.walk_seed(SHA) == int("1a2b3c4d", 16)


@pytest.mark.parametrize("bad", [None, "", "xyz", "1a2b3c4d"])
def test_walk_refuses_without_tag(bad):
    with pytest.raises(SystemExit):
        SH.walk_seed(bad)


def test_list_project_refuses_before_tag():
    with pytest.raises(SystemExit):
        SH.list_project("PXD000000", None)


def test_walk_order_deterministic_and_seed_dependent():
    accs = [f"PXD{i:06d}" for i in range(200)]
    a = SH.walk_order(accs, 123)
    assert a == SH.walk_order(list(reversed(accs)), 123)
    assert a != SH.walk_order(accs, 124)
    assert sorted(a) == sorted(accs)


def test_file_order_seeded_by_accession_digits():
    files = [{"fileName": f"r{i}.raw"} for i in range(10)]
    o1 = [f["fileName"] for f in SH.file_order(files, 7, "PXD000123")]
    assert o1 == [f["fileName"] for f in SH.file_order(list(reversed(files)), 7, "PXD000123")]
    assert o1 != [f["fileName"] for f in SH.file_order(files, 7, "PXD000124")]


def _f(name, size=int(0.5e9), sha1="", cat="RAW"):
    return {"fileName": name, "fileSizeBytes": size, "fileCategory": {"value": cat}, "checksum": sha1,
            "publicFileLocations": [{"value": "ftp://ftp.pride.ebi.ac.uk/x/" + name}]}


def test_usable_raws_rules_and_guard():
    guard = {"sha1": {"dead"}, "name_size": {("known.raw", int(0.5e9))}, "names": {"seen.raw"}, "sha256": set()}
    files = [_f("ok.raw"), _f("small.raw", int(0.05e9)), _f("big.raw", int(2.5e9)), _f("x_TMT_1.raw"),
             _f("dup.raw", sha1="dead"), _f("known.raw"), _f("seen.raw"), _f("a.mzML", cat="PEAK"), _f("b.RAW")]
    ok, rej = SH.usable_raws(files, guard, set())
    assert [f["fileName"] for f in ok] == ["ok.raw", "b.RAW"]
    ok2, _ = SH.usable_raws([_f("ok.raw")], guard, {("ok.raw", int(0.5e9))})
    assert ok2 == []


def test_cluster_slots_and_union():
    c = SH.Clusters()
    for i in range(5):
        assert c.size_if_added({"lab a"}) == i + 1
        c.admit({"lab a"}, 1, 1, 1)
    assert c.size_if_added({"lab a"}) == 6
    assert c.size_if_added({"lab b"}) == 1
    c.admit({"lab b", "pi x"}, 1, 1, 1)
    c.admit({"pi x"}, 1, 1, 1)
    assert c.size_if_added({"lab b"}) == 3          # linked through the shared key
    assert len(c.tot) == 2


def test_v1_minimum_rule_on_capped_counts():
    J = 25
    ok, d = SH.v1_minimum(np.full(J, 60.0), np.full(J, 6000.0), np.full(J, 15000.0))
    assert ok and d["J"] == J and d["capped_A_cand"] == pytest.approx(150000)
    assert not SH.v1_minimum(np.full(19, 100.0), np.full(19, 10000.0), np.full(19, 1e4))[0]      # J < 20
    assert not SH.v1_minimum(np.full(J, 40.0), np.full(J, 6000.0), np.full(J, 1e4))[0]          # proxy 1000
    assert not SH.v1_minimum(np.full(J, 60.0), np.full(J, 4000.0), np.full(J, 1e4))[0]          # A_cand 100k
    # one dominant cluster is capped at 10% of the capped total -> J_eff and the capped sums fall
    cp = np.r_[np.full(24, 10.0), 5000.0]
    ok, d = SH.v1_minimum(cp, np.full(J, 6000.0), np.full(J, 1e4))
    assert not ok and d["capped_proxy"] < cp.sum()
    # A_cand minimum is max(120000, 10 * G1 n_min)
    assert not SH.v1_minimum(np.full(J, 60.0), np.full(J, 6000.0), np.full(J, 1e4), g1_nmin=20000)[0]


def test_cap_keep_deterministic_and_rate():
    u = np.mean([SH.cap_keep("run", i, 0.3) for i in range(20000)])
    assert abs(u - 0.3) < 0.02
    assert SH.cap_keep("r", 5, 0.5) == SH.cap_keep("r", 5, 0.5)
    assert all(SH.cap_keep("r", i, 1.0) for i in range(100))


T0 = dt.datetime(2026, 10, 28, 15, 30, tzinfo=SH.CET)


def _walk(tmp_path, order, keys, yields, now=None, hard_cap=400):
    listed = []

    def list_fn(acc):
        listed.append(acc)
        return [_f(f"{acc}_{i}.raw") for i in range(4)]

    def admit(acc, files, max_runs):
        y = yields.get(acc)
        if y is None:
            return {"runs": [], "n_downloaded": 1, "note": "no run passed"}
        return {"runs": [dict(run=f"{acc}_{i}", **y) for i in range(min(3, max_runs))], "n_downloaded": 3}
    st = SH.run_walk(order, keys, admit, list_fn, S=1, log_path=str(tmp_path / "log.tsv"),
                     now=now or (lambda: T0), guard=GUARD, hard_cap=hard_cap)
    return st, listed


def test_walk_stops_at_v1_minimum_and_logs(tmp_path):
    order = [f"PXD{i:06d}" for i in range(100)]
    keys = {a: {f"lab{i}"} for i, a in enumerate(order)}
    y = {a: {"A": 5000, "A_cand": 2000, "proxy": 25} for a in order}
    st, listed = _walk(tmp_path, order, keys, y)
    assert st["end"] == "v1 minimum" and st["v1_minimum_reached"]
    assert len(st["admitted"]) >= 20
    log = open(tmp_path / "log.tsv").read().splitlines()
    assert log[0].split("\t") == SH.LOG_COLS and len(log) == 1 + st["walked"]
    assert set(listed) == set(st["listed"])


def test_walk_cluster_full_is_skipped_before_listing(tmp_path):
    order = [f"PXD{i:06d}" for i in range(10)]
    keys = {a: {"same lab"} for a in order}
    y = {a: {"A": 10, "A_cand": 5, "proxy": 1} for a in order}
    st, listed = _walk(tmp_path, order, keys, y)
    assert st["cluster_full"] == 5 and listed == order[:5] and st["end"] == "pool exhausted"


def test_walk_hard_cap(tmp_path):
    order = [f"PXD{i:06d}" for i in range(300)]
    keys = {a: {f"lab{i}"} for i, a in enumerate(order)}
    y = {a: {"A": 10, "A_cand": 5, "proxy": 0} for a in order}
    st, _ = _walk(tmp_path, order, keys, y, hard_cap=400)
    assert st["end"] == "hard cap" and st["runs"] == 400


def test_walk_deadline_reached_not_used(tmp_path):
    order = [f"PXD{i:06d}" for i in range(50)]
    keys = {a: {f"lab{i}"} for i, a in enumerate(order)}
    late = SH.DEADLINE + dt.timedelta(seconds=1)
    n_admit = [0]

    def list_fn(acc):
        return [_f(f"{acc}_{i}.raw") for i in range(4)]

    def admit(acc, files, max_runs):
        n_admit[0] += 1
        return {"runs": [{"run": acc, "A": 10, "A_cand": 5, "proxy": 1}], "n_downloaded": 1}
    # the deadline passes while the 3rd project is being admitted
    st = SH.run_walk(order, keys, admit, list_fn, S=1, log_path=str(tmp_path / "log.tsv"),
                     now=lambda: late if n_admit[0] >= 3 else T0, guard=GUARD)
    assert st["end"] == "walk deadline"
    assert st["listed"] == order[:3]
    assert st["projects"][order[2]] == "reached, not used" and order[2] not in st["admitted"]
    assert st["admitted"] == order[:2]


def test_walk_does_not_start_after_deadline(tmp_path):
    late = SH.DEADLINE + dt.timedelta(seconds=1)
    st, listed = _walk(tmp_path, ["PXD000001"], {}, {}, now=lambda: late)
    assert st["end"] == "walk deadline" and listed == []


def test_hn1_manifest_lists_every_listed_project(tmp_path):
    order = [f"PXD{i:06d}" for i in range(10)]
    keys = {a: ({"same lab"} if i < 7 else {f"l{i}"}) for i, a in enumerate(order)}
    y = {a: {"A": 10, "A_cand": 5, "proxy": 1} for a in order if a != "PXD000001"}
    st, listed = _walk(tmp_path, order, keys, y)
    p = tmp_path / "hn1.tsv"
    SH.write_hn1_manifest(st, str(p))
    rows = [r.split("\t") for r in p.read_text().splitlines()[1:]]
    assert [r[0] for r in rows] == listed and all(r[1] == "dev-null" for r in rows)
    assert dict((r[0], r[2]) for r in rows)["PXD000001"] == "rejected"
    assert "PXD000006" not in [r[0] for r in rows]          # cluster full: metadata only, not burnt (A6.3.7)


def test_exclusions_file_is_committed_shape():
    p = SH.EXCL_TSV
    if not os.path.exists(p):
        pytest.skip("exclusions not built")
    lines = open(p).read().splitlines()
    assert lines[0] == "accession\treasons" and len(lines) > 100
    import hashlib
    rec = open(p + ".sha256").read().split()[0]
    assert hashlib.sha256(open(p, "rb").read()).hexdigest() == rec
