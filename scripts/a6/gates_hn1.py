#!/usr/bin/env python3
"""A6 step 5 (agent G only; ONE session): v1 gate rows on HN1 (A6.6.1) and the status file (A6.2).

This is the ONLY A6 module that reads labels of null-copy rows (A6.4.6.3, review B2): per-mixture double-competition
1% thresholds of every pipeline (as scripts/a2/gates_a2.collect), tail membership, target/decoy of winners, the G1
sequence-classifier logit, and the per-label G3 bests.  Statistics: scripts/a6/a6_stats.py.

Refuses to run unless (A6.8.4):
  * the prereg-v0.4.5 tag exists and HEAD contains it; the working tree is CLEAN;
  * config/a6.yaml is committed and unmodified; the v1 scorer pickle and the A2 sequence-classifier hash match;
  * the step-4 record (results/2026-10-28_A6_hn1/step4_record.json) is committed and unmodified and every artefact it
    lists hashes to the recorded sha256;
  * the clock is before the 17:55 CET hard stop.
Outputs: results/2026-10-28_A6_hn1/gates_hn1.json, then a commit, then a6_status.json with that commit sha, and a
second commit.  Status mapping (A6.7): verdict PASS / FAILED / not evaluable / invalid pool; any exception or the hard
stop AFTER the first label read -> "not completed" (completed test); before it -> "label-free stop".

    conda_envs/pe/bin/python scripts/a6/gates_hn1.py --record results/2026-10-28_A6_hn1/step4_record.json

Step-4 record layout (written by agent S; label-free):
  {"tag_commit": sha, "mixtures": {run: {"psms": path, "feat": path, "npz": path, "project": acc, "cluster": id,
    "instrument": class, "capped_A": path-of-npy (bool over the run's A copies, A6.3.5)}}, "sha256": {path: sha}}
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import subprocess
import sys
import traceback

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
for _p in (ROOT, HERE, os.path.join(ROOT, "scripts", "m3"), os.path.join(ROOT, "scripts", "m3e"),
           os.path.join(ROOT, "scripts", "m2")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import status as STATUS  # noqa: E402

TAG = "prereg-v0.4.5"
CET = dt.timezone(dt.timedelta(hours=1))
HARD_STOP = dt.datetime(2026, 10, 28, 17, 55, tzinfo=CET)
OUT_DIR = os.path.join(ROOT, "results", "2026-10-28_A6_hn1")
GATE_OUT = os.path.join(OUT_DIR, "gates_hn1.json")
CONFIG = os.path.join(ROOT, "config", "a6.yaml")
P_A, P_FS = 20, 22
PIPES = ("fixed", "v1", "v1_noRT")
TB_SEED = 11


# ------------------------------------------------------------------------------------------------ guards
def git(*a, check=True):
    return subprocess.run(["git", "-C", ROOT, *a], capture_output=True, text=True, check=check).stdout.strip()


def committed_unmodified(path: str) -> bool:
    rel = os.path.relpath(path, ROOT)
    ls = subprocess.run(["git", "-C", ROOT, "ls-files", "--error-unmatch", rel], capture_output=True)
    st = subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "--", rel], capture_output=True, text=True)
    return ls.returncode == 0 and not st.stdout.strip()


def guards(record_path: str, now=None) -> dict:
    """Every A6.8.4 precondition; raises SystemExit (a label-free refusal) on the first failure."""
    import yaml
    now = now or dt.datetime.now(CET)
    if now >= HARD_STOP:
        raise SystemExit("gates_hn1: past the 17:55 CET hard stop")
    tag = subprocess.run(["git", "-C", ROOT, "rev-parse", "--verify", "--quiet", TAG + "^{commit}"],
                         capture_output=True, text=True).stdout.strip()
    if not tag:
        raise SystemExit(f"gates_hn1: tag {TAG} absent (A6 not tagged)")
    if subprocess.run(["git", "-C", ROOT, "merge-base", "--is-ancestor", tag, "HEAD"]).returncode != 0:
        raise SystemExit(f"gates_hn1: HEAD does not contain {TAG}")
    if git("status", "--porcelain", "--untracked-files=no"):
        raise SystemExit("gates_hn1: working tree is not clean")
    if git("diff", "--name-only", tag, "HEAD", "--", "scripts/a6", "config/a6.yaml", "tests/a6"):
        raise SystemExit("gates_hn1: an A6 file changed after the tag")
    if not committed_unmodified(CONFIG):
        raise SystemExit("gates_hn1: config/a6.yaml not committed/unmodified")
    cfg = yaml.safe_load(open(CONFIG))
    sc = cfg["scorers"]["v1_static_gbm_pred"]
    if STATUS.sha256_file(os.path.join(ROOT, sc["model"])) != sc["sha256"]:
        raise SystemExit("gates_hn1: v1 scorer hash mismatch")
    from vulcan.gates import a2_config
    a2_config.check_seqclf(a2_config.load())
    if not committed_unmodified(record_path):
        raise SystemExit("gates_hn1: step-4 record not committed/unmodified (P2)")
    rec = json.load(open(record_path))
    for p, h in rec["sha256"].items():
        fp = p if os.path.isabs(p) else os.path.join(ROOT, p)
        if not os.path.exists(fp) or STATUS.sha256_file(fp) != h:
            raise SystemExit(f"gates_hn1: step-4 artefact {p} missing or modified")
    return {"tag": tag, "cfg": cfg, "record": rec}


# ------------------------------------------------------------------------------------------------ step 5 (labels)
class LabelRead:
    """Flag: set at the FIRST label read; decides 'not completed' (completed test) vs 'label-free stop' (A6.7)."""
    started = False


def collect(rec: dict, idx: str):  # pragma: no cover (data)
    """Per null copy (A and FS) of every mixture: pool, project, cluster, instrument, capped flag and, per pipeline,
    win_t / score / tail1 / logit; per A copy the G3 bests.  THE ONLY LABEL READS OF NULL ROWS IN A6."""
    import pyarrow.parquet as pq
    import m3_data as MD
    import seqclf
    from e2_correctness import _best_per_label
    from vulcan.gates import stats as GS
    from vulcan.index.build import L1Index
    from vulcan.stats._util import hash_tiebreak
    from vulcan.stats.tdc import double_competition
    l1 = L1Index(idx)
    clf = seqclf.load()
    LabelRead.started = True
    lab = MD.labels(idx)
    dec_all, pair_all = lab["is_decoy"], lab["tdc_pair"]
    D = {k: [] for k in ("pool", "project", "cluster", "instrument", "capped", "run")}
    P = {p: {k: [] for k in ("win_t", "score", "tail1", "logit")} for p in PIPES}
    G3 = {k: [] for k in ("bt", "bd", "kt", "kd")}
    thr_rec = {}
    for run, m in sorted(rec["mixtures"].items()):
        t = pq.read_table(os.path.join(ROOT, m["psms"]))
        psm = {c: t.column(c).to_numpy() for c in ("spec_idx", "pf_row", "pf_id")}
        ft = pq.read_table(os.path.join(ROOT, m["feat"]), columns=["rt_pred"])
        rt_pred = ft.column("rt_pred").to_numpy()
        z = np.load(os.path.join(ROOT, m["npz"]))
        prov = z["prov"]
        src = z["src_store_idx"]
        capA = np.load(os.path.join(ROOT, m["capped_A"]))           # bool over A copies in order
        nA = int((prov == P_A).sum())
        a_pos = np.full(prov.size, -1, np.int64)
        a_pos[prov == P_A] = np.arange(nA)
        a_pos[prov == P_FS] = np.arange(int((prov == P_FS).sum()))     # FS copy k <-> A copy k (same order)
        spec = psm["spec_idx"].astype(np.int64)
        rows = psm["pf_row"].astype(np.int64)
        dec = dec_all[rows]
        tb = hash_tiebreak(psm["pf_id"], TB_SEED)
        w_fixed = z["win_fixed"]
        sw = spec[w_fixed]
        null = np.isin(prov[sw], [P_A, P_FS])
        sel_spec = sw[null]
        D["pool"].append(np.where(prov[sel_spec] == P_A, "A", "FS"))
        D["project"].append(np.full(sel_spec.size, m["project"], object))
        D["cluster"].append(np.full(sel_spec.size, str(m["cluster"]), object))
        D["instrument"].append(np.full(sel_spec.size, m["instrument"], object))
        D["run"].append(np.full(sel_spec.size, run, object))
        D["capped"].append(capA[a_pos[sel_spec]])
        thr_rec[run] = {}
        for p in PIPES:
            sc = z[f"score_{p}"]
            w = z[f"win_{p}"]
            assert np.array_equal(spec[w], sw), "winner spectra differ between pipelines"
            r = double_competition(spec, psm["pf_id"], pair_all[rows], sc, dec, seed=TB_SEED, tiebreak=tb)
            thr = GS.tail_threshold(r.score, r.is_decoy, r.q, 0.01)
            thr_rec[run][p] = thr
            wn = w[null]
            P[p]["win_t"].append(~dec[wn])
            P[p]["score"].append(sc[wn])
            P[p]["tail1"].append(sc[wn] >= thr)
            P[p]["logit"].append(clf.logit_rows(l1, rows[wn], rt_pred[wn]))
        # G3 (A copies only; per A copy, aligned with the A rows of sel_spec)
        bt_x, bd_x = _best_per_label(z["g3_exh_spec"], z["g3_exh_hs"], dec_all[z["g3_exh_row"]], nA)
        bt_k, bd_k = _best_per_label(z["g3_prod_spec"], z["g3_prod_hs"], dec_all[z["g3_prod_row"]], nA)
        ap = a_pos[sel_spec]
        isA = prov[sel_spec] == P_A
        for k, v in (("bt", bt_x), ("bd", bd_x), ("kt", bt_k), ("kd", bd_k)):
            G3[k].append(np.where(isA, v[np.clip(ap, 0, max(nA - 1, 0))] if nA else -np.inf, -np.inf))
    cat = {k: np.concatenate(v) if v else np.zeros(0) for k, v in D.items()}
    out = {"pool": cat["pool"], "project": cat["project"], "cluster": cat["cluster"],
           "instrument": cat["instrument"], "run": cat["run"],
           "fixed": {k: np.concatenate(v) for k, v in P["fixed"].items()},
           "v1": {k: np.concatenate(v) for k, v in P["v1"].items()},
           "v1_noRT": {k: np.concatenate(v) for k, v in P["v1_noRT"].items()},
           "g3": {k: np.concatenate(v) for k, v in G3.items()}}
    return out, cat["capped"].astype(bool), thr_rec


def evaluate(D: dict, capped: np.ndarray) -> dict:
    """Both pools (A6.3.5 intersection-union) for v1; noRT and the rest are reported only."""
    import a6_stats as S
    allm = np.ones(len(D["pool"]), bool)
    res_c = S.evaluate_pool(D, capped)
    res_u = S.evaluate_pool(D, allm)
    verdict = S.v1_verdict(res_c, res_u)
    Dn = dict(D, v1=D["v1_noRT"])
    reported = {"noRT_capped": S.pool_row_verdicts(S.evaluate_pool(Dn, capped)),
                "noRT_uncapped": S.pool_row_verdicts(S.evaluate_pool(Dn, allm))}
    return {"v1": {"capped": res_c, "uncapped": res_u, "verdict": verdict}, "reported": reported}


def commit(paths, msg) -> str:  # pragma: no cover (git)
    git("add", "--", *[os.path.relpath(p, ROOT) for p in paths])
    git("commit", "-m", msg)
    return git("rev-parse", "HEAD")


def main():  # pragma: no cover
    ap = argparse.ArgumentParser()
    ap.add_argument("--record", required=True)
    a = ap.parse_args()
    g = guards(a.record)
    idx = g["record"].get("index", "/tmp/vulcan_e2/index/human_contam_prod")
    try:
        D, capped, thr = collect(g["record"], idx)
        if dt.datetime.now(CET) >= HARD_STOP:
            raise TimeoutError("hard stop reached during step 5")
        res = evaluate(D, capped)
        res["thresholds"] = thr
        res["tag_commit"] = g["tag"]
        os.makedirs(OUT_DIR, exist_ok=True)
        json.dump(res, open(GATE_OUT, "w"), indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
        sha = commit([GATE_OUT], "A6 step 5: v1 gate rows on HN1 (gates_hn1.json)")
        STATUS.write_status(res["v1"]["verdict"]["status"], step5_commit=sha, gate_output=GATE_OUT)
    except BaseException as e:  # noqa: BLE001
        st = "stopped: node time limit" if isinstance(e, TimeoutError) else \
            ("not completed" if LabelRead.started else "label-free stop")
        STATUS.write_status(st, note=f"{type(e).__name__}: {e}; {traceback.format_exc()[-1500:]}")
        commit([STATUS.STATUS_PATH], f"A6 step 5: status {st}")
        raise
    commit([STATUS.STATUS_PATH], f"A6 status: {res['v1']['verdict']['status']}")


if __name__ == "__main__":
    main()
