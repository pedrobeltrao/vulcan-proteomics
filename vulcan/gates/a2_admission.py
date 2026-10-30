"""Amendment A2.2 / A2.3 admission logic (pure functions; I/O and GPU searches live in scripts/a2/*).

Label discipline (P1, P2, A2.9.4):
  * The gate database T∪D is seen only through vulcan.gates.a2_identity.GateDB (label-free) and, for searches, the
    label-free index side (vulcan.index.build.L1Index: peptidoforms, residues, partner).  Nothing in this module or in
    the admission scripts opens the gate index's labels.arrow / proteins.arrow (tests/a2/test_no_label_counts.py).
  * Own-species searches (A2.2a, A2.3 step 3) DO need own-index labels (own target vs own decoy, protein -> source).
    They are read only through :func:`own_labels`, which refuses any index not built by scripts/a2/own_search.py
    (meta.json "a2_role": "own").  Own-search target/decoy counts are not gate per-label counts; still, no admission
    output stores them (only certified / detected counts).
  * Every admission output is a per-spectrum membership code or a label-free count (ADMISSION_COUNT_KEYS).

Membership codes per null spectrum (uint8):
  0 NOT_CERTIFIED  fails (a)                               -> B candidate
  1 A              passes (a), (b), (d), (c-pre)          -> stratum A (pending the final (c) in finalize)
  2 REM_B_F        certified, removed by (b) max f(c) >= 0.20
  3 REM_B_HAM      certified, removed by (b) Hamming/isobaric
  4 REM_D          certified, removed by (d) family screen
  5 REM_C_PRE      certified, removed by (c) on the exhaustive hyperscore top candidate and its pair partner
  codes 2-5 are B candidates.  Order of evaluation (first failing step gives the code): a, b-f, b-ham, d, c-pre.
Final (c) (finalize, after the mixtures are searched, before any gate statistic): a spectrum stays in A / B only if
the matched peptide AND its pair partner of EVERY compared pipeline are not A1.3-near (source proteome ∪ contaminants
∪ Panel H ∪ Panel B); stratum B = B candidates passing the final (c).
"""
from __future__ import annotations

import hashlib
import json
import os

import numpy as np

NOT_CERTIFIED, A, REM_B_F, REM_B_HAM, REM_D, REM_C_PRE = 0, 1, 2, 3, 4, 5
CODE_NAME = {0: "not_certified", 1: "A", 2: "removed_b_f", 3: "removed_b_hamming", 4: "removed_d_family",
             5: "removed_c_pre"}
# The only count keys an admission output may carry (label-free).  Tests assert outputs use only these.
ADMISSION_COUNT_KEYS = ("n_spectra", "n_certified", "removed_b_f", "removed_b_hamming", "removed_d_family",
                        "removed_c_pre", "n_A", "n_A_with_candidate", "tail_proxy", "n_B_candidates",
                        "n_panel_members_detected", "n_trigger_members_detected")

# certification thresholds (A2.2a)
Q_PSM_CERT = 0.001
Q_PEP_CERT = 0.01
Q_PEP_MARKER = 0.01
MARKER_MIN_PEPTIDES = 2
MARKER_EXCLUDE_AT = 2


# ------------------------------------------------------------------------------------------------ label guard
def own_labels(index_dir: str) -> dict:
    """Own-species index labels (A2.2a only).  Refuses any index whose meta.json lacks "a2_role": "own"."""
    meta = json.load(open(os.path.join(index_dir, "meta.json")))
    if meta.get("a2_role") != "own":
        raise PermissionError(f"own_labels: {index_dir} is not an A2 own-species index (a2_role != 'own')")
    import pyarrow as pa
    with pa.OSFile(os.path.join(index_dir, "labels.arrow"), "rb") as src:
        lab = pa.ipc.open_file(src).read_all()
    with pa.OSFile(os.path.join(index_dir, "proteins.arrow"), "rb") as src:
        prot = pa.ipc.open_file(src).read_all()
    acc = open(os.path.join(index_dir, "accessions.txt")).read().split("\n")
    pidx = np.load(os.path.join(index_dir, "protein_idx.npy"))
    return {"is_decoy": lab.column("is_decoy").to_numpy(zero_copy_only=False),
            "tdc_pair": lab.column("tdc_pair").to_numpy(zero_copy_only=False),
            "target_seq": lab.column("target_seq").to_numpy(zero_copy_only=False),
            "prot_peptide": prot.column("peptide").to_numpy(zero_copy_only=False),
            "prot_off": prot.column("prot_off").to_numpy(), "prot_len": prot.column("prot_len").to_numpy(),
            "protein_idx": pidx, "accessions": acc}


def peptide_sources(lab: dict) -> dict:
    """I/L-collapsed target peptide -> set of first-token accessions (OWN|..., CRAP|..., CON|..., PNL|...)."""
    acc = lab["accessions"]
    out = {}
    for pep, o, n in zip(lab["prot_peptide"], lab["prot_off"], lab["prot_len"]):
        out[pep] = {acc[int(i)] for i in lab["protein_idx"][int(o):int(o) + int(n)]}
    return out


def source_class(accs) -> dict:
    """Which sources an accession set touches."""
    return {"own": any(a.startswith("OWN|") for a in accs),
            "contam": any(a.startswith(("CRAP|", "CON|")) for a in accs),
            "panel": any(a.startswith("PNL|") for a in accs)}


# ------------------------------------------------------------------------------------------------ (a) certification
def qvalues_psm(score, is_decoy):
    """PSM-level TDC q over one run's PSM winners (one per spectrum); (D+1)/T, monotone (vulcan.stats.tdc)."""
    from vulcan.stats.tdc import tdc_qvalues
    return tdc_qvalues(np.asarray(score, float), np.asarray(is_decoy, bool))


def certify_run(spec, unit, pair, score, is_decoy, il_peptide, sources, near_contam, tiebreak,
                q_psm: float = Q_PSM_CERT, q_pep: float = Q_PEP_CERT) -> dict:
    """A2.2(a) for one own-search run.

    Inputs per PSM row (all candidates of all hypotheses, production report depth): spectrum id, peptidoform unit id,
    pair key (target unit and its decoy share it), hyperscore, own-index decoy flag, I/L-collapsed stripped peptide,
    tiebreak key (hash of the peptidoform, seed fixed in config).  sources: peptide -> accession set;
    near_contam: callable(list[str]) -> bool array, A1.3-near (6-mer or Hamming <= 2) to cRAP ∪ Frankenfield ∪ Panel H
    ∪ Panel B (vulcan.gates.homolog.HomologFilter.near).
    Returns rows (into the input) of certified rank-1 PSMs and diagnostics that carry no own target/decoy split."""
    from vulcan.stats.tdc import double_competition, psm_competition
    spec = np.asarray(spec)
    score = np.asarray(score, float)
    dec = np.asarray(is_decoy, bool)
    w = psm_competition(spec, score, tiebreak=tiebreak)          # rank-1 PSM per spectrum
    qp = qvalues_psm(score[w], dec[w])
    dc = double_competition(spec, unit, pair, score, dec, tiebreak=tiebreak)
    unit = np.asarray(unit)
    acc_units = np.unique(unit[dc.row[(dc.q <= q_pep) & ~dc.is_decoy]])
    ok = (~dec[w]) & (qp <= q_psm) & np.isin(unit[w], acc_units)
    peps = [il_peptide[int(r)] for r in w]
    own_only = np.array([source_class(sources.get(p, ()))["own"] for p in peps], bool)
    ok &= own_only
    cand = np.flatnonzero(ok)
    if cand.size:
        nr = np.asarray(near_contam([peps[i] for i in cand]), bool)
        ok[cand[nr]] = False
    return {"rows": w[ok], "n_rank1": int(w.size), "n_certified": int(ok.sum())}


# ------------------------------------------------------------------------------------------------ A2.3 step 3
def marker_screen(accepted_target_peptides_by_run, sources, near_own, trigger_panels=("H", "B"),
                  min_peptides: int = MARKER_MIN_PEPTIDES, exclude_at: int = MARKER_EXCLUDE_AT) -> dict:
    """A2.3 step 3, pooled over the project's runs.

    accepted_target_peptides_by_run: {run: iterable of I/L-collapsed stripped target peptides accepted at
    peptide-level q <= 0.01 in that run's own search}.  Pooling = union over runs (stricter than per-run).
    A peptide counts for EVERY panel member it maps to (shared peptides count for each: stricter) unless it is A1.3-near
    the own proteome (near_own: callable(list[str]) -> bool array).  A member is detected with >= 2 distinct peptides.
    Exclude iff >= 2 detected members belong to the trigger panels (H non-keratin, B)."""
    from vulcan.gates.a2_panels import parse_panel_accession
    peps = sorted({p for v in accepted_target_peptides_by_run.values() for p in v})
    panel_peps = [p for p in peps if any(a.startswith("PNL|") for a in sources.get(p, ()))]
    keep = ~np.asarray(near_own(panel_peps), bool) if panel_peps else np.zeros(0, bool)
    per = {}
    for p, k in zip(panel_peps, keep):
        if not k:
            continue
        for a in sources.get(p, ()):
            m = parse_panel_accession(a)
            if m is not None:
                per.setdefault(m, set()).add(p)
    detected = sorted(m for m, s in per.items() if len(s) >= min_peptides)
    trig = [m for m in detected if m.split("|")[0] in trigger_panels]
    return {"detected_members": detected, "peptides_per_detected_member": {m: len(per[m]) for m in detected},
            "trigger_members": trig, "n_panel_members_detected": len(detected),
            "n_trigger_members_detected": len(trig), "exclude": bool(len(trig) >= exclude_at)}


# ------------------------------------------------------------------------------------------------ combination
def combine(n_spectra: int, certified_idx, rem_b_f, rem_b_ham, rem_d, rem_c_pre) -> np.ndarray:
    """Membership codes from the per-certified-spectrum screen results (arrays aligned with certified_idx)."""
    codes = np.zeros(n_spectra, np.uint8)
    ci = np.asarray(certified_idx, np.int64)
    c = np.full(ci.size, A, np.uint8)
    for flag, code in ((rem_c_pre, REM_C_PRE), (rem_d, REM_D), (rem_b_ham, REM_B_HAM), (rem_b_f, REM_B_F)):
        c[np.asarray(flag, bool)] = code          # applied in reverse so the FIRST failing step wins
    codes[ci] = c
    return codes


def label_free_counts(codes: np.ndarray, has_candidate=None, tail_proxy=None) -> dict:
    out = {"n_spectra": int(codes.size), "n_certified": int((codes >= 1).sum()),
           "removed_b_f": int((codes == REM_B_F).sum()), "removed_b_hamming": int((codes == REM_B_HAM).sum()),
           "removed_d_family": int((codes == REM_D).sum()), "removed_c_pre": int((codes == REM_C_PRE).sum()),
           "n_A": int((codes == A).sum()), "n_B_candidates": int((codes != A).sum())}
    if has_candidate is not None:
        out["n_A_with_candidate"] = int(((codes == A) & np.asarray(has_candidate, bool)).sum())
    if tail_proxy is not None:
        out["tail_proxy"] = int(((codes == A) & np.asarray(tail_proxy, bool)).sum())
    assert set(out) <= set(ADMISSION_COUNT_KEYS)
    return out


def membership_hash(codes_by_run: dict) -> str:
    """sha256 over runs sorted by name: name, n, codes bytes."""
    h = hashlib.sha256()
    for run in sorted(codes_by_run):
        c = np.ascontiguousarray(codes_by_run[run], np.uint8)
        h.update(run.encode() + b"\0" + str(c.size).encode() + b"\0" + c.tobytes())
    return h.hexdigest()


# ------------------------------------------------------------------------------------------------ stopping rule
def stopping_state(projects: list[dict], cfg: dict) -> dict:
    """A2.9.3 stopping rule from label-free per-project counts (admitted projects only, walk order).
    Stop when ALL of: >= min_projects admitted, sum n_A >= min_A, sum n_A_with_candidate >= min_A_with_candidate,
    sum tail_proxy >= min_tail_proxy; or when cap admitted projects are reached."""
    n = len(projects)
    tot = {k: int(sum(p[k] for p in projects)) for k in ("n_A", "n_A_with_candidate", "tail_proxy")}
    crit = {"projects": n >= cfg["min_projects"], "n_A": tot["n_A"] >= cfg["min_A"],
            "n_A_with_candidate": tot["n_A_with_candidate"] >= cfg["min_A_with_candidate"],
            "tail_proxy": tot["tail_proxy"] >= cfg["min_tail_proxy"]}
    met = all(crit.values())
    cap = n >= cfg["cap_projects"]
    return {"n_admitted": n, "totals": tot, "criteria": crit, "stop": bool(met or cap),
            "reason": "criteria met" if met else ("cap reached" if cap else "continue")}
