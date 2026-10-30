#!/usr/bin/env python3
"""Append scripts/a3/v1_closure_allow.yaml entries for every changed, not-yet-allowed file of the v1 closure at the
final A3/A5 freeze code commit (D32).  Each entry pins the exact diff (prereg-v0.4.2 -> <commit>) and names the final
freeze evidence (results/2026-10-28_A3A5_final_freeze/dryrun_D1.json: whole driver at prereg-v0.4.2 vs at the freeze
code, v1 artefacts byte-identical, v1 metrics identical, A1 comparator artefacts identical from the same searches) plus
the rehearsal's D0 evidence at e0b94aa, under the D32 deviation (D0 not re-run at the final commit).  The evidence
files are committed AFTER the freeze code commit, so the guard passes only at the evidence commit (by design: it
checks that evidence is committed at the commit it examines).

    python scripts/a5/allow_entries.py --commit <freeze code commit>
"""
import argparse
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts", "a3"))
import v1_closure_guard as G  # noqa: E402

EVIDENCE = ["results/2026-10-28_A3A5_final_freeze/dryrun_D1.json",
            "results/2026-10-28_A3_freeze/D0/v1_artefacts_tag_vs_freeze.json",
            "results/2026-10-28_A3_freeze/D0/a1_comparator_artefacts_tag_vs_freeze.json"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--commit", required=True)
    a = ap.parse_args()
    full = G.git("rev-parse", a.commit).strip()
    # evaluate as if the evidence were committed: list changed files whose (path, diff) has no entry yet
    allow = {(x["path"], x["diff_sha256"]) for x in G.load_allow()}
    cl = sorted(G.closure(G.TAG) | G.closure(full))
    new = []
    for f in cl:
        sha, _ = G.diff_sha(full, f)
        if sha and (f, sha) not in allow:
            new.append((f, sha))
    if not new:
        print("nothing to add")
        return
    txt = ("\n# A3/A5 FINAL freeze (2026-10-28, A5 implementer; D32): changed files of the v1 closure at the freeze code "
           f"commit {full[:7]}.\n# Evidence: the whole driver at prereg-v0.4.2 vs at {full[:7]} (--a3) on D1 (final "
           "dry run) and the rehearsal's D0 byte-identity at e0b94aa; D0 is not re-run at the final commit (D32 "
           "deviation, argued per diff in results/2026-10-28_A3A5_final_freeze/NOTES.md).\n")
    for f, sha in new:
        txt += (f"- path: {f}\n  diff_sha256: {sha}\n  kind: bit_identity\n  commit: {full[:7]}\n"
                f"  evidence: [{', '.join(EVIDENCE)}]\n"
                f"  summary: \"final A3/A5 freeze; D1 end-to-end bit identity at {full[:7]}; D0 per D32 (NOTES)\"\n")
    open(G.ALLOW, "a").write(txt)
    print(f"added {len(new)} entries: " + ", ".join(f for f, _ in new))


if __name__ == "__main__":
    main()
