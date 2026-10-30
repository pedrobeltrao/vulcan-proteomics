#!/usr/bin/env python3
"""Fetch dev set D9 (M6 LFQ dev benchmark): ProteoBench DDA-Astral HYE A/B x 3 + 3 human-only DDA-Astral runs
(MBR false-transfer control, ProteoBench DDA-astral-entrapment) + the ProteoBench HYE FASTA.
Manifest: data/manifests/d9.tsv; hashes: data/manifests/d9.sha256.  Extra args pass through (e.g. --check)."""
import os, runpy, sys
sys.argv = [sys.argv[0], "--manifest", "d9", "--roles", "dev"] + sys.argv[1:]
runpy.run_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fetch_tier_s.py"), run_name="__main__")
