#!/usr/bin/env python3
"""Fetch dev set D12 (M6d MBR false-transfer truth): PXD014415 (Lim, Paulo, Gygi 2019 JPR) Orbitrap Fusion Lumos DDA, 90 min:
6 donor runs of sample HY (human + 10% yeast) + 3 acceptor runs of sample H (human only, yeast absent) + community SDRF.
Manifest: data/manifests/d12.tsv; hashes: data/manifests/d12.sha256.  PRIDE SHA-1 is verified on arrival by the generic fetcher.
Extra args pass through (e.g. --jobs 4, --check)."""
import os, runpy, sys
sys.argv = [sys.argv[0], "--manifest", "d12", "--roles", "dev"] + sys.argv[1:]
runpy.run_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fetch_tier_s.py"), run_name="__main__")
