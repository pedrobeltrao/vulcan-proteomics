#!/usr/bin/env python3
"""Fetch dev set D11 (M6d second LFQ benchmark): PXD003881 IonStar (Shen 2018 PNAS) human + E. coli spike-in, Orbitrap Fusion DDA,
group A (1x, 3% E. coli) x 4 and group E (3x, 9% E. coli) x 4 + the submitter human/E. coli FASTA (reference only).
Manifest: data/manifests/d11.tsv; hashes: data/manifests/d11.sha256.  PRIDE SHA-1 is verified on arrival by the generic fetcher.
Extra args pass through (e.g. --jobs 4, --check)."""
import os, runpy, sys
sys.argv = [sys.argv[0], "--manifest", "d11", "--roles", "dev"] + sys.argv[1:]
runpy.run_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fetch_tier_s.py"), run_name="__main__")
