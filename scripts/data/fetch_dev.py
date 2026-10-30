#!/usr/bin/env python3
"""Fetch the dev sets D2 (PXD000138/PXD009449 dev pools), D3b, D4, D7 and their truth/reference files.
Manifest: data/manifests/dev_sealed.tsv (role=dev rows); hashes: data/manifests/dev_sealed.sha256.  Extra args pass through."""
import os, runpy, sys
sys.argv = [sys.argv[0], "--manifest", "dev_sealed", "--roles", "dev"] + sys.argv[1:]
runpy.run_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fetch_tier_s.py"), run_name="__main__")
