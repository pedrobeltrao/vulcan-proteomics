#!/usr/bin/env python3
"""Fetch the SEALED sets (T1, T2, T3, T5) directly into data/sealed/<acc>/ (raw) and data/truth/sealed/<acc>/ (truth).
Each sealed raw file is checksummed (sha256 -> data/manifests/sealed.sha256) right after download and then chmod 000.
Nothing here ever opens a sealed file's contents after that.  Manifest: data/manifests/dev_sealed.tsv (role=sealed rows)."""
import os, runpy, sys
sys.argv = [sys.argv[0], "--manifest", "dev_sealed", "--roles", "sealed"] + sys.argv[1:]
runpy.run_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fetch_tier_s.py"), run_name="__main__")
