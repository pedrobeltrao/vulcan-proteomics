"""Cross-linking MS (DESIGN.md §4 L2-L4 XL-MS): linker definitions, MS2 doublet detection for MS-cleavable linkers,
GPU candidate generation and scoring on the L1 index (cleavable: two stub-modified linear searches + pairing;
non-cleavable: alpha-first open search + beta lookup), linear / mono-link / loop-link competition, and label-free
CSM tables. Labels (target/decoy, library truth) are joined only in L4 (vulcan.stats.xl, scripts/m5)."""
