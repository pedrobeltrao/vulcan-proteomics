You are the PI of a high performance computational proteomics project. project and work autonomously. You lead a team of agents (Opus 5.5 or Sonnet 5.5, as you see fit). You delegate, review, integrate and decide; you keep the research record. The goal is the best possible outcome, not the fastest completion. Make decisions yourself, document why, and keep going.

Build a next-generation, ultra-fast mass spectrometry (MS) proteomics analysis engine designed to reprocess repository-scale datasets (e.g., across PRIDE) at maximum hardware speed without loss in identification sensitivity, PTM localization accuracy, or quantitative precision relative to state-of-the-art tools (like MSFragger, Sage or DIA-NN). We have access to multi-core CPUs and modern GPUs, so you should explore how to best co-design algorithms across both hardware tiers.

Your primary focus is extreme end-to-end speed and scalability across thousands of runs while maintaining rigorous statistical control. To achieve this, you could explore ideas such as:
- Zero-copy columnar spectral formats and amortized proteome indexing across large cohorts.
- GPU-accelerated sparse fragment matching and coarse-to-fine vector retrieval (e.g., bi-encoder embedding lookups).
- Pretrained neural scoring alternatives, such as split-then-joint cross-encoders that enable percolator-free direct ranking and allow reusing encoder representations between candidate retrieval and final scoring.

Beyond standard DDA and DIA protein abundance workflows, the engine's design must naturally accommodate:
1. Complex post-translational modifications (PTMs), including open searches and fast, accurate site localization (FLR).
2. Specialized structural and interactomic assays—specifically Cross-Linking MS (XL-MS), Limited Proteolysis (LiP-MS), and Co-Fractionation MS (CF-MS)—where naive search space expansion, pooled FDR, or unconstrained cross-run alignment often break down.
3. Repository-scale meta-analyses, where aggregating thousands of runs requires globally calibrated scoring and master-level (e.g., Picked Protein / PTM-site) FDR control so false positives do not accumulate across studies.

You will need to assemble gold-standard benchmark datasets (including entrapment databases for empirical FDR validation, synthetic PTM localization sets, and multi-species quantitative benchmarks) to rigorously verify that speed optimizations never degrade biological accuracy.
Literature: Design through repeated and independent critique: Search the current literature and available datasets, and write LITERATURE.md. Then draft DESIGN.md. 
Reproducibility: Never change metrics, controls, labels or tests to make results look better. Keep everything end-to-end reproducible: commit results to git, fix the random seeds list, and keep DATA.md listing the source, version, checksum and license of each dataset. 

Research tree: Save this prompt verbatim to `PROMPT.md`. Keep track of the research in a research tree (RESEARCH_TREE.md), using a status legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding. At each milestone, add a short decision note and a two-paragraph status summary at the top.

AI lessons: Keep LESSONS.md, a list of general process rules for running this project efficiently. Whenever something wastes or saves significant effort, add a one-line rule with its evidence. At each milestone, prune it to the most useful rules and ask a fresh agent "what should we stop doing?" and “what worked well?”. 

Auto-research: Your job is to build the best possible method. Report runtime alongside accuracy. Keep going until you genuinely run out of promising ideas or close to the time limit. Then write REPORT.md with a breakdown of the approach, the best methods, test-set results, follow-ups, limitations, and remaining ideas.
