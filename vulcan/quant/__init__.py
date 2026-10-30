"""L5 quantification (DESIGN §4 L5, M6): LFQ from MS1 features (link, RT alignment, opt-in MBR with transfer
q-values and provenance, Top-N / MaxLFQ / directLFQ roll-up), TMT reporter extraction (SPS-MS3 linkage,
isolation purity, impurity correction), and the §6.8 quant metrics.  Nothing here reads a target/decoy label
except where stated (MBR decoy-peptide donors are passed in explicitly by the caller, from L4)."""
