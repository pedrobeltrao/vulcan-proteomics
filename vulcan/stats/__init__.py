"""vulcan.stats: FDR estimation and verification (DESIGN.md §4 L4, §6).

Modules
-------
tdc        target-decoy competition, paired double competition, multi-decoy TDC
groupwalk  Group-walk and pre-registered stratum merge rules
reset      RESET (single- and two-decoy) with pluggable semi-supervised learner
entrapment lower-bound / combined / paired entrapment FDP and the §6.3 decision rule
protein    Picked Protein Group FDR (+ classic protein FDR for comparison)
flr        decoy-amino-acid FLR
xl         cross-link FDR (TD-DD)/TT and TD/TT per level, self vs heteromeric
g2         differential G2 null test, diagnostics, power / sample size
familymax  family-max single-pass entrapment reconstruction
global_fdr L6 repository-scale FDR: best evidence across runs -> global paired TDC; naive-merge comparator
sim        Monte-Carlo simulators used by tests/stats and scripts/run_stats_sims.py

See vulcan/stats/README.md for assumptions and citations.
"""
from ._util import CALIB_DECOY_SEED, DECOY_SEEDS, ENTRAP_SEEDS_DEV, ENTRAP_SEEDS_SEALED, SEEDS, make_rng
from .tdc import compete, count_discoveries, double_competition, psm_competition, tdc_qvalues

__all__ = [
    "SEEDS", "DECOY_SEEDS", "ENTRAP_SEEDS_DEV", "ENTRAP_SEEDS_SEALED", "CALIB_DECOY_SEED", "make_rng",
    "compete", "psm_competition", "tdc_qvalues", "double_competition", "count_discoveries",
]
