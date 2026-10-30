# Stats core: Monte-Carlo verification (DESIGN §6.7)

10⁴ replicates per estimator. Mean FDP at α = 1% (all ≤ α):

| Estimator | Mean FDP |
|---|---|
| TDC | 0.0064 |
| TDC, near-exact regime | 0.0097 |
| Paired double competition | 0.0082 |
| Group-walk | 0.0078 (+4.5% IDs vs pooled TDC) |
| RESET single-decoy | 0.0078 |
| RESET two-decoy | 0.0080 |
| Picked protein | 0.0081 |
| Picked protein groups, isoforms + rescue | 0.0060 |
| Multi-decoy (d=3) | 0.0089 |

**Biased null (target-win p = 0.55):**
- TDC FDP = 0.0119 at α = 1% and 0.0609 at α = 5% (≈ p/(1−p) inflation).
- G2 fires in 100% of replicates.

**Design issues:** see design_issues.md. The PI folded them into DESIGN v0.3.1:
- hash tie-breaks;
- G2 one-sided 95% at ≈120k per class;
- the upper-bound estimator bias and the "consistent with control" label;
- the RESET fallback;
- the XL +1 variant and TT ≥ 50;
- dec_d(foreign);
- FLR adjacency stratum.

Full numbers are in metrics.json.
