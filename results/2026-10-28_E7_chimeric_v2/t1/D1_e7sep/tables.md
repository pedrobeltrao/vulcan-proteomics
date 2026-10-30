# T1 driver: D1 (dev), runs=1, decoy seeds [1001, 1002, 1003, 1004, 1005], entrapment seeds [2001, 2002, 2003]

## IDs at own q (production searches; mean over the 5 decoy shuffles, dataset level)

| method | level | 0.1% | 1% | 5% | verdict 1% (a / b) | verdict 5% (a / b) |
|---|---|---|---|---|---|---|
| vulcan | peptidoform | 5886 | 6041 ± 15 | 6284 | supported / supported | inconclusive / inconclusive |
| vulcan | sequence | 5760 | 5908 ± 15 | 6147 | supported / supported | inconclusive / inconclusive |
| vulcan | protein | 0 | 832 ± 3 | 832 | not-evaluable / not-evaluable | supported / supported |
| vulcan_chim_v1_sep | peptidoform | 5747 | 6972 ± 49 | 7386 | supported / supported | inconclusive / inconclusive |
| vulcan_chim_v1_sep | sequence | 5626 | 6802 ± 49 | 7202 | supported / supported | inconclusive / inconclusive |
| vulcan_chim_v1_sep | protein | 0 | 927 ± 6 | 927 | not-evaluable / not-evaluable | supported / supported |

## Entrapment FDP estimates at 1% (mean [95% CI, convention a]; arm P = 15 replicates, arm C = 5)

| method | level | P mean | P CI (a) | P CI (b) | C mean | C CI (a) | LB(P) |
|---|---|---|---|---|---|---|---|
| vulcan | peptidoform | 0.0087 | [0.0073, 0.0101] | [0.0077, 0.0096] | 0.0083 | [0.0065, 0.0100] | 0.0047 |
| vulcan | sequence | 0.0088 | [0.0074, 0.0102] | [0.0079, 0.0097] | 0.0085 | [0.0067, 0.0103] | 0.0048 |
| vulcan | protein | 0.0072 | [0.0036, 0.0107] | [0.0040, 0.0104] | - | - | 0.0050 |
| vulcan_chim_v1_sep | peptidoform | 0.0088 | [0.0075, 0.0102] | [0.0079, 0.0097] | 0.0083 | [0.0069, 0.0097] | 0.0047 |
| vulcan_chim_v1_sep | sequence | 0.0090 | [0.0076, 0.0103] | [0.0081, 0.0099] | 0.0086 | [0.0072, 0.0100] | 0.0048 |
| vulcan_chim_v1_sep | protein | 0.0089 | [0.0053, 0.0124] | [0.0058, 0.0119] | - | - | 0.0066 |

## Non-inferiority (ours / comparator; run-level IDs summed over runs, mean over decoy shuffles; 95% bootstrap CI over runs x decoy shuffles; pass = lower CI >= 0.97)

| comparator | level | own q 1%: ratio [CI] | counts? (comparator verdict) | matched FDP 1%: ratio [CI] | matched FDP 5%: ratio [CI] |
|---|---|---|---|---|---|
| vulcan | peptidoform | 1.154 [1.148, 1.160] PASS | yes (supported) | 1.155 [1.150, 1.160] PASS | 1.176 [1.174, 1.178] PASS |
| vulcan | sequence | 1.151 [1.145, 1.157] PASS | yes (supported) | 1.152 [1.146, 1.157] PASS | 1.172 [1.170, 1.175] PASS |
| vulcan | protein | 1.113 [1.110, 1.117] PASS | no (not-evaluable) | 1.114 [1.111, 1.118] PASS | 1.114 [1.111, 1.118] PASS |

Sage predict_rt check (prod, first decoy seed; Sage fits its RT model only when on): {'predict_rt_true': {'n_confident': 10827, 'frac_predicted_rt_nonzero': 0.9900774654016886, 'sd_predicted_rt': 0.20851611616966198, 'median_abs_delta_rt_model_confident': 0.043854356, 'spearman_aligned_vs_predicted_confident': 0.9426438660677244, 'pin_has_rt_features': True}, 'predict_rt_false': {'n_confident': 10823, 'frac_predicted_rt_nonzero': 0.0, 'sd_predicted_rt': 0.0, 'median_abs_delta_rt_model_confident': 0.999, 'spearman_aligned_vs_predicted_confident': nan, 'pin_has_rt_features': True}}

Sage RT on vs off (prod, first decoy seed): {'peptidoform': {'0.001': {'rt_on': 5725, 'rt_off': 5657}, '0.01': {'rt_on': 6014, 'rt_off': 6008}}, 'sequence': {'0.001': {'rt_on': 5608, 'rt_off': 5537}, '0.01': {'rt_on': 5880, 'rt_off': 5875}}}

## Driver timing: total 214 s wall, 450 core-s

| stage | wall s | core-s (self+children) | loadavg start/end |
|---|---|---|---|
| db | 1 | 0 | 27/27 |
| evaluate | 212 | 450 | 27/29 |
