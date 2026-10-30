# T1 driver: D0 (dev), runs=2, decoy seeds [1001, 1002, 1003, 1004, 1005], entrapment seeds [2001, 2002, 2003]

## IDs at own q (production searches; mean over the 5 decoy shuffles, dataset level)

| method | level | 0.1% | 1% | 5% | verdict 1% (a / b) | verdict 5% (a / b) |
|---|---|---|---|---|---|---|
| vulcan | peptidoform | 23262 | 24328 ± 36 | 25486 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan | sequence | 22014 | 22981 ± 34 | 24107 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan | protein | 6216 | 6429 ± 15 | 6430 | inconclusive / inconclusive | supported / supported |
| vulcan_chim_v1_sep | peptidoform | 26596 | 29056 ± 47 | 30745 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan_chim_v1_sep | sequence | 25121 | 27311 ± 44 | 28942 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan_chim_v1_sep | protein | 6664 | 6951 ± 13 | 6955 | inconclusive / inconclusive | supported / supported |

## Entrapment FDP estimates at 1% (mean [95% CI, convention a]; arm P = 15 replicates, arm C = 5)

| method | level | P mean | P CI (a) | P CI (b) | C mean | C CI (a) | LB(P) |
|---|---|---|---|---|---|---|---|
| vulcan | peptidoform | 0.0101 | [0.0094, 0.0107] | [0.0096, 0.0105] | 0.0083 | [0.0074, 0.0093] | 0.0051 |
| vulcan | sequence | 0.0105 | [0.0099, 0.0112] | [0.0100, 0.0110] | 0.0091 | [0.0081, 0.0101] | 0.0054 |
| vulcan | protein | 0.0122 | [0.0110, 0.0135] | [0.0112, 0.0133] | - | - | 0.0078 |
| vulcan_chim_v1_sep | peptidoform | 0.0100 | [0.0095, 0.0106] | [0.0096, 0.0105] | 0.0088 | [0.0079, 0.0096] | 0.0051 |
| vulcan_chim_v1_sep | sequence | 0.0106 | [0.0100, 0.0111] | [0.0101, 0.0110] | 0.0096 | [0.0087, 0.0105] | 0.0054 |
| vulcan_chim_v1_sep | protein | 0.0129 | [0.0117, 0.0140] | [0.0119, 0.0139] | - | - | 0.0084 |

## Non-inferiority (ours / comparator; run-level IDs summed over runs, mean over decoy shuffles; 95% bootstrap CI over runs x decoy shuffles; pass = lower CI >= 0.97)

| comparator | level | own q 1%: ratio [CI] | counts? (comparator verdict) | matched FDP 1%: ratio [CI] | matched FDP 5%: ratio [CI] |
|---|---|---|---|---|---|
| vulcan | peptidoform | 1.205 [1.203, 1.207] PASS | no (inconclusive) | 1.204 [1.203, 1.208] PASS | 1.216 [1.215, 1.217] PASS |
| vulcan | sequence | 1.202 [1.200, 1.204] PASS | no (inconclusive) | 1.201 [1.199, 1.204] PASS | 1.212 [1.212, 1.214] PASS |
| vulcan | protein | 1.142 [1.128, 1.156] PASS | no (inconclusive) | 1.140 [1.127, 1.156] PASS | 1.142 [1.128, 1.156] PASS |

Sage predict_rt check (prod, first decoy seed; Sage fits its RT model only when on): {'predict_rt_true': {'n_confident': 43967, 'frac_predicted_rt_nonzero': 0.9823871268400998, 'sd_predicted_rt': 0.2359752486561844, 'median_abs_delta_rt_model_confident': 0.05510366, 'spearman_aligned_vs_predicted_confident': 0.9210243368620381, 'pin_has_rt_features': True}, 'predict_rt_false': {'n_confident': 43605, 'frac_predicted_rt_nonzero': 0.0, 'sd_predicted_rt': 0.0, 'median_abs_delta_rt_model_confident': 0.999, 'spearman_aligned_vs_predicted_confident': nan, 'pin_has_rt_features': True}}

Sage RT on vs off (prod, first decoy seed): {'peptidoform': {'0.001': {'rt_on': 20865, 'rt_off': 20621}, '0.01': {'rt_on': 23732, 'rt_off': 23432}}, 'sequence': {'0.001': {'rt_on': 19778, 'rt_off': 19530}, '0.01': {'rt_on': 22431, 'rt_off': 22142}}}

## Driver timing: total 1627 s wall, 1671 core-s

| stage | wall s | core-s (self+children) | loadavg start/end |
|---|---|---|---|
| db | 1 | 0 | 27/27 |
| evaluate | 1625 | 1671 | 27/64 |
