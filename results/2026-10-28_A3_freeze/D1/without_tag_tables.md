# T1 driver: D1 (dev), runs=1, decoy seeds [1001, 1002, 1003, 1004, 1005], entrapment seeds [2001, 2002, 2003]

## IDs at own q (production searches; mean over the 5 decoy shuffles, dataset level)

| method | level | 0.1% | 1% | 5% | verdict 1% (a / b) | verdict 5% (a / b) |
|---|---|---|---|---|---|---|
| vulcan | peptidoform | 5886 | 6041 ± 15 | 6284 | supported / supported | inconclusive / inconclusive |
| vulcan | sequence | 5760 | 5908 ± 15 | 6147 | supported / supported | inconclusive / inconclusive |
| vulcan | protein | 0 | 832 ± 3 | 832 | not-evaluable / not-evaluable | supported / supported |
| sage | peptidoform | 5730 | 6021 ± 5 | 6273 | supported / supported | inconclusive / inconclusive |
| sage | sequence | 5612 | 5887 ± 5 | 6134 | supported / supported | inconclusive / inconclusive |
| sage | protein | 0 | 829 ± 4 | 829 | not-evaluable / not-evaluable | supported / supported |
| sage+mokapot | peptidoform | 5595 | 6036 ± 4 | 6274 | inconclusive / inconclusive | supported / supported |
| sage+mokapot | sequence | 5476 | 5903 ± 3 | 6135 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage+mokapot | protein | 0 | 831 ± 4 | 831 | not-evaluable / not-evaluable | supported / supported |
| sage+pred+mokapot | peptidoform | 5568 | 6050 ± 3 | 6272 | supported / supported | supported / supported |
| sage+pred+mokapot | sequence | 5446 | 5915 ± 3 | 6133 | supported / supported | supported / supported |
| sage+pred+mokapot | protein | 0 | 839 ± 8 | 839 | not-evaluable / not-evaluable | supported / supported |
| comet+mokapot | peptidoform | 5763 | 6044 ± 14 | 6342 | supported / supported | inconclusive / inconclusive |
| comet+mokapot | sequence | 5631 | 5905 ± 14 | 6198 | supported / supported | inconclusive / inconclusive |
| comet+mokapot | protein | 0 | 829 ± 5 | 829 | not-evaluable / not-evaluable | supported / supported |

## Entrapment FDP estimates at 1% (mean [95% CI, convention a]; arm P = 15 replicates, arm C = 5)

| method | level | P mean | P CI (a) | P CI (b) | C mean | C CI (a) | LB(P) |
|---|---|---|---|---|---|---|---|
| vulcan | peptidoform | 0.0087 | [0.0073, 0.0101] | [0.0077, 0.0096] | 0.0083 | [0.0065, 0.0100] | 0.0047 |
| vulcan | sequence | 0.0088 | [0.0074, 0.0102] | [0.0079, 0.0097] | 0.0085 | [0.0067, 0.0103] | 0.0048 |
| vulcan | protein | 0.0072 | [0.0036, 0.0107] | [0.0040, 0.0104] | - | - | 0.0050 |
| sage | peptidoform | 0.0097 | [0.0084, 0.0110] | [0.0087, 0.0107] | 0.0088 | [0.0066, 0.0109] | 0.0050 |
| sage | sequence | 0.0098 | [0.0085, 0.0112] | [0.0089, 0.0108] | 0.0090 | [0.0068, 0.0112] | 0.0051 |
| sage | protein | 0.0084 | [0.0049, 0.0119] | [0.0055, 0.0112] | - | - | 0.0056 |
| sage+mokapot | peptidoform | 0.0094 | [0.0081, 0.0107] | [0.0085, 0.0104] | 0.0101 | [0.0081, 0.0121] | 0.0049 |
| sage+mokapot | sequence | 0.0096 | [0.0083, 0.0109] | [0.0086, 0.0105] | 0.0103 | [0.0083, 0.0124] | 0.0050 |
| sage+mokapot | protein | 0.0112 | [0.0066, 0.0159] | [0.0082, 0.0143] | - | - | 0.0074 |
| sage+pred+mokapot | peptidoform | 0.0091 | [0.0080, 0.0102] | [0.0082, 0.0100] | 0.0087 | [0.0067, 0.0107] | 0.0047 |
| sage+pred+mokapot | sequence | 0.0092 | [0.0081, 0.0103] | [0.0083, 0.0101] | 0.0089 | [0.0070, 0.0109] | 0.0048 |
| sage+pred+mokapot | protein | 0.0105 | [0.0069, 0.0140] | [0.0074, 0.0135] | - | - | 0.0073 |
| comet+mokapot | peptidoform | 0.0082 | [0.0071, 0.0093] | [0.0074, 0.0091] | 0.0088 | [0.0073, 0.0102] | 0.0044 |
| comet+mokapot | sequence | 0.0083 | [0.0071, 0.0094] | [0.0074, 0.0091] | 0.0090 | [0.0074, 0.0105] | 0.0045 |
| comet+mokapot | protein | 0.0082 | [0.0045, 0.0120] | [0.0059, 0.0106] | - | - | 0.0058 |

## Non-inferiority (ours / comparator; run-level IDs summed over runs, mean over decoy shuffles; 95% bootstrap CI over runs x decoy shuffles; pass = lower CI >= 0.97)

| comparator | level | own q 1%: ratio [CI] | counts? (comparator verdict) | matched FDP 1%: ratio [CI] | matched FDP 5%: ratio [CI] |
|---|---|---|---|---|---|
| sage | peptidoform | 1.003 [1.001, 1.005] PASS | yes (supported) | 1.004 [1.003, 1.005] PASS | 1.000 [0.998, 1.001] PASS |
| sage | sequence | 1.004 [1.002, 1.005] PASS | yes (supported) | 1.004 [1.003, 1.006] PASS | 1.000 [0.998, 1.002] PASS |
| sage | protein | 1.004 [1.000, 1.008] PASS | no (not-evaluable) | 1.004 [1.000, 1.008] PASS | 1.004 [1.000, 1.008] PASS |
| sage+mokapot | peptidoform | 1.001 [0.999, 1.003] PASS | no (inconclusive) | 1.001 [0.999, 1.002] PASS | 0.999 [0.998, 1.001] PASS |
| sage+mokapot | sequence | 1.001 [0.999, 1.003] PASS | no (inconclusive) | 1.002 [0.999, 1.003] PASS | 1.000 [0.998, 1.002] PASS |
| sage+mokapot | protein | 1.001 [0.998, 1.005] PASS | no (not-evaluable) | 1.002 [0.999, 1.006] PASS | 1.001 [0.998, 1.005] PASS |
| sage+pred+mokapot | peptidoform | 0.998 [0.996, 1.001] PASS | yes (supported) | 0.998 [0.997, 1.000] PASS | 1.000 [0.999, 1.001] PASS |
| sage+pred+mokapot | sequence | 0.999 [0.996, 1.001] PASS | yes (supported) | 0.999 [0.997, 1.000] PASS | 1.000 [0.999, 1.001] PASS |
| sage+pred+mokapot | protein | 0.993 [0.984, 1.002] PASS | no (not-evaluable) | 0.994 [0.989, 1.002] PASS | 0.993 [0.984, 1.002] PASS |
| comet+mokapot | peptidoform | 0.999 [0.996, 1.003] PASS | yes (supported) | 0.996 [0.994, 0.999] PASS | 0.982 [0.981, 0.984] PASS |
| comet+mokapot | sequence | 1.001 [0.997, 1.004] PASS | yes (supported) | 0.997 [0.995, 1.000] PASS | 0.984 [0.982, 0.986] PASS |
| comet+mokapot | protein | 1.004 [1.000, 1.008] PASS | no (not-evaluable) | 1.003 [0.999, 1.008] PASS | 1.003 [0.999, 1.008] PASS |

Sage predict_rt check (prod, first decoy seed; Sage fits its RT model only when on): {'predict_rt_true': {'n_confident': 10827, 'frac_predicted_rt_nonzero': 0.9900774654016886, 'sd_predicted_rt': 0.20851611616966198, 'median_abs_delta_rt_model_confident': 0.043854356, 'spearman_aligned_vs_predicted_confident': 0.9426438660677244, 'pin_has_rt_features': True}, 'predict_rt_false': None}

## Driver timing: total 543 s wall, 3478 core-s

| stage | wall s | core-s (self+children) | loadavg start/end |
|---|---|---|---|
| db | 1 | 0 | 41/42 |
| cmp_mzml | 0 | 0 | 42/42 |
| cmp_sage_searches | 0 | 0 | 42/42 |
| cmp_sage_pred_features | 114 | 1566 | 42/61 |
| cmp_comet_searches_wait | 0 | 0 | 61/61 |
| cmp_mokapot_wait | 93 | 1564 | 61/71 |
| evaluate | 334 | 348 | 71/42 |
