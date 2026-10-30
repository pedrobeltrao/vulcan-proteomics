# T1 driver: D0 (dev), runs=2, decoy seeds [1001, 1002, 1003, 1004, 1005], entrapment seeds [2001, 2002, 2003]

## IDs at own q (production searches; mean over the 5 decoy shuffles, dataset level)

| method | level | 0.1% | 1% | 5% | verdict 1% (a / b) | verdict 5% (a / b) |
|---|---|---|---|---|---|---|
| vulcan | peptidoform | 23262 | 24328 ± 36 | 25486 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan | sequence | 22014 | 22981 ± 34 | 24107 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan | protein | 6216 | 6429 ± 15 | 6430 | inconclusive / inconclusive | supported / supported |
| sage | peptidoform | 21253 | 23727 ± 44 | 25358 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage | sequence | 20140 | 22434 ± 40 | 23994 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage | protein | 5923 | 6320 ± 12 | 6323 | inconclusive / inconclusive | supported / supported |
| sage+mokapot | peptidoform | 21389 | 23822 ± 28 | 25407 | inconclusive / inconclusive | supported / supported |
| sage+mokapot | sequence | 20247 | 22512 ± 29 | 24032 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage+mokapot | protein | 6006 | 6337 ± 20 | 6339 | inconclusive / inconclusive | supported / supported |
| sage+pred+mokapot | peptidoform | 23530 | 24576 ± 31 | 25656 | inconclusive / inconclusive | supported / supported |
| sage+pred+mokapot | sequence | 22207 | 23180 ± 33 | 24243 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage+pred+mokapot | protein | 6255 | 6465 ± 15 | 6468 | inconclusive / inconclusive | supported / supported |
| comet+mokapot | peptidoform | 20381 | 23443 ± 122 | 25162 | supported / supported | supported / supported |
| comet+mokapot | sequence | 19277 | 22163 ± 118 | 23817 | inconclusive / inconclusive | inconclusive / inconclusive |
| comet+mokapot | protein | 5772 | 6267 ± 16 | 6273 | supported / supported | supported / supported |

## Entrapment FDP estimates at 1% (mean [95% CI, convention a]; arm P = 15 replicates, arm C = 5)

| method | level | P mean | P CI (a) | P CI (b) | C mean | C CI (a) | LB(P) |
|---|---|---|---|---|---|---|---|
| vulcan | peptidoform | 0.0101 | [0.0094, 0.0107] | [0.0096, 0.0105] | 0.0083 | [0.0074, 0.0093] | 0.0051 |
| vulcan | sequence | 0.0105 | [0.0099, 0.0112] | [0.0100, 0.0110] | 0.0091 | [0.0081, 0.0101] | 0.0054 |
| vulcan | protein | 0.0122 | [0.0110, 0.0135] | [0.0112, 0.0133] | - | - | 0.0078 |
| sage | peptidoform | 0.0110 | [0.0103, 0.0117] | [0.0105, 0.0115] | 0.0092 | [0.0083, 0.0101] | 0.0056 |
| sage | sequence | 0.0115 | [0.0108, 0.0122] | [0.0110, 0.0120] | 0.0101 | [0.0090, 0.0111] | 0.0059 |
| sage | protein | 0.0123 | [0.0107, 0.0138] | [0.0113, 0.0132] | - | - | 0.0078 |
| sage+mokapot | peptidoform | 0.0114 | [0.0108, 0.0121] | [0.0109, 0.0119] | 0.0088 | [0.0078, 0.0098] | 0.0058 |
| sage+mokapot | sequence | 0.0120 | [0.0113, 0.0126] | [0.0114, 0.0125] | 0.0096 | [0.0085, 0.0107] | 0.0061 |
| sage+mokapot | protein | 0.0122 | [0.0107, 0.0137] | [0.0112, 0.0132] | - | - | 0.0075 |
| sage+pred+mokapot | peptidoform | 0.0105 | [0.0098, 0.0111] | [0.0100, 0.0109] | 0.0093 | [0.0085, 0.0101] | 0.0053 |
| sage+pred+mokapot | sequence | 0.0110 | [0.0103, 0.0117] | [0.0105, 0.0115] | 0.0102 | [0.0093, 0.0111] | 0.0056 |
| sage+pred+mokapot | protein | 0.0113 | [0.0098, 0.0127] | [0.0102, 0.0123] | - | - | 0.0071 |
| comet+mokapot | peptidoform | 0.0094 | [0.0087, 0.0101] | [0.0089, 0.0098] | 0.0093 | [0.0082, 0.0104] | 0.0047 |
| comet+mokapot | sequence | 0.0098 | [0.0091, 0.0106] | [0.0094, 0.0103] | 0.0102 | [0.0090, 0.0114] | 0.0050 |
| comet+mokapot | protein | 0.0090 | [0.0078, 0.0101] | [0.0081, 0.0098] | - | - | 0.0058 |

## Non-inferiority (ours / comparator; run-level IDs summed over runs, mean over decoy shuffles; 95% bootstrap CI over runs x decoy shuffles; pass = lower CI >= 0.97)

| comparator | level | own q 1%: ratio [CI] | counts? (comparator verdict) | matched FDP 1%: ratio [CI] | matched FDP 5%: ratio [CI] |
|---|---|---|---|---|---|
| sage | peptidoform | 1.025 [1.022, 1.028] PASS | no (inconclusive) | 1.028 [1.026, 1.029] PASS | 1.006 [1.002, 1.010] PASS |
| sage | sequence | 1.024 [1.021, 1.026] PASS | no (inconclusive) | 1.027 [1.025, 1.029] PASS | 1.006 [1.003, 1.010] PASS |
| sage | protein | 1.013 [1.011, 1.016] PASS | no (inconclusive) | 1.014 [1.012, 1.016] PASS | 1.014 [1.011, 1.016] PASS |
| sage+mokapot | peptidoform | 1.021 [1.016, 1.025] PASS | no (inconclusive) | 1.026 [1.023, 1.027] PASS | 1.001 [0.998, 1.006] PASS |
| sage+mokapot | sequence | 1.020 [1.015, 1.024] PASS | no (inconclusive) | 1.026 [1.023, 1.027] PASS | 1.002 [0.999, 1.007] PASS |
| sage+mokapot | protein | 1.011 [1.008, 1.015] PASS | no (inconclusive) | 1.011 [1.009, 1.015] PASS | 1.011 [1.008, 1.015] PASS |
| sage+pred+mokapot | peptidoform | 0.990 [0.988, 0.994] PASS | no (inconclusive) | 0.990 [0.989, 0.992] PASS | 0.992 [0.990, 0.994] PASS |
| sage+pred+mokapot | sequence | 0.992 [0.989, 0.995] PASS | no (inconclusive) | 0.991 [0.990, 0.993] PASS | 0.993 [0.991, 0.995] PASS |
| sage+pred+mokapot | protein | 0.995 [0.993, 0.997] PASS | no (inconclusive) | 0.995 [0.993, 0.997] PASS | 0.995 [0.994, 0.997] PASS |
| comet+mokapot | peptidoform | 1.037 [1.025, 1.048] PASS | yes (supported) | 1.034 [1.026, 1.039] PASS | 1.008 [1.002, 1.013] PASS |
| comet+mokapot | sequence | 1.035 [1.024, 1.046] PASS | no (inconclusive) | 1.033 [1.026, 1.038] PASS | 1.009 [1.002, 1.014] PASS |
| comet+mokapot | protein | 0.997 [0.989, 1.004] PASS | yes (supported) | 0.997 [0.988, 1.004] PASS | 0.997 [0.989, 1.004] PASS |

Sage predict_rt check (prod, first decoy seed; Sage fits its RT model only when on): {'predict_rt_true': {'n_confident': 43967, 'frac_predicted_rt_nonzero': 0.9823871268400998, 'sd_predicted_rt': 0.2359752486561844, 'median_abs_delta_rt_model_confident': 0.05510366, 'spearman_aligned_vs_predicted_confident': 0.9210243368620381, 'pin_has_rt_features': True}, 'predict_rt_false': {'n_confident': 43605, 'frac_predicted_rt_nonzero': 0.0, 'sd_predicted_rt': 0.0, 'median_abs_delta_rt_model_confident': 0.999, 'spearman_aligned_vs_predicted_confident': nan, 'pin_has_rt_features': True}}

Sage RT on vs off (prod, first decoy seed): {'peptidoform': {'0.001': {'rt_on': 20865, 'rt_off': 20621}, '0.01': {'rt_on': 23732, 'rt_off': 23432}}, 'sequence': {'0.001': {'rt_on': 19778, 'rt_off': 19530}, '0.01': {'rt_on': 22431, 'rt_off': 22142}}}

## Driver timing: total 2858 s wall, 2601 core-s

| stage | wall s | core-s (self+children) | loadavg start/end |
|---|---|---|---|
| db | 0 | 0 | 44/44 |
| evaluate | 2852 | 2601 | 44/87 |
