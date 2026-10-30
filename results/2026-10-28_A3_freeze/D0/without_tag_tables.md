# T1 driver: D0 (dev), runs=2, decoy seeds [1001, 1002, 1003, 1004, 1005], entrapment seeds [2001, 2002, 2003]

## IDs at own q (production searches; mean over the 5 decoy shuffles, dataset level)

| method | level | 0.1% | 1% | 5% | verdict 1% (a / b) | verdict 5% (a / b) |
|---|---|---|---|---|---|---|
| vulcan | peptidoform | 23262 | 24328 ± 36 | 25486 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan | sequence | 22014 | 22981 ± 34 | 24107 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan | protein | 6216 | 6429 ± 15 | 6430 | inconclusive / inconclusive | supported / supported |
| sage | peptidoform | 21253 | 23727 ± 44 | 25358 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage | sequence | 20140 | 22434 ± 40 | 23994 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage | protein | 5923 | 6321 ± 12 | 6323 | inconclusive / inconclusive | supported / supported |
| sage+mokapot | peptidoform | 21344 | 23827 ± 17 | 25406 | inconclusive / inconclusive | supported / supported |
| sage+mokapot | sequence | 20205 | 22517 ± 20 | 24031 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage+mokapot | protein | 6005 | 6336 ± 23 | 6337 | inconclusive / inconclusive | supported / supported |
| sage+pred+mokapot | peptidoform | 23539 | 24575 ± 29 | 25655 | inconclusive / inconclusive | supported / supported |
| sage+pred+mokapot | sequence | 22216 | 23179 ± 31 | 24242 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage+pred+mokapot | protein | 6249 | 6465 ± 14 | 6465 | inconclusive / inconclusive | supported / supported |
| comet+mokapot | peptidoform | 20381 | 23443 ± 122 | 25162 | supported / supported | supported / supported |
| comet+mokapot | sequence | 19277 | 22163 ± 118 | 23817 | inconclusive / inconclusive | inconclusive / inconclusive |
| comet+mokapot | protein | 5772 | 6267 ± 16 | 6273 | supported / supported | supported / supported |

## Entrapment FDP estimates at 1% (mean [95% CI, convention a]; arm P = 15 replicates, arm C = 5)

| method | level | P mean | P CI (a) | P CI (b) | C mean | C CI (a) | LB(P) |
|---|---|---|---|---|---|---|---|
| vulcan | peptidoform | 0.0101 | [0.0094, 0.0107] | [0.0096, 0.0105] | 0.0083 | [0.0074, 0.0093] | 0.0051 |
| vulcan | sequence | 0.0105 | [0.0099, 0.0112] | [0.0100, 0.0110] | 0.0091 | [0.0081, 0.0101] | 0.0054 |
| vulcan | protein | 0.0122 | [0.0110, 0.0135] | [0.0112, 0.0133] | - | - | 0.0078 |
| sage | peptidoform | 0.0110 | [0.0103, 0.0116] | [0.0105, 0.0115] | 0.0092 | [0.0083, 0.0101] | 0.0056 |
| sage | sequence | 0.0115 | [0.0108, 0.0122] | [0.0110, 0.0120] | 0.0101 | [0.0090, 0.0111] | 0.0059 |
| sage | protein | 0.0123 | [0.0107, 0.0138] | [0.0113, 0.0132] | - | - | 0.0078 |
| sage+mokapot | peptidoform | 0.0114 | [0.0107, 0.0121] | [0.0109, 0.0119] | 0.0089 | [0.0079, 0.0099] | 0.0058 |
| sage+mokapot | sequence | 0.0120 | [0.0113, 0.0126] | [0.0114, 0.0125] | 0.0097 | [0.0086, 0.0108] | 0.0061 |
| sage+mokapot | protein | 0.0121 | [0.0105, 0.0136] | [0.0110, 0.0131] | - | - | 0.0074 |
| sage+pred+mokapot | peptidoform | 0.0104 | [0.0098, 0.0111] | [0.0099, 0.0109] | 0.0093 | [0.0084, 0.0101] | 0.0053 |
| sage+pred+mokapot | sequence | 0.0109 | [0.0103, 0.0116] | [0.0105, 0.0114] | 0.0102 | [0.0092, 0.0111] | 0.0056 |
| sage+pred+mokapot | protein | 0.0113 | [0.0099, 0.0127] | [0.0103, 0.0123] | - | - | 0.0072 |
| comet+mokapot | peptidoform | 0.0094 | [0.0087, 0.0101] | [0.0089, 0.0098] | 0.0093 | [0.0082, 0.0104] | 0.0047 |
| comet+mokapot | sequence | 0.0098 | [0.0091, 0.0106] | [0.0094, 0.0103] | 0.0102 | [0.0090, 0.0114] | 0.0050 |
| comet+mokapot | protein | 0.0090 | [0.0078, 0.0101] | [0.0081, 0.0098] | - | - | 0.0058 |

## Non-inferiority (ours / comparator; run-level IDs summed over runs, mean over decoy shuffles; 95% bootstrap CI over runs x decoy shuffles; pass = lower CI >= 0.97)

| comparator | level | own q 1%: ratio [CI] | counts? (comparator verdict) | matched FDP 1%: ratio [CI] | matched FDP 5%: ratio [CI] |
|---|---|---|---|---|---|
| sage | peptidoform | 1.025 [1.022, 1.028] PASS | no (inconclusive) | 1.028 [1.026, 1.029] PASS | 1.006 [1.002, 1.010] PASS |
| sage | sequence | 1.024 [1.021, 1.026] PASS | no (inconclusive) | 1.028 [1.025, 1.029] PASS | 1.006 [1.003, 1.010] PASS |
| sage | protein | 1.013 [1.011, 1.016] PASS | no (inconclusive) | 1.014 [1.012, 1.016] PASS | 1.014 [1.011, 1.016] PASS |
| sage+mokapot | peptidoform | 1.021 [1.015, 1.025] PASS | no (inconclusive) | 1.026 [1.022, 1.027] PASS | 1.001 [0.999, 1.006] PASS |
| sage+mokapot | sequence | 1.020 [1.015, 1.024] PASS | no (inconclusive) | 1.026 [1.022, 1.026] PASS | 1.002 [0.999, 1.006] PASS |
| sage+mokapot | protein | 1.011 [1.008, 1.015] PASS | no (inconclusive) | 1.011 [1.008, 1.015] PASS | 1.011 [1.008, 1.015] PASS |
| sage+pred+mokapot | peptidoform | 0.990 [0.988, 0.993] PASS | no (inconclusive) | 0.990 [0.989, 0.992] PASS | 0.992 [0.990, 0.994] PASS |
| sage+pred+mokapot | sequence | 0.992 [0.989, 0.995] PASS | no (inconclusive) | 0.992 [0.990, 0.993] PASS | 0.993 [0.991, 0.995] PASS |
| sage+pred+mokapot | protein | 0.995 [0.994, 0.997] PASS | no (inconclusive) | 0.996 [0.993, 0.997] PASS | 0.996 [0.994, 0.997] PASS |
| comet+mokapot | peptidoform | 1.037 [1.025, 1.048] PASS | yes (supported) | 1.034 [1.026, 1.039] PASS | 1.008 [1.002, 1.013] PASS |
| comet+mokapot | sequence | 1.035 [1.024, 1.046] PASS | no (inconclusive) | 1.033 [1.026, 1.038] PASS | 1.009 [1.002, 1.014] PASS |
| comet+mokapot | protein | 0.997 [0.989, 1.004] PASS | yes (supported) | 0.997 [0.988, 1.004] PASS | 0.997 [0.989, 1.004] PASS |

Sage predict_rt check (prod, first decoy seed; Sage fits its RT model only when on): {'predict_rt_true': {'n_confident': 43967, 'frac_predicted_rt_nonzero': 0.9823871268400998, 'sd_predicted_rt': 0.2359752486561844, 'median_abs_delta_rt_model_confident': 0.05510366, 'spearman_aligned_vs_predicted_confident': 0.9210243368620381, 'pin_has_rt_features': True}, 'predict_rt_false': None}

## Driver timing: total 2003 s wall, 2680 core-s

| stage | wall s | core-s (self+children) | loadavg start/end |
|---|---|---|---|
| db | 1 | 0 | 22/22 |
| evaluate | 2002 | 2680 | 22/7 |
