# T1 driver: D8h (dev), runs=10, decoy seeds [1001, 1002, 1003, 1004, 1005], entrapment seeds [2001, 2002, 2003]

## IDs at own q (production searches; mean over the 5 decoy shuffles, dataset level)

| method | level | 0.1% | 1% | 5% | verdict 1% (a / b) | verdict 5% (a / b) |
|---|---|---|---|---|---|---|
| vulcan | peptidoform | 59638 | 67392 ± 38 | 72732 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan | sequence | 57444 | 64710 ± 40 | 69894 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan | protein | 5316 | 5984 ± 34 | 5990 | supported / supported | supported / supported |
| vulcan_v2_sep | peptidoform | 67442 | 79307 ± 224 | 86697 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan_v2_sep | sequence | 65083 | 76287 ± 213 | 83375 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan_v2_sep | protein | 5669 | 6375 ± 33 | 6387 | inconclusive / inconclusive | supported / supported |
| sage+pred+mokapot | peptidoform | 64486 | 70194 ± 93 | 73754 | supported / supported | supported / supported |
| sage+pred+mokapot | sequence | 61861 | 67361 ± 91 | 70841 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage+pred+mokapot | protein | 5748 | 6149 ± 28 | 6153 | supported / supported | supported / supported |

## Entrapment FDP estimates at 1% (mean [95% CI, convention a]; arm P = 15 replicates, arm C = 5)

| method | level | P mean | P CI (a) | P CI (b) | C mean | C CI (a) | LB(P) |
|---|---|---|---|---|---|---|---|
| vulcan | peptidoform | 0.0102 | [0.0098, 0.0106] | [0.0099, 0.0105] | 0.0104 | [0.0098, 0.0111] | 0.0053 |
| vulcan | sequence | 0.0106 | [0.0101, 0.0110] | [0.0103, 0.0109] | 0.0112 | [0.0105, 0.0119] | 0.0055 |
| vulcan | protein | 0.0100 | [0.0087, 0.0113] | [0.0090, 0.0109] | - | - | 0.0062 |
| vulcan_v2_sep | peptidoform | 0.0106 | [0.0101, 0.0110] | [0.0103, 0.0108] | 0.0087 | [0.0082, 0.0091] | 0.0056 |
| vulcan_v2_sep | sequence | 0.0109 | [0.0105, 0.0114] | [0.0107, 0.0112] | 0.0093 | [0.0088, 0.0098] | 0.0058 |
| vulcan_v2_sep | protein | 0.0102 | [0.0090, 0.0114] | [0.0092, 0.0111] | - | - | 0.0065 |
| sage+pred+mokapot | peptidoform | 0.0097 | [0.0093, 0.0102] | [0.0095, 0.0100] | 0.0099 | [0.0094, 0.0104] | 0.0051 |
| sage+pred+mokapot | sequence | 0.0101 | [0.0096, 0.0105] | [0.0098, 0.0104] | 0.0106 | [0.0101, 0.0112] | 0.0053 |
| sage+pred+mokapot | protein | 0.0098 | [0.0086, 0.0111] | [0.0089, 0.0108] | - | - | 0.0062 |

## Non-inferiority (ours / comparator; run-level IDs summed over runs, mean over decoy shuffles; 95% bootstrap CI over runs x decoy shuffles; pass = lower CI >= 0.97)

| comparator | level | own q 1%: ratio [CI] | counts? (comparator verdict) | matched FDP 1%: ratio [CI] | matched FDP 5%: ratio [CI] |
|---|---|---|---|---|---|
| vulcan | peptidoform | 1.252 [1.153, 1.318] PASS | no (inconclusive) | 1.252 [1.149, 1.321] PASS | 1.263 [1.160, 1.329] PASS |
| vulcan | sequence | 1.256 [1.155, 1.320] PASS | no (inconclusive) | 1.256 [1.152, 1.324] PASS | 1.266 [1.162, 1.331] PASS |
| vulcan | protein | 1.217 [1.153, 1.284] PASS | yes (supported) | 1.221 [1.156, 1.282] PASS | 1.221 [1.156, 1.287] PASS |
| sage+pred+mokapot | peptidoform | 1.220 [1.140, 1.272] PASS | yes (supported) | 1.222 [1.138, 1.279] PASS | 1.253 [1.159, 1.314] PASS |
| sage+pred+mokapot | sequence | 1.225 [1.144, 1.275] PASS | no (inconclusive) | 1.226 [1.142, 1.281] PASS | 1.256 [1.161, 1.316] PASS |
| sage+pred+mokapot | protein | 1.190 [1.131, 1.252] PASS | yes (supported) | 1.192 [1.133, 1.252] PASS | 1.192 [1.133, 1.254] PASS |

Sage predict_rt check (prod, first decoy seed; Sage fits its RT model only when on): {'predict_rt_true': {'n_confident': 160937, 'frac_predicted_rt_nonzero': 0.9993603054787376, 'sd_predicted_rt': 0.17929871775528924, 'median_abs_delta_rt_model_confident': 0.08180553, 'spearman_aligned_vs_predicted_confident': 0.852867597490656, 'pin_has_rt_features': True}, 'predict_rt_false': None}

## Driver timing: total 5947 s wall, 6370 core-s

| stage | wall s | core-s (self+children) | loadavg start/end |
|---|---|---|---|
| db | 1 | 0 | 11/11 |
| evaluate | 5944 | 6370 | 11/8 |
