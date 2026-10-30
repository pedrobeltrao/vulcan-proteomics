# T1 driver: D8h (dev), runs=10, decoy seeds [1001, 1002, 1003, 1004, 1005], entrapment seeds [2001, 2002, 2003]

## IDs at own q (production searches; mean over the 5 decoy shuffles, dataset level)

| method | level | 0.1% | 1% | 5% | verdict 1% (a / b) | verdict 5% (a / b) |
|---|---|---|---|---|---|---|
| vulcan_v2_sep | peptidoform | 67442 | 79307 ± 224 | 86697 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan_v2_sep | sequence | 65083 | 76287 ± 213 | 83375 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan_v2_sep | protein | 5669 | 6375 ± 33 | 6387 | inconclusive / inconclusive | supported / supported |
| sage_chim2 | peptidoform | 50359 | 62757 ± 219 | 71292 | supported / supported | supported / supported |
| sage_chim2 | sequence | 48475 | 60310 ± 203 | 68534 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage_chim2 | protein | 4998 | 5680 ± 35 | 5690 | supported / supported | supported / supported |
| sage_chim2+pred+mokapot | peptidoform | 65540 | 72123 ± 118 | 76467 | inconclusive / inconclusive | supported / supported |
| sage_chim2+pred+mokapot | sequence | 62866 | 69194 ± 108 | 73412 | inconclusive / inconclusive | supported / supported |
| sage_chim2+pred+mokapot | protein | 5748 | 6211 ± 33 | 6215 | supported / supported | supported / supported |

## Entrapment FDP estimates at 1% (mean [95% CI, convention a]; arm P = 15 replicates, arm C = 5)

| method | level | P mean | P CI (a) | P CI (b) | C mean | C CI (a) | LB(P) |
|---|---|---|---|---|---|---|---|
| vulcan_v2_sep | peptidoform | 0.0106 | [0.0101, 0.0110] | [0.0103, 0.0108] | 0.0087 | [0.0082, 0.0091] | 0.0056 |
| vulcan_v2_sep | sequence | 0.0109 | [0.0105, 0.0114] | [0.0107, 0.0112] | 0.0093 | [0.0088, 0.0098] | 0.0058 |
| vulcan_v2_sep | protein | 0.0102 | [0.0090, 0.0114] | [0.0092, 0.0111] | - | - | 0.0065 |
| sage_chim2 | peptidoform | 0.0099 | [0.0095, 0.0103] | [0.0096, 0.0102] | 0.0085 | [0.0078, 0.0092] | 0.0051 |
| sage_chim2 | sequence | 0.0103 | [0.0099, 0.0107] | [0.0100, 0.0106] | 0.0091 | [0.0084, 0.0098] | 0.0053 |
| sage_chim2 | protein | 0.0096 | [0.0080, 0.0111] | [0.0085, 0.0106] | - | - | 0.0060 |
| sage_chim2+pred+mokapot | peptidoform | 0.0097 | [0.0092, 0.0101] | [0.0094, 0.0099] | 0.0102 | [0.0097, 0.0108] | 0.0051 |
| sage_chim2+pred+mokapot | sequence | 0.0100 | [0.0096, 0.0105] | [0.0098, 0.0103] | 0.0110 | [0.0104, 0.0116] | 0.0053 |
| sage_chim2+pred+mokapot | protein | 0.0100 | [0.0086, 0.0114] | [0.0090, 0.0109] | - | - | 0.0062 |

## Non-inferiority (ours / comparator; run-level IDs summed over runs, mean over decoy shuffles; 95% bootstrap CI over runs x decoy shuffles; pass = lower CI >= 0.97)

| comparator | level | own q 1%: ratio [CI] | counts? (comparator verdict) | matched FDP 1%: ratio [CI] | matched FDP 5%: ratio [CI] |
|---|---|---|---|---|---|
| sage_chim2 | peptidoform | 1.258 [1.143, 1.343] PASS | yes (supported) | 1.251 [1.137, 1.329] PASS | 1.218 [1.112, 1.286] PASS |
| sage_chim2 | sequence | 1.257 [1.139, 1.342] PASS | no (inconclusive) | 1.250 [1.134, 1.330] PASS | 1.217 [1.110, 1.286] PASS |
| sage_chim2 | protein | 1.068 [1.011, 1.114] PASS | yes (supported) | 1.070 [1.014, 1.114] PASS | 1.070 [1.014, 1.115] PASS |
| sage_chim2+pred+mokapot | peptidoform | 1.149 [1.074, 1.196] PASS | no (inconclusive) | 1.150 [1.069, 1.203] PASS | 1.173 [1.085, 1.231] PASS |
| sage_chim2+pred+mokapot | sequence | 1.150 [1.075, 1.197] PASS | no (inconclusive) | 1.150 [1.070, 1.203] PASS | 1.173 [1.084, 1.230] PASS |
| sage_chim2+pred+mokapot | protein | 0.997 [0.949, 1.037] fail | yes (supported) | 0.999 [0.951, 1.039] fail | 0.999 [0.951, 1.039] fail |

Sage predict_rt check (prod, first decoy seed; Sage fits its RT model only when on): {'predict_rt_true': {'n_confident': 160937, 'frac_predicted_rt_nonzero': 0.9993603054787376, 'sd_predicted_rt': 0.17929871775528924, 'median_abs_delta_rt_model_confident': 0.08180553, 'spearman_aligned_vs_predicted_confident': 0.852867597490656, 'pin_has_rt_features': True}, 'predict_rt_false': None}

## Driver timing: total 6626 s wall, 7038 core-s

| stage | wall s | core-s (self+children) | loadavg start/end |
|---|---|---|---|
| db | 0 | 0 | 8/8 |
| evaluate | 6625 | 7038 | 8/8 |
