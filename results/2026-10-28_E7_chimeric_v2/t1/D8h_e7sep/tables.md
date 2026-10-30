# T1 driver: D8h (dev), runs=10, decoy seeds [1001, 1002, 1003, 1004, 1005], entrapment seeds [2001, 2002, 2003]

## IDs at own q (production searches; mean over the 5 decoy shuffles, dataset level)

| method | level | 0.1% | 1% | 5% | verdict 1% (a / b) | verdict 5% (a / b) |
|---|---|---|---|---|---|---|
| vulcan | peptidoform | 59638 | 67392 ± 38 | 72732 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan | sequence | 57444 | 64710 ± 40 | 69894 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan | protein | 5316 | 5984 ± 34 | 5990 | supported / supported | supported / supported |
| vulcan_chim_v1_sep | peptidoform | 63955 | 75293 ± 163 | 82592 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan_chim_v1_sep | sequence | 61695 | 72343 ± 157 | 79373 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan_chim_v1_sep | protein | 5497 | 6261 ± 54 | 6272 | inconclusive / inconclusive | supported / supported |

## Entrapment FDP estimates at 1% (mean [95% CI, convention a]; arm P = 15 replicates, arm C = 5)

| method | level | P mean | P CI (a) | P CI (b) | C mean | C CI (a) | LB(P) |
|---|---|---|---|---|---|---|---|
| vulcan | peptidoform | 0.0102 | [0.0098, 0.0106] | [0.0099, 0.0105] | 0.0104 | [0.0098, 0.0111] | 0.0053 |
| vulcan | sequence | 0.0106 | [0.0101, 0.0110] | [0.0103, 0.0109] | 0.0112 | [0.0105, 0.0119] | 0.0055 |
| vulcan | protein | 0.0100 | [0.0087, 0.0113] | [0.0090, 0.0109] | - | - | 0.0062 |
| vulcan_chim_v1_sep | peptidoform | 0.0101 | [0.0097, 0.0105] | [0.0098, 0.0104] | 0.0103 | [0.0097, 0.0108] | 0.0053 |
| vulcan_chim_v1_sep | sequence | 0.0104 | [0.0101, 0.0108] | [0.0102, 0.0107] | 0.0110 | [0.0104, 0.0116] | 0.0055 |
| vulcan_chim_v1_sep | protein | 0.0101 | [0.0086, 0.0116] | [0.0092, 0.0110] | - | - | 0.0062 |

## Non-inferiority (ours / comparator; run-level IDs summed over runs, mean over decoy shuffles; 95% bootstrap CI over runs x decoy shuffles; pass = lower CI >= 0.97)

| comparator | level | own q 1%: ratio [CI] | counts? (comparator verdict) | matched FDP 1%: ratio [CI] | matched FDP 5%: ratio [CI] |
|---|---|---|---|---|---|
| vulcan | peptidoform | 1.190 [1.113, 1.243] PASS | no (inconclusive) | 1.191 [1.113, 1.244] PASS | 1.203 [1.121, 1.257] PASS |
| vulcan | sequence | 1.193 [1.116, 1.245] PASS | no (inconclusive) | 1.194 [1.115, 1.246] PASS | 1.205 [1.123, 1.259] PASS |
| vulcan | protein | 1.197 [1.137, 1.261] PASS | yes (supported) | 1.198 [1.138, 1.261] PASS | 1.198 [1.138, 1.262] PASS |

Sage predict_rt check (prod, first decoy seed; Sage fits its RT model only when on): {'predict_rt_true': {'n_confident': 160937, 'frac_predicted_rt_nonzero': 0.9993603054787376, 'sd_predicted_rt': 0.17929871775528924, 'median_abs_delta_rt_model_confident': 0.08180553, 'spearman_aligned_vs_predicted_confident': 0.852867597490656, 'pin_has_rt_features': True}, 'predict_rt_false': None}

## Driver timing: total 4114 s wall, 4373 core-s

| stage | wall s | core-s (self+children) | loadavg start/end |
|---|---|---|---|
| db | 0 | 0 | 11/11 |
| evaluate | 4111 | 4373 | 11/10 |
