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
| vulcan_v2_sep | peptidoform | 25793 | 29484 ± 28 | 31651 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan_v2_sep | sequence | 24455 | 27746 ± 27 | 29798 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan_v2_sep | protein | 6571 | 6979 ± 8 | 6982 | inconclusive / inconclusive | supported / supported |
| sage_chim2 | peptidoform | 21094 | 23944 ± 43 | 25795 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage_chim2 | sequence | 19999 | 22644 ± 42 | 24420 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage_chim2 | protein | 5884 | 6339 ± 13 | 6342 | inconclusive / inconclusive | supported / supported |
| sage_chim2+pred+mokapot | peptidoform | 24003 | 25292 ± 20 | 26514 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage_chim2+pred+mokapot | sequence | 22650 | 23850 ± 20 | 25046 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage_chim2+pred+mokapot | protein | 6323 | 6553 ± 6 | 6555 | inconclusive / inconclusive | supported / supported |

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
| vulcan_v2_sep | peptidoform | 0.0107 | [0.0101, 0.0114] | [0.0103, 0.0112] | 0.0077 | [0.0070, 0.0084] | 0.0055 |
| vulcan_v2_sep | sequence | 0.0113 | [0.0106, 0.0120] | [0.0108, 0.0118] | 0.0085 | [0.0077, 0.0092] | 0.0058 |
| vulcan_v2_sep | protein | 0.0106 | [0.0093, 0.0120] | [0.0097, 0.0115] | - | - | 0.0068 |
| sage_chim2 | peptidoform | 0.0109 | [0.0102, 0.0115] | [0.0104, 0.0113] | 0.0087 | [0.0079, 0.0096] | 0.0055 |
| sage_chim2 | sequence | 0.0114 | [0.0107, 0.0120] | [0.0109, 0.0119] | 0.0095 | [0.0086, 0.0104] | 0.0058 |
| sage_chim2 | protein | 0.0121 | [0.0106, 0.0136] | [0.0111, 0.0130] | - | - | 0.0077 |
| sage_chim2+pred+mokapot | peptidoform | 0.0102 | [0.0096, 0.0108] | [0.0098, 0.0107] | 0.0087 | [0.0080, 0.0094] | 0.0052 |
| sage_chim2+pred+mokapot | sequence | 0.0108 | [0.0101, 0.0114] | [0.0103, 0.0112] | 0.0095 | [0.0087, 0.0103] | 0.0055 |
| sage_chim2+pred+mokapot | protein | 0.0113 | [0.0099, 0.0126] | [0.0103, 0.0122] | - | - | 0.0074 |

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
| vulcan_v2_sep | peptidoform | 0.817 [0.813, 0.820] fail | no (inconclusive) | 0.818 [0.815, 0.819] fail | 0.797 [0.795, 0.800] fail |
| vulcan_v2_sep | sequence | 0.818 [0.814, 0.820] fail | no (inconclusive) | 0.819 [0.817, 0.821] fail | 0.799 [0.797, 0.803] fail |
| vulcan_v2_sep | protein | 0.870 [0.858, 0.881] fail | no (inconclusive) | 0.870 [0.857, 0.882] fail | 0.870 [0.858, 0.882] fail |
| sage_chim2 | peptidoform | 1.006 [0.999, 1.012] PASS | no (inconclusive) | 1.010 [1.004, 1.014] PASS | 0.982 [0.980, 0.984] PASS |
| sage_chim2 | sequence | 1.003 [0.996, 1.009] PASS | no (inconclusive) | 1.007 [1.002, 1.011] PASS | 0.981 [0.978, 0.983] PASS |
| sage_chim2 | protein | 0.949 [0.936, 0.961] fail | no (inconclusive) | 0.952 [0.938, 0.963] fail | 0.949 [0.936, 0.961] fail |
| sage_chim2+pred+mokapot | peptidoform | 0.954 [0.951, 0.957] fail | no (inconclusive) | 0.954 [0.952, 0.957] fail | 0.952 [0.950, 0.955] fail |
| sage_chim2+pred+mokapot | sequence | 0.954 [0.952, 0.957] fail | no (inconclusive) | 0.954 [0.952, 0.957] fail | 0.952 [0.950, 0.955] fail |
| sage_chim2+pred+mokapot | protein | 0.922 [0.910, 0.932] fail | no (inconclusive) | 0.923 [0.910, 0.933] fail | 0.922 [0.911, 0.932] fail |

Sage predict_rt check (prod, first decoy seed; Sage fits its RT model only when on): {'predict_rt_true': {'n_confident': 43967, 'frac_predicted_rt_nonzero': 0.9823871268400998, 'sd_predicted_rt': 0.2359752486561844, 'median_abs_delta_rt_model_confident': 0.05510366, 'spearman_aligned_vs_predicted_confident': 0.9210243368620381, 'pin_has_rt_features': True}, 'predict_rt_false': None}

## Driver timing: total 10865 s wall, 53732 core-s

| stage | wall s | core-s (self+children) | loadavg start/end |
|---|---|---|---|
| db | 0 | 0 | 17/17 |
| ours_cohort_ingest_and_prod_search | 22 | 35 | 17/20 |
| ours_replicate_searches | 78 | 76 | 20/16 |
| ours_features_and_scoring | 336 | 384 | 16/38 |
| chim_ms1 | 17 | 29 | 38/38 |
| chim_replicate_searches | 90 | 86 | 38/39 |
| chim_features_and_scoring | 630 | 704 | 39/36 |
| cmp_mzml | 223 | 522 | 36/21 |
| cmp_sage_searches | 1461 | 28147 | 21/44 |
| cmp_sage_chim_searches | 1039 | 12692 | 44/15 |
| cmp_sage_chim_pred_features | 717 | 2301 | 15/12 |
| cmp_sage_pred_features | 570 | 1886 | 12/7 |
| cmp_comet_searches_wait | 0 | 0 | 7/7 |
| cmp_mokapot_wait | 91 | 1252 | 7/20 |
| evaluate_a3_secondary | 2260 | 2254 | 6/13 |
| evaluate | 5592 | 5617 | 20/13 |
