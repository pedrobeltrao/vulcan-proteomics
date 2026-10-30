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
| sage+pred+mokapot | peptidoform | 5493 | 6051 ± 2 | 6272 | supported / supported | supported / supported |
| sage+pred+mokapot | sequence | 5373 | 5916 ± 2 | 6133 | supported / supported | supported / supported |
| sage+pred+mokapot | protein | 0 | 838 ± 7 | 838 | not-evaluable / not-evaluable | supported / supported |
| comet+mokapot | peptidoform | 5763 | 6044 ± 14 | 6342 | supported / supported | inconclusive / inconclusive |
| comet+mokapot | sequence | 5631 | 5905 ± 14 | 6198 | supported / supported | inconclusive / inconclusive |
| comet+mokapot | protein | 0 | 829 ± 5 | 829 | not-evaluable / not-evaluable | supported / supported |
| vulcan_v2_sep | peptidoform | 5786 | 6995 ± 126 | 7456 | supported / supported | inconclusive / inconclusive |
| vulcan_v2_sep | sequence | 5668 | 6823 ± 123 | 7273 | supported / supported | inconclusive / inconclusive |
| vulcan_v2_sep | protein | 0 | 951 ± 9 | 951 | not-evaluable / not-evaluable | supported / supported |
| sage_chim2 | peptidoform | 5712 | 6106 ± 15 | 6444 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage_chim2 | sequence | 5594 | 5969 ± 14 | 6295 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage_chim2 | protein | 0 | 836 ± 6 | 836 | not-evaluable / not-evaluable | supported / supported |
| sage_chim2+pred+mokapot | peptidoform | 5995 | 6223 ± 11 | 6485 | supported / supported | supported / supported |
| sage_chim2+pred+mokapot | sequence | 5856 | 6073 ± 10 | 6331 | supported / supported | supported / supported |
| sage_chim2+pred+mokapot | protein | 0 | 871 ± 12 | 871 | not-evaluable / not-evaluable | supported / supported |

## Entrapment FDP estimates at 1% (mean [95% CI, convention a]; arm P = 15 replicates, arm C = 5)

| method | level | P mean | P CI (a) | P CI (b) | C mean | C CI (a) | LB(P) |
|---|---|---|---|---|---|---|---|
| vulcan | peptidoform | 0.0087 | [0.0073, 0.0101] | [0.0077, 0.0096] | 0.0083 | [0.0065, 0.0100] | 0.0047 |
| vulcan | sequence | 0.0088 | [0.0074, 0.0102] | [0.0079, 0.0097] | 0.0085 | [0.0067, 0.0103] | 0.0048 |
| vulcan | protein | 0.0072 | [0.0036, 0.0107] | [0.0040, 0.0104] | - | - | 0.0050 |
| sage | peptidoform | 0.0097 | [0.0084, 0.0110] | [0.0087, 0.0107] | 0.0088 | [0.0066, 0.0109] | 0.0050 |
| sage | sequence | 0.0098 | [0.0085, 0.0112] | [0.0089, 0.0108] | 0.0090 | [0.0068, 0.0112] | 0.0051 |
| sage | protein | 0.0084 | [0.0049, 0.0119] | [0.0055, 0.0112] | - | - | 0.0056 |
| sage+mokapot | peptidoform | 0.0094 | [0.0081, 0.0107] | [0.0084, 0.0104] | 0.0100 | [0.0080, 0.0120] | 0.0049 |
| sage+mokapot | sequence | 0.0095 | [0.0082, 0.0108] | [0.0086, 0.0105] | 0.0103 | [0.0082, 0.0123] | 0.0050 |
| sage+mokapot | protein | 0.0112 | [0.0066, 0.0159] | [0.0082, 0.0143] | - | - | 0.0074 |
| sage+pred+mokapot | peptidoform | 0.0089 | [0.0078, 0.0101] | [0.0080, 0.0099] | 0.0088 | [0.0071, 0.0106] | 0.0046 |
| sage+pred+mokapot | sequence | 0.0090 | [0.0078, 0.0102] | [0.0081, 0.0099] | 0.0091 | [0.0073, 0.0108] | 0.0047 |
| sage+pred+mokapot | protein | 0.0106 | [0.0071, 0.0142] | [0.0076, 0.0137] | - | - | 0.0073 |
| comet+mokapot | peptidoform | 0.0082 | [0.0071, 0.0093] | [0.0074, 0.0091] | 0.0088 | [0.0073, 0.0102] | 0.0044 |
| comet+mokapot | sequence | 0.0083 | [0.0071, 0.0094] | [0.0074, 0.0091] | 0.0090 | [0.0074, 0.0105] | 0.0045 |
| comet+mokapot | protein | 0.0082 | [0.0045, 0.0120] | [0.0059, 0.0106] | - | - | 0.0058 |
| vulcan_v2_sep | peptidoform | 0.0090 | [0.0079, 0.0102] | [0.0082, 0.0099] | 0.0084 | [0.0068, 0.0101] | 0.0047 |
| vulcan_v2_sep | sequence | 0.0092 | [0.0080, 0.0104] | [0.0083, 0.0101] | 0.0086 | [0.0069, 0.0104] | 0.0048 |
| vulcan_v2_sep | protein | 0.0077 | [0.0046, 0.0107] | [0.0050, 0.0103] | - | - | 0.0056 |
| sage_chim2 | peptidoform | 0.0102 | [0.0088, 0.0116] | [0.0092, 0.0112] | 0.0088 | [0.0068, 0.0107] | 0.0053 |
| sage_chim2 | sequence | 0.0103 | [0.0089, 0.0118] | [0.0094, 0.0113] | 0.0090 | [0.0070, 0.0110] | 0.0054 |
| sage_chim2 | protein | 0.0082 | [0.0047, 0.0117] | [0.0055, 0.0108] | - | - | 0.0054 |
| sage_chim2+pred+mokapot | peptidoform | 0.0092 | [0.0079, 0.0105] | [0.0083, 0.0101] | 0.0080 | [0.0057, 0.0104] | 0.0050 |
| sage_chim2+pred+mokapot | sequence | 0.0093 | [0.0080, 0.0106] | [0.0084, 0.0102] | 0.0083 | [0.0059, 0.0106] | 0.0052 |
| sage_chim2+pred+mokapot | protein | 0.0122 | [0.0085, 0.0159] | [0.0093, 0.0151] | - | - | 0.0087 |

## Non-inferiority (ours / comparator; run-level IDs summed over runs, mean over decoy shuffles; 95% bootstrap CI over runs x decoy shuffles; pass = lower CI >= 0.97)

| comparator | level | own q 1%: ratio [CI] | counts? (comparator verdict) | matched FDP 1%: ratio [CI] | matched FDP 5%: ratio [CI] |
|---|---|---|---|---|---|
| sage | peptidoform | 1.003 [1.001, 1.005] PASS | yes (supported) | 1.004 [1.003, 1.005] PASS | 1.000 [0.998, 1.001] PASS |
| sage | sequence | 1.004 [1.002, 1.005] PASS | yes (supported) | 1.004 [1.003, 1.006] PASS | 1.000 [0.998, 1.002] PASS |
| sage | protein | 1.004 [1.000, 1.008] PASS | no (not-evaluable) | 1.004 [1.000, 1.008] PASS | 1.004 [1.000, 1.008] PASS |
| sage+mokapot | peptidoform | 1.001 [0.999, 1.003] PASS | no (inconclusive) | 1.001 [0.999, 1.002] PASS | 0.999 [0.998, 1.001] PASS |
| sage+mokapot | sequence | 1.001 [0.999, 1.003] PASS | no (inconclusive) | 1.002 [0.999, 1.003] PASS | 1.000 [0.998, 1.002] PASS |
| sage+mokapot | protein | 1.001 [0.998, 1.005] PASS | no (not-evaluable) | 1.002 [0.999, 1.006] PASS | 1.001 [0.998, 1.005] PASS |
| sage+pred+mokapot | peptidoform | 0.998 [0.996, 1.000] PASS | yes (supported) | 0.998 [0.996, 0.999] PASS | 1.000 [0.999, 1.001] PASS |
| sage+pred+mokapot | sequence | 0.999 [0.996, 1.001] PASS | yes (supported) | 0.998 [0.997, 1.000] PASS | 1.000 [0.999, 1.001] PASS |
| sage+pred+mokapot | protein | 0.993 [0.985, 1.002] PASS | no (not-evaluable) | 0.995 [0.990, 1.002] PASS | 0.993 [0.985, 1.002] PASS |
| comet+mokapot | peptidoform | 0.999 [0.996, 1.003] PASS | yes (supported) | 0.996 [0.994, 0.999] PASS | 0.982 [0.981, 0.984] PASS |
| comet+mokapot | sequence | 1.001 [0.997, 1.004] PASS | yes (supported) | 0.997 [0.995, 1.000] PASS | 0.984 [0.982, 0.986] PASS |
| comet+mokapot | protein | 1.004 [1.000, 1.008] PASS | no (not-evaluable) | 1.003 [0.999, 1.008] PASS | 1.003 [0.999, 1.008] PASS |
| vulcan_v2_sep | peptidoform | 0.864 [0.853, 0.875] fail | yes (supported) | 0.852 [0.850, 0.864] fail | 0.843 [0.841, 0.844] fail |
| vulcan_v2_sep | sequence | 0.866 [0.855, 0.877] fail | yes (supported) | 0.855 [0.853, 0.871] fail | 0.845 [0.843, 0.846] fail |
| vulcan_v2_sep | protein | 0.875 [0.870, 0.880] fail | no (not-evaluable) | 0.874 [0.869, 0.879] fail | 0.874 [0.869, 0.879] fail |
| sage_chim2 | peptidoform | 0.989 [0.987, 0.992] PASS | no (inconclusive) | 0.991 [0.989, 0.992] PASS | 0.974 [0.972, 0.975] PASS |
| sage_chim2 | sequence | 0.990 [0.987, 0.992] PASS | no (inconclusive) | 0.992 [0.990, 0.993] PASS | 0.975 [0.974, 0.977] PASS |
| sage_chim2 | protein | 0.996 [0.991, 1.000] PASS | no (not-evaluable) | 0.995 [0.991, 1.000] PASS | 0.995 [0.991, 1.000] PASS |
| sage_chim2+pred+mokapot | peptidoform | 0.971 [0.969, 0.972] fail | yes (supported) | 0.970 [0.967, 0.973] fail | 0.965 [0.962, 0.968] fail |
| sage_chim2+pred+mokapot | sequence | 0.973 [0.971, 0.975] PASS | yes (supported) | 0.972 [0.969, 0.975] fail | 0.967 [0.963, 0.970] fail |
| sage_chim2+pred+mokapot | protein | 0.955 [0.942, 0.965] fail | no (not-evaluable) | 0.960 [0.955, 0.966] fail | 0.955 [0.941, 0.964] fail |

Sage predict_rt check (prod, first decoy seed; Sage fits its RT model only when on): {'predict_rt_true': {'n_confident': 10827, 'frac_predicted_rt_nonzero': 0.9900774654016886, 'sd_predicted_rt': 0.20851611616966198, 'median_abs_delta_rt_model_confident': 0.043854356, 'spearman_aligned_vs_predicted_confident': 0.9426438660677244, 'pin_has_rt_features': True}, 'predict_rt_false': None}

## Driver timing: total 2235 s wall, 12009 core-s

| stage | wall s | core-s (self+children) | loadavg start/end |
|---|---|---|---|
| db | 0 | 0 | 19/19 |
| ours_cohort_ingest_and_prod_search | 8 | 8 | 19/18 |
| ours_replicate_searches | 15 | 14 | 18/16 |
| ours_features_and_scoring | 49 | 55 | 16/11 |
| chim_ms1 | 12 | 13 | 11/11 |
| chim_replicate_searches | 12 | 12 | 11/10 |
| chim_features_and_scoring | 53 | 62 | 10/9 |
| cmp_mzml | 52 | 37 | 9/10 |
| cmp_sage_searches | 393 | 5276 | 10/47 |
| cmp_sage_chim_searches | 249 | 2258 | 47/40 |
| cmp_sage_chim_pred_features | 107 | 996 | 40/29 |
| cmp_sage_pred_features | 98 | 1028 | 29/61 |
| cmp_comet_searches_wait | 0 | 0 | 61/61 |
| cmp_mokapot_wait | 69 | 1087 | 61/54 |
| evaluate_a3_secondary | 233 | 246 | 9/8 |
| evaluate_a5_additions | 363 | 374 | 8/13 |
| evaluate | 1117 | 1164 | 54/13 |
