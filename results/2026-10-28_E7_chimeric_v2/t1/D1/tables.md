# T1 driver: D1 (dev), runs=1, decoy seeds [1001, 1002, 1003, 1004, 1005], entrapment seeds [2001, 2002, 2003]

## IDs at own q (production searches; mean over the 5 decoy shuffles, dataset level)

| method | level | 0.1% | 1% | 5% | verdict 1% (a / b) | verdict 5% (a / b) |
|---|---|---|---|---|---|---|
| vulcan | peptidoform | 5886 | 6041 ± 15 | 6284 | supported / supported | inconclusive / inconclusive |
| vulcan | sequence | 5760 | 5908 ± 15 | 6147 | supported / supported | inconclusive / inconclusive |
| vulcan | protein | 0 | 832 ± 3 | 832 | not-evaluable / not-evaluable | supported / supported |
| vulcan_chim_v1 | peptidoform | 6662 | 7066 ± 29 | 7418 | supported / supported | supported / supported |
| vulcan_chim_v1 | sequence | 6511 | 6890 ± 29 | 7231 | supported / supported | supported / supported |
| vulcan_chim_v1 | protein | 0 | 927 ± 6 | 927 | not-evaluable / not-evaluable | supported / supported |
| vulcan_chim_v1_sep | peptidoform | 5747 | 6972 ± 49 | 7386 | supported / supported | inconclusive / inconclusive |
| vulcan_chim_v1_sep | sequence | 5626 | 6802 ± 49 | 7202 | supported / supported | inconclusive / inconclusive |
| vulcan_chim_v1_sep | protein | 0 | 927 ± 6 | 927 | not-evaluable / not-evaluable | supported / supported |
| vulcan_v2 | peptidoform | 6746 | 7165 ± 22 | 7511 | supported / supported | supported / supported |
| vulcan_v2 | sequence | 6598 | 6985 ± 20 | 7322 | supported / supported | supported / supported |
| vulcan_v2 | protein | 0 | 951 ± 9 | 951 | not-evaluable / not-evaluable | supported / supported |
| vulcan_v2_sep | peptidoform | 5786 | 6995 ± 126 | 7456 | supported / supported | inconclusive / inconclusive |
| vulcan_v2_sep | sequence | 5668 | 6823 ± 123 | 7273 | supported / supported | inconclusive / inconclusive |
| vulcan_v2_sep | protein | 0 | 951 ± 9 | 951 | not-evaluable / not-evaluable | supported / supported |
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
| vulcan_chim_v1 | peptidoform | 0.0088 | [0.0076, 0.0101] | [0.0080, 0.0096] | 0.0079 | [0.0065, 0.0093] | 0.0048 |
| vulcan_chim_v1 | sequence | 0.0090 | [0.0077, 0.0103] | [0.0081, 0.0099] | 0.0081 | [0.0067, 0.0096] | 0.0049 |
| vulcan_chim_v1 | protein | 0.0089 | [0.0053, 0.0124] | [0.0058, 0.0119] | - | - | 0.0066 |
| vulcan_chim_v1_sep | peptidoform | 0.0088 | [0.0075, 0.0102] | [0.0079, 0.0097] | 0.0083 | [0.0069, 0.0097] | 0.0047 |
| vulcan_chim_v1_sep | sequence | 0.0090 | [0.0076, 0.0103] | [0.0081, 0.0099] | 0.0086 | [0.0072, 0.0100] | 0.0048 |
| vulcan_chim_v1_sep | protein | 0.0089 | [0.0053, 0.0124] | [0.0058, 0.0119] | - | - | 0.0066 |
| vulcan_v2 | peptidoform | 0.0097 | [0.0084, 0.0109] | [0.0088, 0.0105] | 0.0089 | [0.0076, 0.0103] | 0.0051 |
| vulcan_v2 | sequence | 0.0099 | [0.0086, 0.0111] | [0.0090, 0.0108] | 0.0092 | [0.0078, 0.0106] | 0.0052 |
| vulcan_v2 | protein | 0.0077 | [0.0046, 0.0107] | [0.0050, 0.0103] | - | - | 0.0056 |
| vulcan_v2_sep | peptidoform | 0.0090 | [0.0079, 0.0102] | [0.0082, 0.0099] | 0.0084 | [0.0068, 0.0101] | 0.0047 |
| vulcan_v2_sep | sequence | 0.0092 | [0.0080, 0.0104] | [0.0083, 0.0101] | 0.0086 | [0.0069, 0.0104] | 0.0048 |
| vulcan_v2_sep | protein | 0.0077 | [0.0046, 0.0107] | [0.0050, 0.0103] | - | - | 0.0056 |
| sage | peptidoform | 0.0097 | [0.0084, 0.0110] | [0.0087, 0.0107] | 0.0088 | [0.0066, 0.0109] | 0.0050 |
| sage | sequence | 0.0098 | [0.0085, 0.0112] | [0.0089, 0.0108] | 0.0090 | [0.0068, 0.0112] | 0.0051 |
| sage | protein | 0.0084 | [0.0049, 0.0119] | [0.0055, 0.0112] | - | - | 0.0056 |
| sage+mokapot | peptidoform | 0.0094 | [0.0081, 0.0107] | [0.0085, 0.0104] | 0.0102 | [0.0081, 0.0122] | 0.0049 |
| sage+mokapot | sequence | 0.0096 | [0.0083, 0.0109] | [0.0086, 0.0105] | 0.0104 | [0.0083, 0.0125] | 0.0050 |
| sage+mokapot | protein | 0.0112 | [0.0066, 0.0159] | [0.0082, 0.0143] | - | - | 0.0074 |
| sage+pred+mokapot | peptidoform | 0.0090 | [0.0078, 0.0101] | [0.0081, 0.0099] | 0.0088 | [0.0071, 0.0106] | 0.0046 |
| sage+pred+mokapot | sequence | 0.0090 | [0.0079, 0.0102] | [0.0081, 0.0100] | 0.0091 | [0.0073, 0.0108] | 0.0047 |
| sage+pred+mokapot | protein | 0.0106 | [0.0071, 0.0141] | [0.0076, 0.0137] | - | - | 0.0073 |
| comet+mokapot | peptidoform | 0.0082 | [0.0071, 0.0093] | [0.0074, 0.0091] | 0.0088 | [0.0073, 0.0102] | 0.0044 |
| comet+mokapot | sequence | 0.0083 | [0.0071, 0.0094] | [0.0074, 0.0091] | 0.0090 | [0.0074, 0.0105] | 0.0045 |
| comet+mokapot | protein | 0.0082 | [0.0045, 0.0120] | [0.0059, 0.0106] | - | - | 0.0058 |

## Non-inferiority (ours / comparator; run-level IDs summed over runs, mean over decoy shuffles; 95% bootstrap CI over runs x decoy shuffles; pass = lower CI >= 0.97)

| comparator | level | own q 1%: ratio [CI] | counts? (comparator verdict) | matched FDP 1%: ratio [CI] | matched FDP 5%: ratio [CI] |
|---|---|---|---|---|---|
| vulcan | peptidoform | 1.158 [1.143, 1.173] PASS | yes (supported) | 1.173 [1.158, 1.176] PASS | 1.187 [1.185, 1.189] PASS |
| vulcan | sequence | 1.155 [1.140, 1.170] PASS | yes (supported) | 1.169 [1.147, 1.172] PASS | 1.184 [1.182, 1.186] PASS |
| vulcan | protein | 1.143 [1.136, 1.149] PASS | no (not-evaluable) | 1.144 [1.137, 1.151] PASS | 1.144 [1.137, 1.151] PASS |
| vulcan_chim_v1 | peptidoform | 0.990 [0.978, 1.002] PASS | yes (supported) | 1.002 [0.988, 1.003] PASS | 1.003 [1.002, 1.005] PASS |
| vulcan_chim_v1 | sequence | 0.990 [0.978, 1.002] PASS | yes (supported) | 1.002 [0.982, 1.003] PASS | 1.004 [1.003, 1.007] PASS |
| vulcan_chim_v1 | protein | 1.027 [1.020, 1.033] PASS | no (not-evaluable) | 1.027 [1.020, 1.034] PASS | 1.027 [1.020, 1.033] PASS |
| vulcan_chim_v1_sep | peptidoform | 1.003 [0.992, 1.013] PASS | yes (supported) | 1.016 [1.001, 1.018] PASS | 1.009 [1.008, 1.010] PASS |
| vulcan_chim_v1_sep | sequence | 1.003 [0.992, 1.013] PASS | yes (supported) | 1.015 [0.996, 1.018] PASS | 1.010 [1.009, 1.011] PASS |
| vulcan_chim_v1_sep | protein | 1.027 [1.020, 1.033] PASS | no (not-evaluable) | 1.027 [1.020, 1.034] PASS | 1.027 [1.020, 1.033] PASS |
| vulcan_v2 | peptidoform | 0.976 [0.962, 0.991] fail | yes (supported) | 0.990 [0.974, 0.993] PASS | 0.990 [0.989, 0.992] PASS |
| vulcan_v2 | sequence | 0.977 [0.962, 0.991] fail | yes (supported) | 0.990 [0.969, 0.993] fail | 0.991 [0.990, 0.993] PASS |
| vulcan_v2 | protein | 1.000 [1.000, 1.000] PASS | no (not-evaluable) | 1.000 [1.000, 1.000] PASS | 1.000 [1.000, 1.000] PASS |
| sage | peptidoform | 1.162 [1.145, 1.178] PASS | yes (supported) | 1.178 [1.161, 1.181] PASS | 1.186 [1.186, 1.187] PASS |
| sage | sequence | 1.159 [1.143, 1.175] PASS | yes (supported) | 1.174 [1.151, 1.178] PASS | 1.183 [1.182, 1.184] PASS |
| sage | protein | 1.148 [1.136, 1.158] PASS | no (not-evaluable) | 1.149 [1.137, 1.159] PASS | 1.149 [1.137, 1.159] PASS |
| sage+mokapot | peptidoform | 1.159 [1.143, 1.175] PASS | no (inconclusive) | 1.175 [1.157, 1.178] PASS | 1.186 [1.186, 1.187] PASS |
| sage+mokapot | sequence | 1.156 [1.140, 1.172] PASS | no (inconclusive) | 1.171 [1.147, 1.174] PASS | 1.183 [1.182, 1.185] PASS |
| sage+mokapot | protein | 1.144 [1.134, 1.154] PASS | no (not-evaluable) | 1.147 [1.141, 1.155] PASS | 1.146 [1.136, 1.155] PASS |
| sage+pred+mokapot | peptidoform | 1.156 [1.140, 1.172] PASS | yes (supported) | 1.171 [1.155, 1.175] PASS | 1.186 [1.185, 1.188] PASS |
| sage+pred+mokapot | sequence | 1.154 [1.138, 1.170] PASS | yes (supported) | 1.168 [1.145, 1.171] PASS | 1.183 [1.182, 1.185] PASS |
| sage+pred+mokapot | protein | 1.135 [1.125, 1.148] PASS | no (not-evaluable) | 1.138 [1.129, 1.150] PASS | 1.136 [1.126, 1.150] PASS |
| comet+mokapot | peptidoform | 1.157 [1.139, 1.175] PASS | yes (supported) | 1.169 [1.152, 1.174] PASS | 1.166 [1.164, 1.168] PASS |
| comet+mokapot | sequence | 1.156 [1.138, 1.174] PASS | yes (supported) | 1.166 [1.144, 1.171] PASS | 1.164 [1.162, 1.166] PASS |
| comet+mokapot | protein | 1.148 [1.139, 1.156] PASS | no (not-evaluable) | 1.148 [1.139, 1.157] PASS | 1.148 [1.139, 1.157] PASS |

Sage predict_rt check (prod, first decoy seed; Sage fits its RT model only when on): {'predict_rt_true': {'n_confident': 10827, 'frac_predicted_rt_nonzero': 0.9900774654016886, 'sd_predicted_rt': 0.20851611616966198, 'median_abs_delta_rt_model_confident': 0.043854356, 'spearman_aligned_vs_predicted_confident': 0.9426438660677244, 'pin_has_rt_features': True}, 'predict_rt_false': {'n_confident': 10823, 'frac_predicted_rt_nonzero': 0.0, 'sd_predicted_rt': 0.0, 'median_abs_delta_rt_model_confident': 0.999, 'spearman_aligned_vs_predicted_confident': nan, 'pin_has_rt_features': True}}

Sage RT on vs off (prod, first decoy seed): {'peptidoform': {'0.001': {'rt_on': 5725, 'rt_off': 5657}, '0.01': {'rt_on': 6014, 'rt_off': 6008}}, 'sequence': {'0.001': {'rt_on': 5608, 'rt_off': 5537}, '0.01': {'rt_on': 5880, 'rt_off': 5875}}}

## Driver timing: total 611 s wall, 1725 core-s

| stage | wall s | core-s (self+children) | loadavg start/end |
|---|---|---|---|
| db | 0 | 0 | 16/16 |
| evaluate | 610 | 1725 | 16/16 |
