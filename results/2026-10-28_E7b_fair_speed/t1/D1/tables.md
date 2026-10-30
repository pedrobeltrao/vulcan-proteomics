# T1 driver: D1 (dev), runs=1, decoy seeds [1001, 1002, 1003, 1004, 1005], entrapment seeds [2001, 2002, 2003]

## IDs at own q (production searches; mean over the 5 decoy shuffles, dataset level)

| method | level | 0.1% | 1% | 5% | verdict 1% (a / b) | verdict 5% (a / b) |
|---|---|---|---|---|---|---|
| vulcan_v2_sep | peptidoform | 5786 | 6995 ± 126 | 7456 | supported / supported | inconclusive / inconclusive |
| vulcan_v2_sep | sequence | 5668 | 6823 ± 123 | 7273 | supported / supported | inconclusive / inconclusive |
| vulcan_v2_sep | protein | 0 | 951 ± 9 | 951 | not-evaluable / not-evaluable | supported / supported |
| sage+pred+mokapot | peptidoform | 5568 | 6050 ± 3 | 6272 | supported / supported | supported / supported |
| sage+pred+mokapot | sequence | 5446 | 5915 ± 3 | 6133 | supported / supported | supported / supported |
| sage+pred+mokapot | protein | 0 | 839 ± 8 | 839 | not-evaluable / not-evaluable | supported / supported |
| sage_chim2 | peptidoform | 5712 | 6106 ± 15 | 6444 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage_chim2 | sequence | 5594 | 5969 ± 14 | 6295 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage_chim2 | protein | 0 | 836 ± 6 | 836 | not-evaluable / not-evaluable | supported / supported |
| sage_chim2+pred+mokapot | peptidoform | 6000 | 6225 ± 10 | 6484 | supported / supported | supported / supported |
| sage_chim2+pred+mokapot | sequence | 5861 | 6075 ± 10 | 6330 | supported / supported | supported / supported |
| sage_chim2+pred+mokapot | protein | 0 | 871 ± 11 | 871 | not-evaluable / not-evaluable | supported / supported |
| sage_chim5 | peptidoform | 5699 | 6120 ± 16 | 6449 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage_chim5 | sequence | 5585 | 5983 ± 16 | 6299 | inconclusive / inconclusive | inconclusive / violated |
| sage_chim5 | protein | 0 | 842 ± 6 | 842 | not-evaluable / not-evaluable | supported / supported |
| sage_chim5+pred+mokapot | peptidoform | 6020 | 6232 ± 11 | 6486 | supported / supported | supported / supported |
| sage_chim5+pred+mokapot | sequence | 5880 | 6082 ± 10 | 6332 | supported / supported | supported / supported |
| sage_chim5+pred+mokapot | protein | 0 | 871 ± 17 | 871 | not-evaluable / not-evaluable | supported / supported |

## Entrapment FDP estimates at 1% (mean [95% CI, convention a]; arm P = 15 replicates, arm C = 5)

| method | level | P mean | P CI (a) | P CI (b) | C mean | C CI (a) | LB(P) |
|---|---|---|---|---|---|---|---|
| vulcan_v2_sep | peptidoform | 0.0090 | [0.0079, 0.0102] | [0.0082, 0.0099] | 0.0084 | [0.0068, 0.0101] | 0.0047 |
| vulcan_v2_sep | sequence | 0.0092 | [0.0080, 0.0104] | [0.0083, 0.0101] | 0.0086 | [0.0069, 0.0104] | 0.0048 |
| vulcan_v2_sep | protein | 0.0077 | [0.0046, 0.0107] | [0.0050, 0.0103] | - | - | 0.0056 |
| sage+pred+mokapot | peptidoform | 0.0090 | [0.0078, 0.0101] | [0.0081, 0.0099] | 0.0088 | [0.0071, 0.0106] | 0.0046 |
| sage+pred+mokapot | sequence | 0.0090 | [0.0079, 0.0102] | [0.0081, 0.0100] | 0.0091 | [0.0073, 0.0108] | 0.0047 |
| sage+pred+mokapot | protein | 0.0106 | [0.0071, 0.0141] | [0.0076, 0.0137] | - | - | 0.0073 |
| sage_chim2 | peptidoform | 0.0102 | [0.0088, 0.0116] | [0.0092, 0.0112] | 0.0088 | [0.0068, 0.0107] | 0.0053 |
| sage_chim2 | sequence | 0.0103 | [0.0089, 0.0118] | [0.0094, 0.0113] | 0.0090 | [0.0070, 0.0110] | 0.0054 |
| sage_chim2 | protein | 0.0082 | [0.0047, 0.0117] | [0.0055, 0.0108] | - | - | 0.0054 |
| sage_chim2+pred+mokapot | peptidoform | 0.0092 | [0.0079, 0.0104] | [0.0083, 0.0101] | 0.0085 | [0.0060, 0.0109] | 0.0050 |
| sage_chim2+pred+mokapot | sequence | 0.0093 | [0.0080, 0.0105] | [0.0084, 0.0102] | 0.0087 | [0.0062, 0.0112] | 0.0051 |
| sage_chim2+pred+mokapot | protein | 0.0122 | [0.0083, 0.0162] | [0.0090, 0.0154] | - | - | 0.0086 |
| sage_chim5 | peptidoform | 0.0107 | [0.0094, 0.0119] | [0.0097, 0.0116] | 0.0084 | [0.0060, 0.0107] | 0.0056 |
| sage_chim5 | sequence | 0.0108 | [0.0096, 0.0121] | [0.0099, 0.0118] | 0.0086 | [0.0062, 0.0110] | 0.0057 |
| sage_chim5 | protein | 0.0083 | [0.0047, 0.0119] | [0.0057, 0.0108] | - | - | 0.0055 |
| sage_chim5+pred+mokapot | peptidoform | 0.0091 | [0.0079, 0.0103] | [0.0082, 0.0100] | 0.0076 | [0.0058, 0.0094] | 0.0050 |
| sage_chim5+pred+mokapot | sequence | 0.0092 | [0.0080, 0.0104] | [0.0083, 0.0101] | 0.0078 | [0.0059, 0.0096] | 0.0051 |
| sage_chim5+pred+mokapot | protein | 0.0110 | [0.0076, 0.0144] | [0.0084, 0.0136] | - | - | 0.0082 |

## Non-inferiority (ours / comparator; run-level IDs summed over runs, mean over decoy shuffles; 95% bootstrap CI over runs x decoy shuffles; pass = lower CI >= 0.97)

| comparator | level | own q 1%: ratio [CI] | counts? (comparator verdict) | matched FDP 1%: ratio [CI] | matched FDP 5%: ratio [CI] |
|---|---|---|---|---|---|
| sage+pred+mokapot | peptidoform | 1.156 [1.140, 1.172] PASS | yes (supported) | 1.171 [1.155, 1.175] PASS | 1.186 [1.185, 1.188] PASS |
| sage+pred+mokapot | sequence | 1.154 [1.138, 1.170] PASS | yes (supported) | 1.168 [1.145, 1.171] PASS | 1.183 [1.182, 1.185] PASS |
| sage+pred+mokapot | protein | 1.135 [1.125, 1.148] PASS | no (not-evaluable) | 1.138 [1.129, 1.150] PASS | 1.136 [1.126, 1.150] PASS |
| sage_chim2 | peptidoform | 1.146 [1.129, 1.163] PASS | no (inconclusive) | 1.163 [1.145, 1.166] PASS | 1.155 [1.154, 1.157] PASS |
| sage_chim2 | sequence | 1.143 [1.126, 1.160] PASS | no (inconclusive) | 1.160 [1.136, 1.163] PASS | 1.154 [1.152, 1.156] PASS |
| sage_chim2 | protein | 1.138 [1.128, 1.149] PASS | no (not-evaluable) | 1.139 [1.129, 1.149] PASS | 1.139 [1.129, 1.149] PASS |
| sage_chim2+pred+mokapot | peptidoform | 1.124 [1.109, 1.138] PASS | yes (supported) | 1.137 [1.120, 1.140] PASS | 1.145 [1.142, 1.148] PASS |
| sage_chim2+pred+mokapot | sequence | 1.123 [1.109, 1.138] PASS | yes (supported) | 1.136 [1.114, 1.139] PASS | 1.144 [1.141, 1.147] PASS |
| sage_chim2+pred+mokapot | protein | 1.092 [1.084, 1.100] PASS | no (not-evaluable) | 1.100 [1.094, 1.106] PASS | 1.093 [1.084, 1.101] PASS |
| sage_chim5 | peptidoform | 1.143 [1.127, 1.160] PASS | no (inconclusive) | 1.161 [1.144, 1.164] PASS | 1.156 [1.155, 1.158] PASS |
| sage_chim5 | sequence | 1.140 [1.124, 1.157] PASS | no (inconclusive) | 1.158 [1.134, 1.161] PASS | 1.155 [1.154, 1.157] PASS |
| sage_chim5 | protein | 1.129 [1.119, 1.140] PASS | no (not-evaluable) | 1.131 [1.120, 1.141] PASS | 1.131 [1.120, 1.141] PASS |
| sage_chim5+pred+mokapot | peptidoform | 1.123 [1.107, 1.138] PASS | yes (supported) | 1.138 [1.121, 1.141] PASS | 1.144 [1.142, 1.146] PASS |
| sage_chim5+pred+mokapot | sequence | 1.122 [1.107, 1.137] PASS | yes (supported) | 1.136 [1.114, 1.140] PASS | 1.144 [1.141, 1.146] PASS |
| sage_chim5+pred+mokapot | protein | 1.092 [1.078, 1.103] PASS | no (not-evaluable) | 1.095 [1.084, 1.106] PASS | 1.092 [1.077, 1.104] PASS |

Sage predict_rt check (prod, first decoy seed; Sage fits its RT model only when on): {'predict_rt_true': {'n_confident': 10827, 'frac_predicted_rt_nonzero': 0.9900774654016886, 'sd_predicted_rt': 0.20851611616966198, 'median_abs_delta_rt_model_confident': 0.043854356, 'spearman_aligned_vs_predicted_confident': 0.9426438660677244, 'pin_has_rt_features': True}, 'predict_rt_false': {'n_confident': 10823, 'frac_predicted_rt_nonzero': 0.0, 'sd_predicted_rt': 0.0, 'median_abs_delta_rt_model_confident': 0.999, 'spearman_aligned_vs_predicted_confident': nan, 'pin_has_rt_features': True}}

Sage RT on vs off (prod, first decoy seed): {'peptidoform': {'0.001': {'rt_on': 5725, 'rt_off': 5657}, '0.01': {'rt_on': 6014, 'rt_off': 6008}}, 'sequence': {'0.001': {'rt_on': 5608, 'rt_off': 5537}, '0.01': {'rt_on': 5880, 'rt_off': 5875}}}

## Driver timing: total 506 s wall, 876 core-s

| stage | wall s | core-s (self+children) | loadavg start/end |
|---|---|---|---|
| db | 0 | 0 | 45/45 |
| evaluate | 505 | 876 | 45/28 |
