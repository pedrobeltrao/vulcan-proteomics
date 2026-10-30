# T1 driver: D0 (dev), runs=2, decoy seeds [1001, 1002, 1003, 1004, 1005], entrapment seeds [2001, 2002, 2003]

## IDs at own q (production searches; mean over the 5 decoy shuffles, dataset level)

| method | level | 0.1% | 1% | 5% | verdict 1% (a / b) | verdict 5% (a / b) |
|---|---|---|---|---|---|---|
| vulcan_v2_sep | peptidoform | 25793 | 29484 ± 28 | 31651 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan_v2_sep | sequence | 24455 | 27746 ± 27 | 29798 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan_v2_sep | protein | 6571 | 6979 ± 8 | 6982 | inconclusive / inconclusive | supported / supported |
| sage+pred+mokapot | peptidoform | 23530 | 24576 ± 31 | 25656 | inconclusive / inconclusive | supported / supported |
| sage+pred+mokapot | sequence | 22207 | 23180 ± 33 | 24243 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage+pred+mokapot | protein | 6255 | 6465 ± 15 | 6468 | inconclusive / inconclusive | supported / supported |
| sage_chim2 | peptidoform | 21094 | 23944 ± 43 | 25795 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage_chim2 | sequence | 19999 | 22644 ± 42 | 24420 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage_chim2 | protein | 5884 | 6339 ± 13 | 6342 | inconclusive / inconclusive | supported / supported |
| sage_chim2+pred+mokapot | peptidoform | 24010 | 25293 ± 23 | 26514 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage_chim2+pred+mokapot | sequence | 22656 | 23852 ± 22 | 25046 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage_chim2+pred+mokapot | protein | 6325 | 6553 ± 7 | 6555 | inconclusive / inconclusive | supported / supported |
| sage_chim5 | peptidoform | 20702 | 23822 ± 50 | 25753 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage_chim5 | sequence | 19630 | 22534 ± 48 | 24390 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage_chim5 | protein | 5830 | 6321 ± 18 | 6328 | inconclusive / inconclusive | supported / supported |
| sage_chim5+pred+mokapot | peptidoform | 24084 | 25458 ± 27 | 26795 | inconclusive / inconclusive | supported / supported |
| sage_chim5+pred+mokapot | sequence | 22725 | 24008 ± 28 | 25311 | inconclusive / inconclusive | inconclusive / inconclusive |
| sage_chim5+pred+mokapot | protein | 6329 | 6571 ± 12 | 6574 | inconclusive / inconclusive | supported / supported |

## Entrapment FDP estimates at 1% (mean [95% CI, convention a]; arm P = 15 replicates, arm C = 5)

| method | level | P mean | P CI (a) | P CI (b) | C mean | C CI (a) | LB(P) |
|---|---|---|---|---|---|---|---|
| vulcan_v2_sep | peptidoform | 0.0107 | [0.0101, 0.0114] | [0.0103, 0.0112] | 0.0077 | [0.0070, 0.0084] | 0.0055 |
| vulcan_v2_sep | sequence | 0.0113 | [0.0106, 0.0120] | [0.0108, 0.0118] | 0.0085 | [0.0077, 0.0092] | 0.0058 |
| vulcan_v2_sep | protein | 0.0106 | [0.0093, 0.0120] | [0.0097, 0.0115] | - | - | 0.0068 |
| sage+pred+mokapot | peptidoform | 0.0105 | [0.0098, 0.0111] | [0.0100, 0.0109] | 0.0093 | [0.0085, 0.0101] | 0.0053 |
| sage+pred+mokapot | sequence | 0.0110 | [0.0103, 0.0117] | [0.0105, 0.0115] | 0.0102 | [0.0093, 0.0111] | 0.0056 |
| sage+pred+mokapot | protein | 0.0113 | [0.0098, 0.0127] | [0.0102, 0.0123] | - | - | 0.0071 |
| sage_chim2 | peptidoform | 0.0109 | [0.0102, 0.0115] | [0.0104, 0.0113] | 0.0087 | [0.0079, 0.0096] | 0.0055 |
| sage_chim2 | sequence | 0.0114 | [0.0107, 0.0120] | [0.0109, 0.0119] | 0.0095 | [0.0086, 0.0104] | 0.0058 |
| sage_chim2 | protein | 0.0121 | [0.0106, 0.0136] | [0.0111, 0.0130] | - | - | 0.0077 |
| sage_chim2+pred+mokapot | peptidoform | 0.0102 | [0.0096, 0.0109] | [0.0098, 0.0107] | 0.0086 | [0.0078, 0.0093] | 0.0052 |
| sage_chim2+pred+mokapot | sequence | 0.0108 | [0.0101, 0.0114] | [0.0103, 0.0112] | 0.0094 | [0.0086, 0.0102] | 0.0055 |
| sage_chim2+pred+mokapot | protein | 0.0113 | [0.0100, 0.0126] | [0.0103, 0.0122] | - | - | 0.0073 |
| sage_chim5 | peptidoform | 0.0104 | [0.0098, 0.0111] | [0.0099, 0.0109] | 0.0086 | [0.0078, 0.0094] | 0.0053 |
| sage_chim5 | sequence | 0.0109 | [0.0103, 0.0116] | [0.0104, 0.0114] | 0.0094 | [0.0085, 0.0102] | 0.0056 |
| sage_chim5 | protein | 0.0116 | [0.0102, 0.0129] | [0.0106, 0.0125] | - | - | 0.0075 |
| sage_chim5+pred+mokapot | peptidoform | 0.0107 | [0.0100, 0.0113] | [0.0102, 0.0112] | 0.0084 | [0.0077, 0.0091] | 0.0055 |
| sage_chim5+pred+mokapot | sequence | 0.0112 | [0.0106, 0.0119] | [0.0107, 0.0117] | 0.0092 | [0.0084, 0.0100] | 0.0058 |
| sage_chim5+pred+mokapot | protein | 0.0108 | [0.0096, 0.0119] | [0.0099, 0.0117] | - | - | 0.0072 |

## Non-inferiority (ours / comparator; run-level IDs summed over runs, mean over decoy shuffles; 95% bootstrap CI over runs x decoy shuffles; pass = lower CI >= 0.97)

| comparator | level | own q 1%: ratio [CI] | counts? (comparator verdict) | matched FDP 1%: ratio [CI] | matched FDP 5%: ratio [CI] |
|---|---|---|---|---|---|
| sage+pred+mokapot | peptidoform | 1.212 [1.205, 1.220] PASS | no (inconclusive) | 1.211 [1.208, 1.214] PASS | 1.245 [1.240, 1.247] PASS |
| sage+pred+mokapot | sequence | 1.213 [1.206, 1.220] PASS | no (inconclusive) | 1.211 [1.208, 1.214] PASS | 1.242 [1.236, 1.246] PASS |
| sage+pred+mokapot | protein | 1.144 [1.129, 1.160] PASS | no (inconclusive) | 1.144 [1.129, 1.161] PASS | 1.144 [1.129, 1.161] PASS |
| sage_chim2 | peptidoform | 1.232 [1.228, 1.237] PASS | no (inconclusive) | 1.235 [1.230, 1.240] PASS | 1.232 [1.227, 1.236] PASS |
| sage_chim2 | sequence | 1.227 [1.222, 1.232] PASS | no (inconclusive) | 1.230 [1.225, 1.236] PASS | 1.227 [1.220, 1.231] PASS |
| sage_chim2 | protein | 1.090 [1.088, 1.092] PASS | no (inconclusive) | 1.094 [1.090, 1.097] PASS | 1.090 [1.088, 1.091] PASS |
| sage_chim2+pred+mokapot | peptidoform | 1.168 [1.162, 1.175] PASS | no (inconclusive) | 1.166 [1.164, 1.171] PASS | 1.195 [1.190, 1.198] PASS |
| sage_chim2+pred+mokapot | sequence | 1.167 [1.162, 1.173] PASS | no (inconclusive) | 1.165 [1.163, 1.168] PASS | 1.191 [1.186, 1.195] PASS |
| sage_chim2+pred+mokapot | protein | 1.060 [1.057, 1.062] PASS | no (inconclusive) | 1.060 [1.058, 1.064] PASS | 1.059 [1.057, 1.063] PASS |
| sage_chim5 | peptidoform | 1.227 [1.219, 1.234] PASS | no (inconclusive) | 1.228 [1.218, 1.237] PASS | 1.223 [1.213, 1.230] PASS |
| sage_chim5 | sequence | 1.220 [1.211, 1.228] PASS | no (inconclusive) | 1.221 [1.210, 1.231] PASS | 1.215 [1.203, 1.223] PASS |
| sage_chim5 | protein | 1.051 [1.044, 1.059] PASS | no (inconclusive) | 1.057 [1.049, 1.064] PASS | 1.050 [1.042, 1.058] PASS |
| sage_chim5+pred+mokapot | peptidoform | 1.150 [1.147, 1.154] PASS | no (inconclusive) | 1.149 [1.147, 1.152] PASS | 1.171 [1.165, 1.175] PASS |
| sage_chim5+pred+mokapot | sequence | 1.147 [1.145, 1.150] PASS | no (inconclusive) | 1.146 [1.143, 1.150] PASS | 1.166 [1.158, 1.170] PASS |
| sage_chim5+pred+mokapot | protein | 1.015 [1.010, 1.021] PASS | no (inconclusive) | 1.016 [1.009, 1.021] PASS | 1.015 [1.010, 1.021] PASS |

Sage predict_rt check (prod, first decoy seed; Sage fits its RT model only when on): {'predict_rt_true': {'n_confident': 43967, 'frac_predicted_rt_nonzero': 0.9823871268400998, 'sd_predicted_rt': 0.2359752486561844, 'median_abs_delta_rt_model_confident': 0.05510366, 'spearman_aligned_vs_predicted_confident': 0.9210243368620381, 'pin_has_rt_features': True}, 'predict_rt_false': {'n_confident': 43605, 'frac_predicted_rt_nonzero': 0.0, 'sd_predicted_rt': 0.0, 'median_abs_delta_rt_model_confident': 0.999, 'spearman_aligned_vs_predicted_confident': nan, 'pin_has_rt_features': True}}

Sage RT on vs off (prod, first decoy seed): {'peptidoform': {'0.001': {'rt_on': 20865, 'rt_off': 20621}, '0.01': {'rt_on': 23732, 'rt_off': 23432}}, 'sequence': {'0.001': {'rt_on': 19778, 'rt_off': 19530}, '0.01': {'rt_on': 22431, 'rt_off': 22142}}}

## Driver timing: total 6942 s wall, 7747 core-s

| stage | wall s | core-s (self+children) | loadavg start/end |
|---|---|---|---|
| db | 1 | 0 | 12/12 |
| evaluate | 6939 | 7747 | 12/11 |
