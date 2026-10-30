# T1 driver: D0 (dev), runs=2, decoy seeds [1001, 1002, 1003, 1004, 1005], entrapment seeds [2001, 2002, 2003]

## IDs at own q (production searches; mean over the 5 decoy shuffles, dataset level)

| method | level | 0.1% | 1% | 5% | verdict 1% (a / b) | verdict 5% (a / b) |
|---|---|---|---|---|---|---|
| vulcan | peptidoform | 23262 | 24328 ± 36 | 25486 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan | sequence | 22014 | 22981 ± 34 | 24107 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan | protein | 6216 | 6429 ± 15 | 6430 | inconclusive / inconclusive | supported / supported |
| vulcan_chim_v1 | peptidoform | 27273 | 29165 ± 35 | 30597 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan_chim_v1 | sequence | 25729 | 27411 ± 34 | 28812 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan_chim_v1 | protein | 6664 | 6951 ± 13 | 6955 | inconclusive / inconclusive | supported / supported |
| vulcan_chim_v1_sep | peptidoform | 26596 | 29056 ± 47 | 30745 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan_chim_v1_sep | sequence | 25121 | 27311 ± 44 | 28942 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan_chim_v1_sep | protein | 6664 | 6951 ± 13 | 6955 | inconclusive / inconclusive | supported / supported |
| vulcan_v2 | peptidoform | 26539 | 29640 ± 56 | 31642 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan_v2 | sequence | 25148 | 27884 ± 43 | 29790 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan_v2 | protein | 6571 | 6979 ± 8 | 6982 | inconclusive / inconclusive | supported / supported |
| vulcan_v2_sep | peptidoform | 25793 | 29484 ± 28 | 31651 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan_v2_sep | sequence | 24455 | 27746 ± 27 | 29798 | inconclusive / inconclusive | inconclusive / inconclusive |
| vulcan_v2_sep | protein | 6571 | 6979 ± 8 | 6982 | inconclusive / inconclusive | supported / supported |
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
| vulcan_chim_v1 | peptidoform | 0.0105 | [0.0099, 0.0111] | [0.0101, 0.0109] | 0.0087 | [0.0078, 0.0095] | 0.0054 |
| vulcan_chim_v1 | sequence | 0.0110 | [0.0104, 0.0117] | [0.0106, 0.0115] | 0.0095 | [0.0085, 0.0104] | 0.0057 |
| vulcan_chim_v1 | protein | 0.0129 | [0.0117, 0.0140] | [0.0119, 0.0139] | - | - | 0.0084 |
| vulcan_chim_v1_sep | peptidoform | 0.0100 | [0.0095, 0.0106] | [0.0096, 0.0105] | 0.0088 | [0.0079, 0.0096] | 0.0051 |
| vulcan_chim_v1_sep | sequence | 0.0106 | [0.0100, 0.0111] | [0.0101, 0.0110] | 0.0096 | [0.0087, 0.0105] | 0.0054 |
| vulcan_chim_v1_sep | protein | 0.0129 | [0.0117, 0.0140] | [0.0119, 0.0139] | - | - | 0.0084 |
| vulcan_v2 | peptidoform | 0.0105 | [0.0099, 0.0111] | [0.0101, 0.0110] | 0.0082 | [0.0074, 0.0090] | 0.0054 |
| vulcan_v2 | sequence | 0.0111 | [0.0105, 0.0117] | [0.0106, 0.0115] | 0.0090 | [0.0081, 0.0098] | 0.0057 |
| vulcan_v2 | protein | 0.0106 | [0.0093, 0.0120] | [0.0097, 0.0115] | - | - | 0.0068 |
| vulcan_v2_sep | peptidoform | 0.0107 | [0.0101, 0.0114] | [0.0103, 0.0112] | 0.0077 | [0.0070, 0.0084] | 0.0055 |
| vulcan_v2_sep | sequence | 0.0113 | [0.0106, 0.0120] | [0.0108, 0.0118] | 0.0085 | [0.0077, 0.0092] | 0.0058 |
| vulcan_v2_sep | protein | 0.0106 | [0.0093, 0.0120] | [0.0097, 0.0115] | - | - | 0.0068 |
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
| vulcan | peptidoform | 1.225 [1.220, 1.230] PASS | no (inconclusive) | 1.222 [1.220, 1.226] PASS | 1.255 [1.250, 1.258] PASS |
| vulcan | sequence | 1.223 [1.219, 1.228] PASS | no (inconclusive) | 1.221 [1.219, 1.225] PASS | 1.251 [1.245, 1.254] PASS |
| vulcan | protein | 1.149 [1.135, 1.166] PASS | no (inconclusive) | 1.149 [1.134, 1.167] PASS | 1.149 [1.134, 1.166] PASS |
| vulcan_chim_v1 | peptidoform | 1.013 [1.010, 1.016] PASS | no (inconclusive) | 1.012 [1.011, 1.015] PASS | 1.036 [1.030, 1.041] PASS |
| vulcan_chim_v1 | sequence | 1.014 [1.011, 1.017] PASS | no (inconclusive) | 1.014 [1.012, 1.016] PASS | 1.035 [1.028, 1.041] PASS |
| vulcan_chim_v1 | protein | 1.007 [1.004, 1.009] PASS | no (inconclusive) | 1.008 [1.005, 1.011] PASS | 1.006 [1.003, 1.009] PASS |
| vulcan_chim_v1_sep | peptidoform | 1.016 [1.013, 1.020] PASS | no (inconclusive) | 1.015 [1.013, 1.017] PASS | 1.032 [1.028, 1.035] PASS |
| vulcan_chim_v1_sep | sequence | 1.018 [1.015, 1.021] PASS | no (inconclusive) | 1.017 [1.014, 1.020] PASS | 1.032 [1.027, 1.035] PASS |
| vulcan_chim_v1_sep | protein | 1.007 [1.004, 1.009] PASS | no (inconclusive) | 1.008 [1.005, 1.011] PASS | 1.006 [1.003, 1.009] PASS |
| vulcan_v2 | peptidoform | 0.995 [0.993, 0.996] PASS | no (inconclusive) | 0.995 [0.992, 0.997] PASS | 1.000 [0.996, 1.003] PASS |
| vulcan_v2 | sequence | 0.995 [0.994, 0.996] PASS | no (inconclusive) | 0.995 [0.993, 0.997] PASS | 1.000 [0.996, 1.003] PASS |
| vulcan_v2 | protein | 1.000 [1.000, 1.000] PASS | no (inconclusive) | 1.000 [1.000, 1.000] PASS | 1.000 [1.000, 1.000] PASS |
| sage | peptidoform | 1.255 [1.252, 1.257] PASS | no (inconclusive) | 1.257 [1.253, 1.261] PASS | 1.262 [1.259, 1.264] PASS |
| sage | sequence | 1.252 [1.249, 1.255] PASS | no (inconclusive) | 1.255 [1.250, 1.259] PASS | 1.259 [1.256, 1.260] PASS |
| sage | protein | 1.165 [1.150, 1.182] PASS | no (inconclusive) | 1.166 [1.150, 1.183] PASS | 1.165 [1.150, 1.181] PASS |
| sage+mokapot | peptidoform | 1.250 [1.248, 1.253] PASS | no (inconclusive) | 1.254 [1.249, 1.257] PASS | 1.256 [1.254, 1.260] PASS |
| sage+mokapot | sequence | 1.248 [1.246, 1.250] PASS | no (inconclusive) | 1.253 [1.247, 1.255] PASS | 1.254 [1.251, 1.257] PASS |
| sage+mokapot | protein | 1.162 [1.148, 1.179] PASS | no (inconclusive) | 1.162 [1.147, 1.183] PASS | 1.162 [1.147, 1.179] PASS |
| sage+pred+mokapot | peptidoform | 1.212 [1.205, 1.220] PASS | no (inconclusive) | 1.211 [1.208, 1.214] PASS | 1.245 [1.240, 1.247] PASS |
| sage+pred+mokapot | sequence | 1.213 [1.206, 1.220] PASS | no (inconclusive) | 1.211 [1.208, 1.214] PASS | 1.242 [1.236, 1.246] PASS |
| sage+pred+mokapot | protein | 1.144 [1.129, 1.160] PASS | no (inconclusive) | 1.144 [1.129, 1.161] PASS | 1.144 [1.129, 1.161] PASS |
| comet+mokapot | peptidoform | 1.269 [1.260, 1.279] PASS | yes (supported) | 1.264 [1.257, 1.271] PASS | 1.265 [1.254, 1.274] PASS |
| comet+mokapot | sequence | 1.266 [1.256, 1.276] PASS | no (inconclusive) | 1.261 [1.255, 1.268] PASS | 1.262 [1.248, 1.271] PASS |
| comet+mokapot | protein | 1.146 [1.135, 1.160] PASS | yes (supported) | 1.146 [1.133, 1.159] PASS | 1.146 [1.133, 1.159] PASS |

Sage predict_rt check (prod, first decoy seed; Sage fits its RT model only when on): {'predict_rt_true': {'n_confident': 43967, 'frac_predicted_rt_nonzero': 0.9823871268400998, 'sd_predicted_rt': 0.2359752486561844, 'median_abs_delta_rt_model_confident': 0.05510366, 'spearman_aligned_vs_predicted_confident': 0.9210243368620381, 'pin_has_rt_features': True}, 'predict_rt_false': {'n_confident': 43605, 'frac_predicted_rt_nonzero': 0.0, 'sd_predicted_rt': 0.0, 'median_abs_delta_rt_model_confident': 0.999, 'spearman_aligned_vs_predicted_confident': nan, 'pin_has_rt_features': True}}

Sage RT on vs off (prod, first decoy seed): {'peptidoform': {'0.001': {'rt_on': 20865, 'rt_off': 20621}, '0.01': {'rt_on': 23732, 'rt_off': 23432}}, 'sequence': {'0.001': {'rt_on': 19778, 'rt_off': 19530}, '0.01': {'rt_on': 22431, 'rt_off': 22142}}}

## Driver timing: total 4086 s wall, 5292 core-s

| stage | wall s | core-s (self+children) | loadavg start/end |
|---|---|---|---|
| db | 0 | 0 | 16/16 |
| evaluate | 4085 | 5291 | 16/9 |
