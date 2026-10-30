# Rebuilt node-local mokapot env (D38/D39)

The node wipe on 2026-10-29 deleted `/tmp/vulcan_m1e/envs/mokapot`, the mokapot env whose digest is frozen in `config/a5.yaml` (3264719f…).
- **Rebuild:** it was rebuilt on 2026-10-30 with the committed `scripts/baselines/setup_mokapot_env.sh`.
- **Versions:** mokapot 0.10.0, numpy 1.26.4, pandas 2.1.4, scikit-learn 1.5.2, scipy 1.13.1.
- **Digest:** b5ecb13ab10ff68e7671c14a5025f53d06aa4962a40be202cb42008d81191989, which differs from the frozen value.
- **What is recorded here:** the full manifest and package list, so any future rebuild can be compared file by file.
- **D1 comparator check: not done.** It was started and then stopped before any comparison (D39), so no claim is made about equivalence with the old env.
