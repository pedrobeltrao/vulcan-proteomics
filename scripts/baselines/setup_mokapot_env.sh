#!/usr/bin/env bash
# mokapot 0.10.0 (latest on PyPI) is incompatible with the pe env (numpy 2.5 removed np.float_; pandas 3 removed
# to_numeric(errors="ignore")). Build a small pinned env on local disk (the project Lustre quota is nearly full).
set -euo pipefail
P=${VULCAN_MOKAPOT_ENV:-/tmp/vulcan_m1e/envs/mokapot}
[ -x "$P/bin/mokapot" ] && exit 0
. /opt/conda/etc/profile.d/conda.sh
CONDA_PKGS_DIRS=/tmp/vulcan_m1e/conda_pkgs mamba create -y -p "$P" -c conda-forge python=3.11 "numpy<2" "pandas<2.2" \
  "scikit-learn<1.6" "scipy<1.14" matplotlib-base numba lxml joblib pip
"$P/bin/pip" install "mokapot==0.10.0"
"$P/bin/python" -c "import mokapot, numpy, pandas, sklearn; print('mokapot', mokapot.__version__, 'numpy', numpy.__version__, 'pandas', pandas.__version__, 'sklearn', sklearn.__version__)"
