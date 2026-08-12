#!/usr/bin/env bash
set -euo pipefail
BASE="${SC_PCQTL_GIMAP_MIXING_ROOT:?Set SC_PCQTL_GIMAP_MIXING_ROOT.}"
VENV="${BASE}/env/tensorqtl_venv"
PY="${VENV}/bin/python"

if [[ ! -x "${PY}" ]]; then
  if command -v uv >/dev/null 2>&1; then
    uv venv "${VENV}" --python python3
  else
    echo "ERROR: uv is required because system python3 has no pip. Add uv to PATH or create ${VENV} manually." >&2
    exit 1
  fi
fi

if ! "${PY}" - <<'PY' >/dev/null 2>&1
import tensorqtl, pandas, numpy, torch, pandas_plink
PY
then
  uv pip install --python "${PY}" tensorqtl pandas numpy pyarrow scipy statsmodels matplotlib seaborn pandas-plink
fi

"${PY}" - <<'PY'
import tensorqtl, pandas, numpy, torch, pandas_plink
print('TensorQTL environment ready')
print('tensorqtl', getattr(tensorqtl, '__version__', 'unknown'))
print('pandas', pandas.__version__)
print('numpy', numpy.__version__)
print('torch', torch.__version__)
PY
