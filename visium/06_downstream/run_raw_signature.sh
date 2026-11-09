#!/usr/bin/env bash
set -euo pipefail
: "${PYTHON:?Set PYTHON to the pinned Scanpy 1.10.0rc2 environment Python}"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0
export OMP_NUM_THREADS=64 OPENBLAS_NUM_THREADS=64 MKL_NUM_THREADS=64
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec "$PYTHON" -s "$script_dir/03_raw_to_signature.py" "$@"
