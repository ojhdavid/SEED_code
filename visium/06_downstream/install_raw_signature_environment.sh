#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 1 ]]; then echo "Usage: BASE_PYTHON=python3.10 BOOTSTRAP_PYTHON=python3 install_raw_signature_environment.sh NEW_PREFIX" >&2; exit 2; fi
prefix="$1"
if [[ -e "$prefix" ]]; then echo "Refusing to modify an existing environment: $prefix" >&2; exit 2; fi
base="${BASE_PYTHON:-python3.10}"
bootstrap="${BOOTSTRAP_PYTHON:-$base}"
"$base" -c 'import sys; assert sys.version_info[:3] == (3,10,12), sys.version'
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONDONTWRITEBYTECODE=1
"$base" -m venv --without-pip "$prefix"
"$bootstrap" -m pip install --isolated --index-url https://pypi.org/simple --no-deps --no-compile --target "$prefix/lib/python3.10/site-packages" pip==24.0
export PYTHONNOUSERSITE=1
"$prefix/bin/python" -s -m pip install --isolated --index-url https://pypi.org/simple --no-compile -r "$script_dir/environment-raw-signature.lock.txt"
"$prefix/bin/python" -s -m pip check
