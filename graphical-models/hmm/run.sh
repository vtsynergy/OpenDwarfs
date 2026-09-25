#!/usr/bin/env bash
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
device="${1:?usage: $0 <device> [make overrides...]}"
shift
"${script_dir}/../../common/run_benchmark.sh" "${script_dir}" "${device}" "$@"
