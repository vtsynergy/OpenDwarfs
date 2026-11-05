#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "usage: $0 <benchmark-dir> <device> [make overrides...]" >&2
  echo "  e.g.: $0 . a100 RUN_SIZE=4096" >&2
  exit 2
fi

bench_dir="$1"
device="$2"
shift 2

make -C "${bench_dir}" run DEVICE="${device}" "$@"
