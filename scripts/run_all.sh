#!/bin/bash
# =============================================================================
# One entry point for all three implementations.
#   scripts/run_all.sh <opencl|omp|openacc|all> [args passed to the model's run_all.sh]
# Examples:
#   scripts/run_all.sh opencl --vendor nvidia --size small
#   scripts/run_all.sh omp --device v100 --size medium --repeats 5
#   scripts/run_all.sh openacc --run
# The model scripts live in <model>/run_all.sh; GPU/arch/prefix settings are in
# toolchains/env.sh (GPU_INDEX=... to pick another GPU).
# =============================================================================
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
m=${1:-all}; shift || true
case $m in
  opencl|omp|openacc) exec "$ROOT/$m/run_all.sh" "$@";;
  all) rc=0; for x in opencl omp openacc; do echo "##### $x"; "$ROOT/$x/run_all.sh" "$@" || rc=1; done; exit $rc;;
  *) echo "usage: $0 <opencl|omp|openacc|all> [args]"; exit 2;;
esac
