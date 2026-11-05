#!/bin/bash
# =============================================================================
# Build (and run) the OpenACC benchmarks with the NVIDIA HPC SDK.
#   ./run_all.sh [--target gpu|multicore] [--run] [--outdir DIR]
# By default only builds and reports BUILD OK/FAIL per benchmark; --run executes each
# successfully built binary on GPU_INDEX. The tree is the 2014 PGI/CoreTSAR code,
# kept unmodified; benchmarks using obsolete syntax are reported as BUILD FAIL.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../toolchains/env.sh"
source "$ACC_PREFIX/env.sh" || { echo "run toolchains/install-openacc-nvhpc-module.sh first"; exit 1; }
TARGET=gpu; DO_RUN=0; OUT="$HERE/results/$(date +%Y%m%d-%H%M%S)"
while [ $# -gt 0 ]; do case $1 in
  --target) TARGET=$2; shift 2;; --run) DO_RUN=1; shift;; --outdir) OUT=$2; shift 2;;
  *) echo "unknown arg $1"; exit 2;; esac; done
mkdir -p "$OUT"
make -C "$HERE" clean >/dev/null 2>&1
make -C "$HERE" libctsar ACC_TARGET=$TARGET GPU_CC=$ACC_GPU_CC > "$OUT/libctsar.log" 2>&1 || { echo "libctsar failed, see $OUT/libctsar.log"; exit 1; }
echo "benchmark,target,build,run_rc,wall_s" > "$OUT/summary.csv"
for d in $(grep -o '^[[:space:]][a-zA-Z/_-]*' "$HERE/benchmarks.mk" | tr -d ' \t' | tr '\n' ' '); do
  n=$(echo "$d" | sed 's|/|_|g')
  if make -C "$HERE/$d" build ACC_TARGET=$TARGET GPU_CC=$ACC_GPU_CC > "$OUT/$n.build.log" 2>&1; then b=OK; else b=FAIL; fi
  rc=NA; wall=NA
  if [ $b = OK ] && [ $DO_RUN = 1 ]; then
    t0=$(date +%s.%N)
    make -C "$HERE/$d" run ACC_TARGET=$TARGET GPU_INDEX=$GPU_INDEX > "$OUT/$n.run.log" 2>&1; rc=$?
    wall=$(echo "$(date +%s.%N) - $t0" | bc -l | xargs printf %.3f)
  fi
  echo "$d,$TARGET,$b,$rc,$wall" | tee -a "$OUT/summary.csv"
done
