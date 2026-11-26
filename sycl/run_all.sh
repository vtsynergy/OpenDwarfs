#!/bin/bash
# =============================================================================
# Run every SYCL benchmark once (or REPEATS times) and log results.
#
#   ./run_all.sh [--device cpu|gpu] [--size small|medium] [--repeats N] [--only app1,app2]
#                [--build] [--outdir DIR] [--cpuset LIST]
#
# Environment: GPU_INDEX (default 0, toolchains/env.sh), ONEAPI_ROOT (default $INTEL_ONEAPI_ROOT).
# Output: <outdir>/summary.csv (app,device,size,rep,status,wall_s) + per-run logs.
# OK means: exit status 0, a DWARF_RESULT line with valid=1, and the requested device.
# For comparison of computed outputs with the OpenCL tree use tools/validate.sh.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(dirname "$HERE")"
source "$REPO/toolchains/env.sh"
export ONEAPI_ROOT="${ONEAPI_ROOT:-$INTEL_ONEAPI_ROOT}"
source "$HERE/benchmarks.sh"
DEVICE=cpu; SIZE=small; REPEATS=1; ONLY=""; DO_BUILD=0; CPUSET=""
OUT="${HERE}/results/$(date +%Y%m%d-%H%M%S)"
while [ $# -gt 0 ]; do case $1 in
  --device) DEVICE=$2; shift 2;; --size) SIZE=$2; shift 2;; --repeats) REPEATS=$2; shift 2;;
  --only) ONLY=$2; shift 2;; --build) DO_BUILD=1; shift;; --outdir) OUT=$2; shift 2;; --cpuset) CPUSET=$2; shift 2;;
  *) echo "unknown arg $1"; exit 2;; esac; done
[ -n "$ONLY" ] && APPS="${ONLY//,/ }"
PIN=""; [ -n "$CPUSET" ] && PIN="taskset -c $CPUSET"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"; echo "app,device,size,rep,status,wall_s" > "$OUT/summary.csv"
fail=0
for a in $APPS; do
  d="${DIR[$a]:-}"; [ -n "$d" ] || { echo "unknown app $a"; continue; }
  if [ "$DO_BUILD" = 1 ] || [ ! -x "$HERE/$d/build/${BINNAME[$a]}" ]; then
    if ! make -C "$HERE/$d" build -j8 >"$OUT/$a.build.log" 2>&1; then echo "$a,$DEVICE,$SIZE,0,BLOCKED_BUILD," | tee -a "$OUT/summary.csv"; fail=1; continue; fi
  fi
  eval "args=\${ARGS_${SIZE}[$a]}"
  for r in $(seq 1 "$REPEATS"); do
    log="$OUT/$a.$DEVICE.$SIZE.$r.log"; t0=$(date +%s.%N)
    $PIN make -s -C "$HERE/$d" run SYCL_DEVICE=$DEVICE GPU_INDEX=$GPU_INDEX RUN_ARGS="$args" >"$log" 2>&1
    rc=$?; t1=$(date +%s.%N)
    line=$(grep -a '^DWARF_RESULT' "$log" | tail -1)
    if [ -z "$line" ]; then st=BLOCKED_RUNTIME
    elif ! grep -a -q "SYCL requested=$DEVICE " "$log" || ! grep -a -q ",device=$DEVICE," <<< "$line"; then st=WRONG_DEVICE
    elif [ $rc != 0 ] || ! grep -a -q ',valid=1' <<< "$line"; then st=FAIL:$rc
    else st=OK; fi
    [ $st = OK ] || fail=1
    echo "$a,$DEVICE,$SIZE,$r,$st,$(echo "$t1 - $t0" | bc -l | xargs printf %.3f)" | tee -a "$OUT/summary.csv"
  done
done
exit $fail
