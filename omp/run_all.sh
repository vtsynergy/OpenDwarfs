#!/bin/bash
# =============================================================================
# Run every OpenMP benchmark and collect the DWARF_RESULT lines into a CSV.
#
#   ./run_all.sh [--device cpu|v100] [--size small|medium] [--repeats N] [--warmup N]
#                [--only app1,app2] [--build] [--outdir DIR]
#
# Uses the Clang installed by toolchains/install-omp-llvm-offload.sh
# ($OMP_PREFIX), CUDA from the NVHPC module, GPU_INDEX from toolchains/env.sh.
# GPU runs set OMP_TARGET_OFFLOAD=MANDATORY so a silent host fallback fails the run.
# Output: <outdir>/summary.csv + per-app logs.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../toolchains/env.sh"
DEVICE=cpu; SIZE=small; REPEATS=5; WARMUP=1; ONLY=""; DO_BUILD=0
OUT="$HERE/results/$(date +%Y%m%d-%H%M%S)"
while [ $# -gt 0 ]; do case $1 in
  --device) DEVICE=$2; shift 2;; --size) SIZE=$2; shift 2;; --repeats) REPEATS=$2; shift 2;;
  --warmup) WARMUP=$2; shift 2;; --only) ONLY=$2; shift 2;; --build) DO_BUILD=1; shift;;
  --outdir) OUT=$2; shift 2;; *) echo "unknown arg $1"; exit 2;; esac; done
[ "$DEVICE" = gpu ] && DEVICE=$GPU_NAME
BENCHES=$(grep -o '^[[:space:]][a-zA-Z/_-]*' "$HERE/common/benchmarks.mk" | tr -d ' \t' | tr '\n' ' ')
mkdir -p "$OUT"
export OMP_TARGET_OFFLOAD=$([ "$DEVICE" = cpu ] && echo DISABLED || echo MANDATORY)
[ "$DEVICE" != cpu ] && export CUDA_VISIBLE_DEVICES=$GPU_INDEX
T="$HERE/test"
# per-app extra make overrides for each size class (sizes are app-specific; see README)
declare -A SMALL MEDIUM
set_sz() { SMALL[$1]="$2"; MEDIUM[$1]="$3"; }
set_sz dense-linear-algebra/gemm "RUN_SIZE=512"  "RUN_SIZE=2048"
set_sz dense-linear-algebra/atax "RUN_SIZE=2048" "RUN_SIZE=8192"
set_sz dense-linear-algebra/cg   "RUN_SIZE=1024" "RUN_SIZE=4096"
set_sz dense-linear-algebra/lud  "RUN_INPUT=$T/dense-linear-algebra/lud/256.dat" "RUN_INPUT=$T/dense-linear-algebra/lud/2048.dat"
set_sz dense-linear-algebra/kmeans "RUN_INPUT=$T/dense-linear-algebra/kmeans/204800.txt" "RUN_INPUT=$T/dense-linear-algebra/kmeans/819200.txt"
set_sz branch-and-bound/nqueens "RUN_SIZE=12" "RUN_SIZE=14"
set_sz dynamic-programming/nw   "RUN_SIZE=1024" "RUN_SIZE=4096"
set_sz dynamic-programming/swat "RUN_INPUT=$T/dynamic-programming/swat/query1K1 RUN_INPUT2=$T/dynamic-programming/swat/sampledb1K1" "RUN_INPUT=$T/dynamic-programming/swat/query3K1 RUN_INPUT2=$T/dynamic-programming/swat/sampledb3K1"
set_sz finite-state-machine/tdm "RUN_INPUT=$T/finite-state-machine/tdm/sim-64-size50.csv RUN_INPUT2=$T/finite-state-machine/tdm/ivl.txt RUN_INPUT3=$T/finite-state-machine/tdm/30-episodes.txt" "RUN_INPUT=$T/finite-state-machine/tdm/sim-64-size500.csv RUN_INPUT2=$T/finite-state-machine/tdm/ivl.txt RUN_INPUT3=$T/finite-state-machine/tdm/30-episodes.txt"
set_sz graphical-models/hmm     "RUN_SIZE=60" "RUN_SIZE=120"
set_sz structured-grids/fdtd2d "RUN_SIZE=1024" "RUN_SIZE=4096"
set_sz structured-grids/srad    "RUN_SIZE=512" "RUN_SIZE=2048"
set_sz sparse-linear-algebra/spmv "RUN_SIZE=2048" "RUN_SIZE=8192"
set_sz spectral-methods/fft     "RUN_SIZE=65536" "RUN_SIZE=4194304"
set_sz graph-traversal/bfs      "RUN_INPUT=$T/graph-traversal/bfs/graph4096.txt" "RUN_INPUT=$T/graph-traversal/bfs/graph65536.txt"
set_sz n-body-methods/gem       "RUN_INPUT=$T/n-body-methods/gem/1uwo_A.pqr RUN_INPUT2=$T/n-body-methods/gem/1uwo_A.vert" "RUN_INPUT=$T/n-body-methods/gem/Mb.Hhelix.bondi.pqr RUN_INPUT2=$T/n-body-methods/gem/Mb.Hhelix.bondi.vert"
set_sz unstructured-grids/cfd   "RUN_INPUT=$T/unstructured-grids/cfd/fvcorr.domn.097K" "RUN_INPUT=$T/unstructured-grids/cfd/fvcorr.domn.193K"
set_sz mapreduce/pagerank       "RUN_SIZE=4096" "RUN_SIZE=65536"
echo "benchmark,device,size_class,rc,valid,mean_ms,median_ms,min_ms,max_ms,overrides" > "$OUT/summary.csv"
for b in $BENCHES; do
  short=$(basename "$b")
  if [ -n "$ONLY" ] && ! [[ ",$ONLY," == *",$short,"* ]]; then continue; fi
  ov="${SMALL[$b]:-}"; [ "$SIZE" = medium ] && ov="${MEDIUM[$b]:-}"
  log="$OUT/$short.$DEVICE.$SIZE.log"
  # shellcheck disable=SC2086
  make -C "$HERE/$b" run DEVICE=$DEVICE RUN_REPEATS=$REPEATS RUN_WARMUP=$WARMUP CUDA_PATH="${CUDA_PATH:-}" $ov >"$log" 2>&1
  rc=$?
  line=$(grep '^DWARF_RESULT,' "$log" | tail -1)
  f() { echo "$line" | tr ',' '\n' | sed -n "s/^$1=//p"; }
  echo "$short,$DEVICE,$SIZE,$rc,$(f valid),$(f mean_ms),$(f median_ms),$(f min_ms),$(f max_ms),\"$ov\"" | tee -a "$OUT/summary.csv"
done
