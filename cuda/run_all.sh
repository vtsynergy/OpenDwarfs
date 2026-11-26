#!/bin/bash
# =============================================================================
# Build and run every CUDA benchmark once (or REPEATS times) and log results.  Counterpart of ../opencl/run_all.sh.
#
#   ./run_all.sh [--size small|medium] [--repeats N] [--only app1,app2] [--outdir DIR]
#
# Environment: GPU_INDEX (default 0, see toolchains/env.sh), CUDACC (nvcc-compatible driver; default NVIDIA nvcc),
#              BUILD_DIR (per-app build directory name, default build; use another one per compiler).
# Output: <outdir>/summary.csv (app,size,rep,status) + per-run logs.  status:
#   PASS            exit 0, DWARF_RESULT valid=1 (the application's own validator)
#   FAIL            the program ran but validation failed or it exited non-zero after reporting
#   BLOCKED_BUILD   make build failed
#   BLOCKED_RUNTIME no DWARF_RESULT line (crash, launch error, timeout)
#   WRONG_DEVICE    the chosen device is not an NVIDIA GPU
# These are validator results only; comparison of outputs with OpenCL is described in README.md.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../toolchains/env.sh"
SIZE=small; REPEATS=1; ONLY=""
OUT="${HERE}/results/$(date +%Y%m%d-%H%M%S)"
BUILD_DIR="${BUILD_DIR:-build}"
while [ $# -gt 0 ]; do case $1 in
  --size) SIZE=$2; shift 2;; --repeats) REPEATS=$2; shift 2;; --only) ONLY=$2; shift 2;; --outdir) OUT=$2; shift 2;;
  *) echo "unknown arg $1"; exit 2;; esac; done
T="$HERE/../opencl/test"          # input files are shared with the OpenCL tree (read only)
mkdir -p "$OUT/data"; OUT="$(cd "$OUT" && pwd)"
echo "app,size,rep,status" > "$OUT/summary.csv"
MK=(BUILD_DIR="$BUILD_DIR" GPU_INDEX="$GPU_INDEX"); [ -n "${CUDACC:-}" ] && MK+=(CUDACC="$CUDACC")

declare -A DIR ARGS_small ARGS_medium
reg() { DIR[$1]=$2; ARGS_small[$1]="$3"; ARGS_medium[$1]="$4"; }
reg nqueens branch-and-bound/nqueens        "12"                                   "14"
reg lud     dense-linear-algebra/lud        "-i $T/dense-linear-algebra/lud/256.dat -v"  "-i $T/dense-linear-algebra/lud/2048.dat -v"
reg kmeans  dense-linear-algebra/kmeans     "-i $T/dense-linear-algebra/kmeans/204800.txt" "-i $T/dense-linear-algebra/kmeans/819200.txt"
reg nw      dynamic-programming/nw          "1024 10"                              "4096 10"
reg swat    dynamic-programming/swat        "$T/dynamic-programming/swat/query1K1 $T/dynamic-programming/swat/sampledb1K1" "$T/dynamic-programming/swat/query3K1 $T/dynamic-programming/swat/sampledb3K1"
reg tdm     finite-state-machine/tdm        "$T/finite-state-machine/tdm/sim-64-size50.csv $T/finite-state-machine/tdm/ivl.txt $T/finite-state-machine/tdm/30-episodes.txt 128" "$T/finite-state-machine/tdm/sim-64-size500.csv $T/finite-state-machine/tdm/ivl.txt $T/finite-state-machine/tdm/30-episodes.txt 128"
reg bfs     graph-traversal/bfs             "$T/graph-traversal/bfs/graph4096.txt" "$T/graph-traversal/bfs/graph65536.txt"
reg hmm     graphical-models/hmm            "-n 60 -v n"                           "-n 120 -v n"
reg gem     n-body-methods/gem              "$OUT/data/1uwo_A 80 1 0"              "$OUT/data/Mb.Hhelix.bondi 80 1 0"
reg srad    structured-grids/srad           "512 512 0 127 0 127 0.5 2"            "2048 2048 0 127 0 127 0.5 2"
reg cfd     unstructured-grids/cfd          "$T/unstructured-grids/cfd/fvcorr.domn.097K" "$T/unstructured-grids/cfd/fvcorr.domn.193K"
reg crc     combinational-logic/crc         "-i $OUT/data/crc_small.dat -a"        "-i $OUT/data/crc_medium.dat -a"
reg spmv    sparse-linear-algebra/SPMV      "-i $OUT/data/csr_small.dat -a"        "-i $OUT/data/csr_medium.dat -a"
reg fft     spectral-methods/fft            "--pts 1 --passes 5"                    "--pts 3 --passes 5"

APPS="nqueens crc lud kmeans nw swat tdm bfs hmm gem spmv srad cfd fft"
[ -n "$ONLY" ] && APPS="${ONLY//,/ }"
cp -f "$T"/n-body-methods/gem/1uwo_A.* "$T"/n-body-methods/gem/Mb.Hhelix.bondi.* "$OUT/data/"   # gem gets private copies
for a in $APPS; do
  d="${DIR[$a]:-}"; [ -n "$d" ] || { echo "unknown app $a"; continue; }
  if ! make -C "$HERE/$d" build -j8 "${MK[@]}" >"$OUT/$a.build.log" 2>&1; then
    echo "$a,$SIZE,0,BLOCKED_BUILD" | tee -a "$OUT/summary.csv"; continue
  fi
  B="$HERE/$d/$BUILD_DIR"
  [ $a = crc ]  && [ ! -f "$OUT/data/crc_$SIZE.dat" ] && "$B/createcrc" -R -n $([ $SIZE = small ] && echo 1024 || echo 16384) -s 1024 -f "$OUT/data/crc_$SIZE.dat" >/dev/null
  [ $a = spmv ] && [ ! -f "$OUT/data/csr_$SIZE.dat" ] && "$B/createcsr" -R -n $([ $SIZE = small ] && echo 2048 || echo 8192) -d 5000 -f "$OUT/data/csr_$SIZE.dat" >/dev/null
  eval "args=\${ARGS_${SIZE}[$a]}"
  for r in $(seq 1 "$REPEATS"); do
    log="$OUT/$a.$SIZE.$r.log"
    timeout 900 make -C "$HERE/$d" run "${MK[@]}" RUN_ARGS="$args" >"$log" 2>&1; rc=$?
    line=$(grep -a '^DWARF_RESULT' "$log" | tail -1)
    if ! grep -a -m1 'Device Chosen' "$log" | grep -qi nvidia; then st=$([ -n "$line" ] && echo WRONG_DEVICE || echo BLOCKED_RUNTIME)
    elif [ -z "$line" ]; then st=BLOCKED_RUNTIME
    elif [ $rc = 0 ] && echo "$line" | grep -q 'valid=1'; then st=PASS
    else st=FAIL; fi
    echo "$a,$SIZE,$r,$st" | tee -a "$OUT/summary.csv"
  done
done
echo "summary: $OUT/summary.csv"
! grep -qE ',(FAIL|BLOCKED_BUILD|BLOCKED_RUNTIME|WRONG_DEVICE)$' "$OUT/summary.csv"
