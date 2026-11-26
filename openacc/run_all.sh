#!/bin/bash
# =============================================================================
# Build and run the 14 OpenACC translations of the OpenCL dwarfs once each, with
# each application's own validator. No timing claims are derived from this script.
#
#   ./run_all.sh [--target gpu|multicore|both] [--size small|medium] [--only a,b] [--outdir DIR]
#
# Environment: GPU_INDEX (default from toolchains/env.sh). For multicore runs set ACC_NUM_CORES and,
# if wanted, pin with taskset around this script.
# Output: <outdir>/summary.csv (app,target,size,status,rc) + logs.
#   PASS            DWARF_RESULT valid=1, exit 0, no FAILURE line
#   FAIL            the program ran and reported invalid / non-zero exit
#   BLOCKED_BUILD   build failed        BLOCKED_RUNTIME  no DWARF_RESULT line
# The legacy 2014 PGI/CoreTSAR code is not run here (see README.md: make legacy-build).
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../toolchains/env.sh"
TARGET=both; SIZE=small; ONLY=""; OUT="$HERE/results/$(date +%Y%m%d-%H%M%S)"
while [ $# -gt 0 ]; do case $1 in
  --target) TARGET=$2; shift 2;; --size) SIZE=$2; shift 2;; --only) ONLY=$2; shift 2;; --outdir) OUT=$2; shift 2;;
  *) echo "unknown arg $1"; exit 2;; esac; done
[ "$TARGET" = both ] && TARGET="gpu multicore"
mkdir -p "$OUT"; echo "app,target,size,status,rc" > "$OUT/summary.csv"
T="$HERE/../opencl/test"; D="$HERE/../data/inputs"
declare -A DIR ARGS_small ARGS_medium
reg() { DIR[$1]=$2; ARGS_small[$1]="$3"; ARGS_medium[$1]="$4"; }
# same arguments as opencl/run_all.sh
reg lud     dense-linear-algebra/lud     "-i $T/dense-linear-algebra/lud/256.dat -v"  "-i $T/dense-linear-algebra/lud/2048.dat"
reg spmv    sparse-linear-algebra/spmv   "-i $D/csr_small.dat -a"   "-i $D/csr_medium.dat -a"
reg bfs     graph-traversal/bfs          "$T/graph-traversal/bfs/graph4096.txt"  "$T/graph-traversal/bfs/graph65536.txt"
reg crc     combinational-logic/crc      "-i $D/crc_small.dat -a"   "-i $D/crc_medium.dat -a"
reg nw      dynamic-programming/nw       "1024 10"  "4096 10"
reg srad    structured-grids/srad        "512 512 0 127 0 127 0.5 2"  "2048 2048 0 127 0 127 0.5 2"
reg cfd     unstructured-grids/cfd       "$T/unstructured-grids/cfd/fvcorr.domn.097K"  "$T/unstructured-grids/cfd/fvcorr.domn.193K"
reg kmeans  dense-linear-algebra/kmeans  "-i $T/dense-linear-algebra/kmeans/204800.txt -r -o"  "-i $T/dense-linear-algebra/kmeans/819200.txt -r -o"
reg nqueens branch-and-bound/nqueens     "12"  "14"
reg hmm     graphical-models/hmm         "-n 60 -v n"  "-n 120 -v n"
reg swat    dynamic-programming/swat     "$T/dynamic-programming/swat/query1K1 $T/dynamic-programming/swat/sampledb1K1"  "$T/dynamic-programming/swat/query3K1 $T/dynamic-programming/swat/sampledb3K1"
reg tdm     finite-state-machine/tdm     "$T/finite-state-machine/tdm/sim-64-size50.csv $T/finite-state-machine/tdm/ivl.txt $T/finite-state-machine/tdm/30-episodes.txt 128"  "$T/finite-state-machine/tdm/sim-64-size500.csv $T/finite-state-machine/tdm/ivl.txt $T/finite-state-machine/tdm/30-episodes.txt 128"
reg gem     n-body-methods/gem-dwarf     "$T/n-body-methods/gem/1uwo_A 80 1 0"  "$T/n-body-methods/gem/Mb.Hhelix.bondi 80 1 0"
reg fft     spectral-methods/fft         "--pts 1 --passes 5"  "--pts 4 --passes 5"
APPS="lud spmv bfs crc nw srad cfd kmeans nqueens hmm swat tdm gem fft"
[ -n "$ONLY" ] && APPS="${ONLY//,/ }"
overall=0
for tgt in $TARGET; do
  for a in $APPS; do
    d="${DIR[$a]:-}"; [ -n "$d" ] || { echo "unknown app $a"; overall=1; continue; }
    eval "args=\${ARGS_${SIZE}[$a]}"
    log="$OUT/$a.$tgt.$SIZE.log"; rc=NA
    if ! make -C "$HERE/$d" build -j8 ACC_TARGET=$tgt > "$OUT/$a.$tgt.build.log" 2>&1; then st=BLOCKED_BUILD
    else
      make -s -C "$HERE/$d" run ACC_TARGET=$tgt GPU_INDEX=$GPU_INDEX RUN_ARGS="$args" > "$log" 2>&1; rc=$?
      line=$(grep -a '^DWARF_RESULT' "$log" | tail -1)
      if [ -z "$line" ]; then st=BLOCKED_RUNTIME
      elif [ $rc = 0 ] && echo "$line" | grep -q 'valid=1' && ! grep -a -q '^FAILURE' "$log"; then st=PASS
      else st=FAIL; fi
    fi
    [ $st = PASS ] || overall=1
    echo "$a,$tgt,$SIZE,$st,$rc" | tee -a "$OUT/summary.csv"
  done
done
exit $overall
