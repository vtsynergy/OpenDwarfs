#!/bin/bash
# =============================================================================
# Run every OpenCL benchmark once (or REPEATS times) and log results.
#
#   ./run_all.sh [--vendor nvidia|pocl|intel] [--size small|medium] [--repeats N]
#                [--only app1,app2] [--build] [--outdir DIR]
#
# Environment: GPU_INDEX (default 0, see toolchains/env.sh), OCL_PREFIX.
# Output: <outdir>/summary.csv (app,vendor,size,rep,status,wall_s) + per-run logs.
# Every app is launched through its own Makefile (`make run`), so the exact same
# path is used when you run a single benchmark by hand.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../toolchains/env.sh"
VENDOR=nvidia; SIZE=small; REPEATS=1; ONLY=""; DO_BUILD=0
OUT="${HERE}/results/$(date +%Y%m%d-%H%M%S)"
while [ $# -gt 0 ]; do case $1 in
  --vendor) VENDOR=$2; shift 2;; --size) SIZE=$2; shift 2;; --repeats) REPEATS=$2; shift 2;;
  --only) ONLY=$2; shift 2;; --build) DO_BUILD=1; shift;; --outdir) OUT=$2; shift 2;;
  *) echo "unknown arg $1"; exit 2;; esac; done
T="$HERE/test"
mkdir -p "$OUT"; echo "app,vendor,size,rep,status,wall_s" > "$OUT/summary.csv"

# app | directory | "args for make RUN_ARGS" per size.  Paths are absolute.
# Sizes chosen to match omp/ where an omp counterpart exists (see ../scripts/compare_map.csv).
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
reg gem     n-body-methods/gem              "$T/n-body-methods/gem/1uwo_A 80 1 0"  "$T/n-body-methods/gem/Mb.Hhelix.bondi 80 1 0"
reg srad    structured-grids/srad           "512 512 0 127 0 127 0.5 2"            "2048 2048 0 127 0 127 0.5 2"
reg cfd     unstructured-grids/cfd          "$T/unstructured-grids/cfd/fvcorr.domn.097K" "$T/unstructured-grids/cfd/fvcorr.domn.193K"
# crc and spmv need generated inputs (deterministic, -R): @D@ = $OUT/data
reg crc     combinational-logic/crc         "-i @D@/crc_small.dat -a"              "-i @D@/crc_medium.dat -a"
reg spmv    sparse-linear-algebra/SPMV      "-i @D@/csr_small.dat -a"              "-i @D@/csr_medium.dat -a"
reg fft     spectral-methods/fft            "--pts 1 --passes 5"                    "--pts 3 --passes 5"

APPS="nqueens crc lud kmeans nw swat tdm bfs hmm gem spmv srad cfd fft"
[ -n "$ONLY" ] && APPS="${ONLY//,/ }"
gen_inputs() {  # builds generators on demand, creates inputs once
  mkdir -p "$OUT/data"
  make -C "$HERE/combinational-logic/crc" build -j8 >/dev/null 2>&1
  make -C "$HERE/sparse-linear-algebra/SPMV" build -j8 >/dev/null 2>&1
  [ -f "$OUT/data/crc_small.dat" ]  || "$HERE/combinational-logic/crc/build/createcrc" -R -n 1024  -s 1024 -f "$OUT/data/crc_small.dat"  >/dev/null
  [ -f "$OUT/data/crc_medium.dat" ] || "$HERE/combinational-logic/crc/build/createcrc" -R -n 16384 -s 1024 -f "$OUT/data/crc_medium.dat" >/dev/null
  [ -f "$OUT/data/csr_small.dat" ]  || "$HERE/sparse-linear-algebra/SPMV/build/createcsr" -R -n 2048 -d 5000 -f "$OUT/data/csr_small.dat"  >/dev/null
  [ -f "$OUT/data/csr_medium.dat" ] || "$HERE/sparse-linear-algebra/SPMV/build/createcsr" -R -n 8192 -d 5000 -f "$OUT/data/csr_medium.dat" >/dev/null
}
gen_inputs
for a in $APPS; do
  d="${DIR[$a]:-}"; [ -n "$d" ] || { echo "unknown app $a"; continue; }
  [ "$DO_BUILD" = 1 ] && make -C "$HERE/$d" build -j8 >"$OUT/$a.build.log" 2>&1
  eval "args=\${ARGS_${SIZE}[$a]}"; args="${args//@D@/$OUT/data}"
  for r in $(seq 1 "$REPEATS"); do
    t0=$(date +%s.%N)
    make -C "$HERE/$d" run OCL_VENDOR=$VENDOR GPU_INDEX=$GPU_INDEX OCL_PREFIX=$OCL_PREFIX RUN_ARGS="$args" >"$OUT/$a.$VENDOR.$SIZE.$r.log" 2>&1
    rc=$?; t1=$(date +%s.%N)
    # guard: the implementation must be the one asked for (a CPU run must not land on the GPU platform or vice versa)
    dev=$(grep -a -m1 'Device Chosen' "$OUT/$a.$VENDOR.$SIZE.$r.log")
    if [ $VENDOR != nvidia ]; then echo "$dev" | grep -q -i 'NVIDIA' && rc=97; else echo "$dev" | grep -q -i 'NVIDIA' || rc=97; fi
    echo "$a,$VENDOR,$SIZE,$r,$([ $rc = 0 ] && echo OK || { [ $rc = 97 ] && echo WRONG_DEVICE || echo FAIL:$rc; }),$(echo "$t1 - $t0" | bc -l | xargs printf %.3f)" | tee -a "$OUT/summary.csv"
  done
done
