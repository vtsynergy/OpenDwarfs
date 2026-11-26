#!/bin/bash
# Validate every SYCL app against the OpenCL tree (ground truth) with small inputs (those of opencl/run_all.sh).
#   sycl/tools/validate.sh [--devices "cpu gpu"] [--only a,b] [--size small|medium] [--outdir DIR] [--cpuset 64-95] [--no-ref]
# Per app and device: (1) SYCL build, (2) SYCL run: the app's own validator (valid=1 in DWARF_RESULT),
# (3) OpenCL reference run on the same device from a COPY of the prebuilt opencl/<app>/build (nothing is built or
# written under opencl/), with tools/ocl_dump.so recording every read-back, (4) comparison of the two streams.
# CPU reference = Intel CPU OpenCL runtime (same runtime as the SYCL CPU backend); GPU reference = NVIDIA OpenCL.
# Status: PASS (validator and comparison), PASS_VALIDATOR_ONLY (no reference), DIFF, FAIL, BLOCKED_BUILD, BLOCKED_RUNTIME.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; REPO="$(dirname "$HERE")"
source "$REPO/toolchains/env.sh" >/dev/null
export ONEAPI_ROOT="${ONEAPI_ROOT:-$INTEL_ONEAPI_ROOT}"
source "$HERE/benchmarks.sh"
[ -n "${BENCH_EXTRA:-}" ] && source "$BENCH_EXTRA"   # optional file overriding ARGS_/CMP entries (extra sizes)
DEVICES="cpu gpu"; ONLY=""; SIZE=small; OUT="$HERE/results/validate-$(date +%Y%m%d-%H%M%S)"; CPUSET=""; NOREF=0
while [ $# -gt 0 ]; do case $1 in
  --devices) DEVICES=$2; shift 2;; --only) ONLY=$2; shift 2;; --size) SIZE=$2; shift 2;; --outdir) OUT=$2; shift 2;;
  --cpuset) CPUSET=$2; shift 2;; --no-ref) NOREF=1; shift;; *) echo "unknown arg $1"; exit 2;; esac; done
[ -n "$ONLY" ] && APPS="${ONLY//,/ }"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
PIN=""; [ -n "$CPUSET" ] && PIN="taskset -c $CPUSET"
gcc -O2 -fPIC -shared "$HERE/tools/ocl_dump.c" -o "$OUT/ocl_dump.so" -ldl || exit 1
[ -s "$OUT/summary.csv" ] || echo "app,device,status,sycl_exit,ref_exit,comparison" > "$OUT/summary.csv"
for a in $APPS; do
  eval "args=\${ARGS_${SIZE}[$a]}"; d=${DIR[$a]}
  if ! make -C "$HERE/$d" build -j8 > "$OUT/$a.build.log" 2>&1; then
    for dv in $DEVICES; do echo "$a,$dv,BLOCKED_BUILD,,," | tee -a "$OUT/summary.csv"; done; continue; fi
  for dv in $DEVICES; do
    sd="$OUT/$a.$dv.sycl.dump"; mkdir -p "$sd"
    $PIN timeout 600 make -s -C "$HERE/$d" run SYCL_DEVICE=$dv GPU_INDEX=$GPU_INDEX DWARF_DUMP_DIR="$sd" RUN_ARGS="$args" > "$OUT/$a.$dv.sycl.log" 2>&1; rc=$?
    line=$(grep -a '^DWARF_RESULT' "$OUT/$a.$dv.sycl.log" | tail -1)
    refrc=""; cmp=""
    if [ -z "$line" ]; then st=BLOCKED_RUNTIME
    elif ! grep -a -q ",device=$dv," <<< "$line" || ! grep -a -q "SYCL requested=$dv " "$OUT/$a.$dv.sycl.log"; then st=FAIL; cmp="wrong device"
    elif ! grep -a -q ',valid=1' <<< "$line" || [ $rc != 0 ]; then st=FAIL
    else st=PASS_VALIDATOR_ONLY; fi
    if [ $st = PASS_VALIDATOR_ONLY ] && [ $NOREF = 0 ] && [ -d "$REPO/opencl/${OCLDIR[$a]}/build" ]; then
      rd="$OUT/$a.$dv.ocl.run"; od="$OUT/$a.$dv.ocl.dump"; mkdir -p "$rd/icd" "$od"
      find "$REPO/opencl/${OCLDIR[$a]}/build" -maxdepth 1 -type f -exec cp {} "$rd/" \;
      if [ $dv = cpu ]; then echo "$INTEL_OCL_LIB/libintelocl.so" > "$rd/icd/ref.icd"; lp="$INTEL_OCL_LIB:$INTEL_TBB_LIB"
      else echo libnvidia-opencl.so.1 > "$rd/icd/ref.icd"; lp=""; fi
      ( cd "$rd" && env OCL_ICD_VENDORS="$rd/icd" LD_LIBRARY_PATH="$lp" CUDA_VISIBLE_DEVICES=$GPU_INDEX DWARF_DUMP_DIR="$od" LD_PRELOAD="$OUT/ocl_dump.so" \
          $PIN timeout 600 ./"${BINNAME[$a]}" -p 0 -d 0 -- $args > "$OUT/$a.$dv.ocl.log" 2>&1 ); refrc=$?
      chosen=$(grep -a -m1 'Device Chosen' "$OUT/$a.$dv.ocl.log")
      isnv=0; grep -qi nvidia <<< "$chosen" && isnv=1
      if [ $refrc != 0 ] || ! grep -a -q '^DWARF_RESULT.*valid=1' "$OUT/$a.$dv.ocl.log"; then cmp="reference run failed"
      elif { [ $dv = cpu ] && [ $isnv = 1 ]; } || { [ $dv = gpu ] && [ $isnv = 0 ]; }; then cmp="reference on wrong device"
      elif [[ ${CMP[$a]} == "block "* ]]; then
        # comparison of whole output blocks (sed address range), e.g. swat's printed alignments and scores
        re=${CMP[$a]#block }; x=$(tr -d '\000' < "$OUT/$a.$dv.ocl.log" | sed -n "$re"); y=$(tr -d '\000' < "$OUT/$a.$dv.sycl.log" | sed -n "$re")
        if [ -n "$x" ] && [ "$x" = "$y" ]; then cmp="MATCH $(wc -l <<< "$y") identical output lines (sed -n '$re')"; st=PASS
        else cmp="DIFF output block"; st=DIFF; fi
      elif [[ ${CMP[$a]} == "text "* ]]; then
        # comparison of selected output lines (apps whose raw read-backs contain uninitialised device memory or print their result)
        re=${CMP[$a]#text }; x=$(grep -a -o -E "$re" "$OUT/$a.$dv.ocl.log"); y=$(grep -a -o -E "$re" "$OUT/$a.$dv.sycl.log")
        if [ -n "$x" ] && [ "$x" = "$y" ]; then cmp="MATCH output lines: $(head -c 200 <<< "$y" | tr '\n' ';')"; st=PASS
        else cmp="DIFF output lines"; st=DIFF; fi
      else
        cmp=$(python3 -I "$HERE/tools/compare_readback.py" "$od/readback.bin" "$sd/readback.bin" ${CMP[$a]}); crc=$?
        if [ $crc = 0 ]; then st=PASS; else st=DIFF; fi
      fi
    fi
    echo "$a,$dv,$st,$rc,$refrc,\"$cmp\"" | tee -a "$OUT/summary.csv"
  done
done
