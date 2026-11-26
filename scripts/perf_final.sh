#!/bin/bash
# FINAL cross-model timed comparison (pinned configuration; omp/results/PROTOCOL-final.md).
#   scripts/perf_final.sh --outdir DIR [--cfgs "omp-clang omp-nvhpc opencl opencl-intel"] [--only a,b] [--launches N] [--devices "cpu gpu"] [--end-epoch SECS] [--reserve SECS]
# * CPU: every configuration (omp-clang, omp-nvhpc, PoCL) uses the SAME 128 threads on the SAME cpuset (64 physical cores x 2 SMT).
# * GPU: every configuration uses the pinned launch configuration of common/dwarf_launch.h (see LAUNCH-CONFIGS.md).
# * Resumable: an app whose marker DIR/done/<app> exists is skipped; raw.csv is appended to. Exit code 75 = time ran out, work remains
#   (a Slurm chain job resubmits itself on 75). --end-epoch is the job end time; no new app is started with less than --reserve seconds left.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$HERE/toolchains/env.sh"
OUT="$HERE/omp/results/perf-final-$(date +%Y%m%d-%H%M%S)"; ONLY=""; LAUNCHES=3; DEVICES="cpu gpu"; DEADLINE=""; END_EPOCH=""; RESERVE=900; CFGS_ARG=""; THREADS=128; REPS=10; WARM=2; TMO=900
while [ $# -gt 0 ]; do case $1 in --outdir) OUT=$2; shift 2;; --only) ONLY=$2; shift 2;; --launches) LAUNCHES=$2; shift 2;; --devices) DEVICES=$2; shift 2;; --deadline) DEADLINE=$2; shift 2;; --end-epoch) END_EPOCH=$2; shift 2;; --reserve) RESERVE=$2; shift 2;; --cfgs) CFGS_ARG=$2; shift 2;; --threads) THREADS=$2; shift 2;; --reps) REPS=$2; shift 2;; --warmup) WARM=$2; shift 2;; --timeout) TMO=$2; shift 2;; *) shift;; esac; done
mkdir -p "$OUT/done"; cp "$HERE/omp/results/PROTOCOL-final.md" "$OUT/" 2>/dev/null; cp "$HERE/omp/results/LAUNCH-CONFIGS.md" "$OUT/" 2>/dev/null
NTHR=$THREADS; CPUSET=$(python3 -I -c "import sys; n=int(sys.argv[1]); c=[x for a,b in ((0,63),(96,159)) for x in range(a,b+1)]; print(','.join(map(str,c[:n])) if n<128 else '0-63,96-159')" $NTHR)
TO="$HERE/omp/test"; TC="$HERE/opencl/test"; D="$HERE/data/inputs"; OCLP="$OMP_PREFIX"
declare -A ODIR CDIR OARGS CARGS
reg() { ODIR[$1]=$2; CDIR[$1]=$3; OARGS[$1]="$4"; CARGS[$1]="$5"; }
reg lud     dense-linear-algebra/lud       dense-linear-algebra/lud      "RUN_INPUT=$TO/dense-linear-algebra/lud/2048.dat"   "-i $TC/dense-linear-algebra/lud/2048.dat"
reg spmv    sparse-linear-algebra/spmv     sparse-linear-algebra/SPMV    "RUN_INPUT=$D/csr_medium.dat"                        "-i $D/csr_medium.dat -a"
reg bfs     graph-traversal/bfs            graph-traversal/bfs           "RUN_INPUT=$TO/graph-traversal/bfs/graph65536.txt"   "$TC/graph-traversal/bfs/graph65536.txt"
reg crc     combinational-logic/crc        combinational-logic/crc       "RUN_INPUT=$D/crc_medium.dat"                        "-i $D/crc_medium.dat -a"
reg nw      dynamic-programming/nw         dynamic-programming/nw        "RUN_SIZE=4096"                                       "4096 10"
reg srad    structured-grids/srad          structured-grids/srad         "RUN_SIZE=2048"                                       "2048 2048 0 127 0 127 0.5 2"
reg cfd     unstructured-grids/cfd         unstructured-grids/cfd        "RUN_INPUT=$TO/unstructured-grids/cfd/fvcorr.domn.193K" "$TC/unstructured-grids/cfd/fvcorr.domn.193K"
reg kmeans  dense-linear-algebra/kmeans    dense-linear-algebra/kmeans   "RUN_INPUT=$TO/dense-linear-algebra/kmeans/819200.txt"  "-i $TC/dense-linear-algebra/kmeans/819200.txt -r -o"
reg nqueens branch-and-bound/nqueens       branch-and-bound/nqueens      "RUN_SIZE=14"                                         "14"
reg hmm     graphical-models/hmm           graphical-models/hmm          "RUN_SIZE=120"                                        "-n 120 -v n"
reg swat    dynamic-programming/swat       dynamic-programming/swat      "RUN_INPUT=$TO/dynamic-programming/swat/query3K1 RUN_INPUT2=$TO/dynamic-programming/swat/sampledb3K1" "$TC/dynamic-programming/swat/query3K1 $TC/dynamic-programming/swat/sampledb3K1"
reg tdm     finite-state-machine/tdm       finite-state-machine/tdm      "RUN_INPUT=$TO/finite-state-machine/tdm/sim-64-size500.csv" "$TC/finite-state-machine/tdm/sim-64-size500.csv $TC/finite-state-machine/tdm/ivl.txt $TC/finite-state-machine/tdm/30-episodes.txt 128"
reg gem     n-body-methods/gem             n-body-methods/gem            "RUN_INPUT=$TO/n-body-methods/gem/Mb.Hhelix.bondi"    "$TC/n-body-methods/gem/Mb.Hhelix.bondi 80 1 0"
reg fft     spectral-methods/fft           spectral-methods/fft          "RUN_SIZE=4"                                          "--pts 4 --passes 5"
# lud last: its x86 OpenMP cells take minutes per launch (scheduling only; see the protocol addendum)
APPS="spmv bfs crc nw srad cfd kmeans nqueens hmm swat tdm gem fft lud"; [ -n "$ONLY" ] && APPS="${ONLY//,/ }"
GPU=$GPU_INDEX
# ---- environment snapshot ----
{ echo "threads=$NTHR cpuset=$CPUSET launch_LS=128 reps=$REPS warmup=$WARM timeout=$TMO cfgs=${CFGS_ARG:-default}"; ls "$INTEL_OCL_LIB"/libintelocl.so.* 2>/dev/null | head -1 | sed "s/^/intel-ocl: /"; date; hostname; nvidia-smi --query-gpu=name,driver_version,compute_cap --format=csv,noheader | head -2; lscpu | grep -E 'Model name|^CPU\(s\)|Socket|Thread|L3'
  $OMP_PREFIX/bin/clang --version | head -1; $NVHPC_ROOT/compilers/bin/nvc++ --version | sed -n 2p; echo "PoCL: $OCL_PREFIX"; ls $OCL_PREFIX/lib | grep 'libpocl.so\.' | head -1
  echo "GPU_NAME=$GPU_NAME GPU_INDEX=$GPU_INDEX NVHPC=$NVHPC_ROOT"; } >> "$OUT/environment.txt" 2>&1
snap() { local l f o; o=$("$HERE/scripts/system_clear.sh" --gpu $GPU --no-load); l=$(echo "$o" | sed -n 's/.*foreign_cores=\([0-9.]*\).*/\1/p'); f=$(echo "$o" | sed -n 's/.*foreign_procs=\([0-9]*\).*/\1/p'); echo "$l:$f"; }
[ -s "$OUT/raw.csv" ] || echo "app,device,config,launch,status,median_ms,mean_ms,min_ms,max_ms,foreign_cores_before,foreign_gpu_before,foreign_cores_after,foreign_gpu_after" > "$OUT/raw.csv"
# ---- build everything first ----
for a in $APPS; do
  [ -e "$OUT/done/$a" ] && continue
  for tc in clang nvhpc; do for dev in cpu $GPU_NAME; do make -C "$HERE/omp/${ODIR[$a]}" build DEVICE=$dev TOOLCHAIN=$tc >"$OUT/$a.omp-$tc.$dev.build.log" 2>&1 || echo "BLOCKED_BUILD omp-$tc $a $dev" | tee -a "$OUT/blocked.txt"; done; done
  make -C "$HERE/opencl/${CDIR[$a]}" build >"$OUT/$a.opencl.build.log" 2>&1 || echo "BLOCKED_BUILD opencl $a" | tee -a "$OUT/blocked.txt"
done
skip=" "
run_cfg() {  # app dev(cpu|gpu) cfg launch
  local a=$1 dv=$2 cfg=$3 L=$4 key="$1/$2/$3" log="$OUT/$1.$3.$2.$4.log" rc dev b af line st
  [ "$cfg" = opencl-intel ] && [ "$dv" = gpu ] && return   # Intel CPU OpenCL runtime: CPU cells only
  case "$skip" in *" $key "*) return;; esac
  if [ -n "$DEADLINE" ] && [ "$(date +%s)" -ge "$(date -d "$DEADLINE" +%s)" ]; then echo "$a,$dv,$cfg,$L,SKIPPED_TIME,,,,,,,," >> "$OUT/raw.csv"; return; fi
  grep -q "BLOCKED_BUILD ${cfg%%.*} $a" "$OUT/blocked.txt" 2>/dev/null && return
  local waited=0; until "$HERE/scripts/system_clear.sh" --gpu $GPU --no-load --max-foreign-cores 8 --quiet; do sleep 20; waited=$((waited+20)); [ $waited -ge 600 ] && break; done
  if [ $waited -ge 600 ]; then echo "$a,$dv,$cfg,$L,SKIPPED_BUSY,,,,,,,," >> "$OUT/raw.csv"; return; fi
  b=$(snap)
  local pin="taskset -c $CPUSET"
  case $cfg in
    omp-clang|omp-nvhpc)
      local tc=${cfg#omp-} dd=cpu off=DISABLED env_extra="OMP_NUM_THREADS=$NTHR OMP_PROC_BIND=close OMP_PLACES=threads"
      [ $dv = gpu ] && { dd=$GPU_NAME; off=MANDATORY; env_extra="OMP_NUM_THREADS=$NTHR"; }
      env CUDA_VISIBLE_DEVICES=$GPU OMP_TARGET_OFFLOAD=$off OMP_THREAD_LIMIT=2048 $env_extra $pin timeout $TMO make -C "$HERE/omp/${ODIR[$a]}" run DEVICE=$dd TOOLCHAIN=$tc RUN_REPEATS=$REPS RUN_WARMUP=$WARM ${OARGS[$a]} >"$log" 2>&1; rc=$?;;
    opencl|opencl-intel)
      local vend=pocl; [ $dv = gpu ] && vend=nvidia; [ $cfg = opencl-intel ] && vend=intel
      env POCL_MAX_PTHREAD_COUNT=$NTHR OMP_NUM_THREADS=$NTHR $pin timeout $TMO make -C "$HERE/opencl/${CDIR[$a]}" run OCL_VENDOR=$vend GPU_INDEX=$GPU DWARF_REPEATS=$REPS DWARF_WARMUP=$WARM RUN_ARGS="${CARGS[$a]}" >"$log" 2>&1; rc=$?
      dev=$(grep -a -m1 'Device Chosen' "$log")
      if [ $vend != nvidia ]; then echo "$dev" | grep -q -i NVIDIA && rc=97; else echo "$dev" | grep -q -i NVIDIA || rc=97; fi;;
  esac
  af=$(snap); line=$(grep -a '^DWARF_RESULT' "$log" | tail -1)
  f() { echo "$line" | tr ',' '\n' | sed -n "s/^$1=//p"; }
  if [ $rc = 97 ]; then st=WRONG_DEVICE; skip="$skip$key "
  elif [ $rc = 124 ]; then st=TIMEOUT; skip="$skip$key "
  elif [ -z "$line" ]; then st=BLOCKED_RUNTIME; skip="$skip$key "
  elif [ "$(f valid)" != 1 ]; then st=FAIL_VALIDATION
  else st=OK; fi
  local lb=${b%%:*} fb=${b##*:} la=${af%%:*} fa=${af##*:}
  if [ $st = OK ] && { [ "${fa:-0}" -gt 0 ] || awk -v l="$la" 'BEGIN{exit !(l>=8)}'; }; then st=CONTAMINATED; fi
  echo "$a,$dv,$cfg,$L,$st,$(f median_ms),$(f mean_ms),$(f min_ms),$(f max_ms),$lb,$fb,$la,$fa" | tee -a "$OUT/raw.csv"
}
CFGS=(omp-clang omp-nvhpc opencl); [ -n "$CFGS_ARG" ] && read -r -a CFGS <<< "$CFGS_ARG"
rc_all=0
for a in $APPS; do
  [ -e "$OUT/done/$a" ] && continue
  if [ -n "$END_EPOCH" ] && [ $(( END_EPOCH - $(date +%s) )) -lt "$RESERVE" ]; then echo "TIME: not starting $a (less than ${RESERVE}s left)"; rc_all=75; continue; fi
  # an app without a done marker was interrupted: drop its partial rows so no launch is counted twice
  grep -v "^$a," "$OUT/raw.csv" > "$OUT/raw.csv.tmp"; mv "$OUT/raw.csv.tmp" "$OUT/raw.csv"
  for dv in $DEVICES; do for L in $(seq 1 $LAUNCHES); do
    k=$(( L % ${#CFGS[@]} )); order=("${CFGS[@]:$k}" "${CFGS[@]:0:$k}")
    for cfg in "${order[@]}"; do run_cfg $a $dv $cfg $L; done
  done; done
  touch "$OUT/done/$a"
done
python3 -I "$HERE/scripts/perf_final_summarize.py" "$OUT" || true
exit $rc_all
