#!/usr/bin/env bash
# Generic multi-device sweep over the 7-tier cache-bounding problem sizes.
# Builds (or reuses) each dwarf's binary for the requested targets and runs
# it at every requested tier, appending DWARF_RESULT lines into a
# summary.csv. This is the "generate the data" half of the 7-tier
# methodology: run compute_cache_tier_sizes.py first to produce the sizes
# manifest, then point this script at it.
#
# No hostnames or hardcoded toolchain paths -- point CXX/ROCM_PATH/etc. at
# whatever this machine's offload compiler is via environment variables.
#
# Usage:
#   CXX=/path/to/clang++ \
#   DWARF_TIERS_MANIFEST=/path/to/cache_tier_sizes_7tier.json \
#     bash run_cache_tier_sweep.sh
#
#   DWARF_TARGETS="cpu a100"           bash run_cache_tier_sweep.sh   # subset of targets
#   DWARF_TIERS="l1lo l1hi stream"     bash run_cache_tier_sweep.sh   # subset of tiers
#   DWARF_FILTER="atax spmv"           bash run_cache_tier_sweep.sh   # subset of dwarfs
#   DWARF_BUILD_DIR=build-custom      bash run_cache_tier_sweep.sh   # per-target build dir
#   DWARF_GPU_DEV=1                    bash run_cache_tier_sweep.sh   # pick a GPU index
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
host="$(hostname -s)"
timestamp="$(date +%Y%m%d_%H%M%S_%N)"
run_id="${timestamp}__cache_tier_sweep__${host}"
out_dir="${root_dir}/results/${run_id}"
mkdir -p "${out_dir}"

: "${CXX:?set CXX to this machine OpenMP-target-offload clang++}"
DWARF_TIERS_MANIFEST="${DWARF_TIERS_MANIFEST:?set DWARF_TIERS_MANIFEST to a manifest from compute_cache_tier_sizes.py}"
[[ -f "${DWARF_TIERS_MANIFEST}" ]] || { echo "ERROR: manifest ${DWARF_TIERS_MANIFEST} not found" >&2; exit 1; }

DWARF_BUILD_DIR="${DWARF_BUILD_DIR:-build}"
DWARF_TIERS="${DWARF_TIERS:-l1lo l1hi l2lo l2hi llclo llchi stream}"
DWARF_TARGETS="${DWARF_TARGETS:-cpu}"
DWARF_FILTER="${DWARF_FILTER:-}"
DWARF_WARMUP="${DWARF_WARMUP:-1}"
DWARF_REPEATS="${DWARF_REPEATS:-3}"
DWARF_GPU_DEV="${DWARF_GPU_DEV:-0}"
DWARF_TIMEOUT_S="${DWARF_TIMEOUT_S:-600}"

export CXX

# ---------- all supported benchmarks (including empirical seeds) ------------
benchmarks=(
  "branch-and-bound/nqueens:nqueens"
  "combinational-logic/crc:crc"
  "dense-linear-algebra/gemm:gemm"
  "dense-linear-algebra/atax:atax"
  "dense-linear-algebra/cg:cg"
  "dense-linear-algebra/kmeans:kmeans"
  "dense-linear-algebra/lud:lud"
  "dynamic-programming/nw:nw"
  "dynamic-programming/swat:swat"
  "finite-state-machine/tdm:tdm"
  "graphical-models/hmm:hmm"
  "structured-grids/fdtd2d:fdtd2d"
  "structured-grids/srad:srad"
  "sparse-linear-algebra/spmv:spmv"
  "spectral-methods/fft:fft"
  "graph-traversal/bfs:bfs"
  "n-body-methods/gem:gem"
  "unstructured-grids/cfd:cfd"
  "mapreduce/pagerank:pagerank"
)

# ---------- load sizes from manifest -----------------------------------------
declare -A SIZE SEC ITERS
assignments="$(python3 "${root_dir}/scripts/_manifest_to_bash_arrays.py" "${DWARF_TIERS_MANIFEST}")"
eval "${assignments}"

for target in ${DWARF_TARGETS}; do
  case "$target" in cpu|a100|h100|h200|mi100|mi210|mi300a) ;;
    *) echo "ERROR: unknown target: $target" >&2; exit 2 ;;
  esac
done
for tier in ${DWARF_TIERS}; do
  case "$tier" in l1lo|l1hi|l2lo|l2hi|llclo|llchi|stream) ;;
    *) echo "ERROR: unknown tier: $tier" >&2; exit 2 ;;
  esac
done
for name in ${DWARF_FILTER}; do
  found=0
  for entry in "${benchmarks[@]}"; do
    [[ "${entry##*:}" == "$name" ]] && found=1
  done
  if [[ "$found" == 0 ]]; then
    echo "ERROR: unknown benchmark: $name" >&2; exit 2
  fi
done
failed=0
run_count=0

{
  echo "run_id=${run_id}"
  echo "host=${host}"
  echo "cxx=${CXX}"
  echo "build_dir=${DWARF_BUILD_DIR}"
  echo "tiers=${DWARF_TIERS}"
  echo "targets=${DWARF_TARGETS}"
  echo "filter=${DWARF_FILTER}"
  echo "warmup_runs=${DWARF_WARMUP}"
  echo "measured_runs=${DWARF_REPEATS}"
  echo "manifest=${DWARF_TIERS_MANIFEST}"
  echo "date=$(date --iso-8601=seconds)"
} > "${out_dir}/manifest.txt"

summary="${out_dir}/summary.csv"
echo "tier,benchmark,motif,target,size,secondary_size,iterations,repeats,mean_ms,median_ms,valid" > "${summary}"

device_env_for() {
  case "$1" in
    a100|h100|h200) echo "CUDA_VISIBLE_DEVICES=${DWARF_GPU_DEV}" ;;
    mi100|mi210|mi300a) echo "ROCR_VISIBLE_DEVICES=${DWARF_GPU_DEV} HIP_VISIBLE_DEVICES=${DWARF_GPU_DEV}" ;;
    *) echo "" ;;
  esac
}

run_one() {
  local tier="$1" bench_path="$2" bench_short="$3" target="$4" p_size="$5" p_sec="$6" p_it="$7"
  run_count=$((run_count + 1))
  local log="${out_dir}/${tier}__${bench_short}__${target}.log"
  echo "  [${tier}] ${bench_short} ${target} (size=${p_size} sec=${p_sec} iters=${p_it})" | tee -a "${out_dir}/run.log"

  local dev_env; dev_env="$(device_env_for "${target}")"
  # shellcheck disable=SC2086
  if env ${dev_env} timeout "${DWARF_TIMEOUT_S}" make -C "${root_dir}/${bench_path}" run DEVICE="${target}" \
      BUILD_DIR="${DWARF_BUILD_DIR}" CXX="${CXX}" \
      RUN_SIZE="${p_size}" RUN_SECONDARY_SIZE="${p_sec}" RUN_ITERATIONS="${p_it}" \
      RUN_WARMUP="${DWARF_WARMUP}" RUN_REPEATS="${DWARF_REPEATS}" \
      > "${log}" 2>&1; then
    if ! grep -Eq '^DWARF_RESULT,.*valid=1(,|$)' "${log}"; then
      echo "    FAILED validation; see ${log}" | tee -a "${out_dir}/run.log"
      failed=1
      echo "${tier},${bench_short},??,${target},${p_size},${p_sec},${p_it},${DWARF_REPEATS},,,0" >> "${summary}"
      return
    fi
    grep '^DWARF_RESULT,' "${log}" | sed 's/^DWARF_RESULT,//' | awk -F, -v t="${tier}" '
      { delete v; for (i = 1; i <= NF; ++i) { split($i, kv, "="); v[kv[1]] = kv[2] }
        print t "," v["benchmark"] "," v["motif"] "," v["target"] "," v["size"] "," \
              v["secondary_size"] "," v["iterations"] "," v["repeats"] "," \
              v["mean_ms"] "," v["median_ms"] "," v["valid"] }' >> "${summary}"
  else
    local rc=$?; failed=1; local reason="FAILED"; [[ $rc -eq 124 ]] && reason="TIMEOUT(${DWARF_TIMEOUT_S}s)"
    echo "    ${reason} ${bench_short} ${target} [${tier}]; see ${log}" | tee -a "${out_dir}/run.log"
    echo "${tier},${bench_short},??,${target},${p_size},${p_sec},${p_it},${DWARF_REPEATS},,,0" >> "${summary}"
  fi
}

for tier in ${DWARF_TIERS}; do
  echo "=== Cache Tier: ${tier} ===" | tee -a "${out_dir}/run.log"
  for entry in "${benchmarks[@]}"; do
    bench_path="${entry%%:*}"; bench_short="${entry##*:}"

    if [[ -n "${DWARF_FILTER}" ]]; then
      match=0
      for f in ${DWARF_FILTER}; do [[ "${bench_short}" == "${f}" ]] && match=1 && break; done
      [[ "${match}" -eq 0 ]] && continue
    fi

    p_size="${SIZE[${bench_short}:${tier}]:-}"
    p_sec="${SEC[${bench_short}:${tier}]:-0}"
    p_it="${ITERS[${bench_short}:${tier}]:-1}"

    if [[ -z "${p_size}" ]]; then
      echo "  [${tier}] ${bench_path} -- no size in manifest, skipping" | tee -a "${out_dir}/run.log"
      continue
    fi

    for target in ${DWARF_TARGETS}; do
      run_one "${tier}" "${bench_path}" "${bench_short}" "${target}" "${p_size}" "${p_sec}" "${p_it}"
    done
  done
done

echo ""
echo "Run complete. Results: ${out_dir}"
echo "${out_dir}"

if [[ "$run_count" == 0 ]]; then
  echo "ERROR: no matching benchmark sizes in manifest" >&2
  exit 1
fi
exit "$failed"
