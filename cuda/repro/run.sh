#!/bin/bash
# CUDA-WA-01 reproducer.  ./run.sh [outdir]   (default ./build; NVCC and GPU_INDEX overridable)
# Expected with nvcc 12.9.86: -O3 FAILS; -O3 -Xcicc -O0 and -O3 -DWORKAROUND PASS.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NVCC="${NVCC:-/opt/nvidia/hpc_sdk/Linux_x86_64/25.7/cuda/12.9/bin/nvcc}"
OUT="${1:-$HERE/build}"; mkdir -p "$OUT"
"$NVCC" --version | tail -2
i=0
for flags in "-O3" "-O3 -Xcicc -O0" "-O3 -Xptxas -O0" "-O3 -DWORKAROUND"; do
  i=$((i+1)); echo "== nvcc -w -gencode arch=compute_90,code=sm_90 $flags cuda_wa01.cu"
  "$NVCC" -w -gencode arch=compute_90,code=sm_90 $flags "$HERE/cuda_wa01.cu" -o "$OUT/wa01_$i" && CUDA_VISIBLE_DEVICES="${GPU_INDEX:-0}" "$OUT/wa01_$i"
done
"$NVCC" -w -O3 -arch=sm_90 -ptx "$HERE/cuda_wa01.cu" -o "$OUT/cuda_wa01.ptx"
echo "== PTX evidence ($OUT/cuda_wa01.ptx): stores with a +64 byte offset, and the +16 element base they are derived from"
grep -n 'add.s32.*, 16;\|st.global.f32.*+64\|ld.global.f32.*+64' "$OUT/cuda_wa01.ptx"
