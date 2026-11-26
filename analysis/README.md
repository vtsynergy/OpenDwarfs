# analysis/ - comparative toolchain analysis (proof of concept)

Goal: when two toolchains that are nominally "the same thing" (Clang vs NVHPC OpenMP offload; PoCL vs Intel OpenCL on CPU) run the same program at different speeds,
find out **why**, from evidence, and keep the method reusable as toolchains change.

Status: proof of concept (GPU OpenMP pair validated on three apps; OpenCL interposer and CPU pair being trialled). Nothing here is a published result; see "Limits".

## Method (four steps)
1. **Decompose** each (app, toolchain) run into a normalised *time budget* (ms per repetition): kernel (device) time, launch, host-to-device copy, device-to-host copy, sync/wait,
   allocation, free, event bookkeeping, dispatch latency (OpenCL), init (excluded). Same buckets for every runtime, so budgets are comparable across toolchains.
2. **Attribute** the wall-time gap to buckets (delta per bucket) and run the **rule engine** (`compare.py`): each rule tests a measured signature (e.g. allocation calls per launch,
   sync time per launch, executed-instruction ratio, launch-shape mismatch) and emits a finding with its evidence, its contribution in ms/rep, and the intervention that would confirm it.
3. **Drill down** into the dominant bucket with lowered-code / counter evidence: per-kernel Nsight Compute counters (instructions executed, warps active, registers, DRAM bytes) for GPU;
   launch shape (grid x block / global x local) for any runtime. (Static SASS/PTX diff and CPU assembly/vectorisation diff are planned; `cuobjdump` returns nothing for Clang's offload images.)
4. **Validate by intervention** (not yet automated): a finding stays a *hypothesis* until the confirming intervention is run (change one knob / source setting) and the predicted gap closes.
   Controls: negative (same toolchain vs itself must give no cause findings), positive (known causes must be rediscovered blind).

## Files
| file | role |
|---|---|
| `apps.py` | app registry; re-uses run arguments from `scripts/perf_final.sh` |
| `gpu_collect.py` | run one OpenMP/CUDA-backed config under Nsight Systems (`--ncu`: per-kernel Nsight Compute counters) -> evidence bundle JSON |
| `clshim.c` / `ocl_collect.py` | LD_PRELOAD interposer for the OpenCL API (any ICD) + collector: launch shapes, per-API host time, event-profiled kernel time and queued->start latency |
| `budget.py` / `mark.c` | additive partition of the timed window; the timer-hook library (region marks + profiler capture) |
| `cpu_lowered.py` | PoCL vs Intel lowered-code diff for one app: runtime code dumps (PoCL kernel cache, Intel `CL_CONFIG_DUMP_ASM`), static x86 features per kernel, `perf stat` counters (cycles user/kernel, instructions, FP micro-ops by width) restricted to the timed repetitions via the timer hook |
| `crc_launchsize.py` | launch-size intervention on crc (pages per launch varied at fixed total work) for PoCL / Intel / Intel with the vectorizer off |
| `pocl_knobs.py` | PoCL code-generation knobs (vector width, work-group method, target CPU, specialization, vectorizer remarks) on kmeans and lud |
| `lud_memory.py` / `lud_chain.py` | lud cache-counter comparison; scratch-kernel dependency-chain / loop-structure test |
| `sweep_summary.py` | CPU thread-sweep summary (PROTOCOL-final.md Addendum 3) |
| `compare.py` | two bundles -> bucket deltas + ranked findings (markdown + JSON) |
| `launch_audit.py` | compares the grid x block actually launched by omp-clang and omp-nvhpc for every kernel |
| `lowered_gpu.py` | static SASS features from a binary (cuobjdump) - currently only works for NVHPC binaries |

## Evidence bundle (`<app>.<config>.bundle.json`)
`wall_ms` (DWARF timer median) | `buckets_ms_per_rep`, `bucket_counts_per_rep` | `kernel_ms_per_rep`, `launches_per_rep` | `kernels{name: ms_per_rep, grids, blocks, regs, ncu counters}` |
`api{name: count, total_ms_per_rep}` | `memcpy` | `lowered` | `toolchain`, `date`. Per-repetition values are process totals divided by (warm-up + repetitions); `init` calls are excluded from the budget.

## Rules implemented (`compare.py`)
K0/K1 kernel time same/different | C1 launch shape differs | L1 launch count differs | A1 alloc/free per launch in one toolchain only | E1 event bookkeeping | S1 sync time | D1 D2H call cost
(with/without pinned staging) | U1 GPU busy fraction | K2 executed instructions differ (Nsight Compute) | Q1 dispatch latency (OpenCL) | H1 host launch cost | O1 overhead factor.

## Results so far (2026-11-22, H100, `omp/results/analysis-poc-2026-11-22/`)
- crc, Clang vs NVHPC: kernel time equal (24.1 vs 21.9 ms/rep); the 102 ms gap is host-side: D2H copy calls +40 ms, sync +32 ms, alloc/free +9.5 ms, event records +4.8 ms (2049 vs 0 per rep);
  only NVHPC makes a pinned host allocation. Rediscovered the hand analysis blind.
- Negative control (crc, Clang vs Clang): 0.4% gap, no cause findings.
- swat: kernel time equal (29.1 vs 27.1 ms/rep); gap is sync (241 vs 41 ms/rep). This replaced an earlier hand-made claim that Clang's `k_match` kernel was slower (that was a launch-shape artefact before pinning).
- **Launch audit:** with the nominal 256-thread pin, 40 of 41 kernels launched different shapes under Clang and NVHPC (NVHPC caps flat kernels at 128 threads/block, rounds barrier kernels up to a warp multiple and adds a
  master warp). The pin was changed to 128; 35 of 41 now match (the rest are barrier kernels in lud, nw, fft, swat traceback).

- **CPU pair (PoCL vs Intel OpenCL, via `clshim`):** on crc and nw the gap is inside the event-measured kernel time (crc 309 vs 15 ms/rep; nw 212 vs 18 ms/rep), plus host launch cost on nw (161 vs 6 us per launch).
  Thread-scaling intervention (`omp/results/analysis-poc-2026-11-22/thread-scaling_pocl-vs-intel.md`): PoCL's per-launch kernel time grows with worker threads (crc 9 us at 1 thread -> 275 us at 128; wall 16.7 -> 381 ms) while Intel's stays flat (13-15 us),
  so a fixed 128-thread setting penalises PoCL on small launches. nqueens is the reverse case (Intel 2x slower): switching Intel's implicit vectorizer off (`CL_CONFIG_USE_VECTORIZER=false`) cut Intel's wall from 297 to 146 ms
  (kernel 282 -> 124 ms), closing the gap - the first *confirmed* attribution (hypothesis -> single-knob intervention -> predicted effect observed); crc was unaffected by the same knob (40 -> 36 ms).
- **Accounting fixed (2026-11-22):** the first version summed per-API-class times, which include the kernel time the host is waiting for, so deltas overlapped (crc PoCL-vs-Intel: sync +405 ms vs a +369 ms gap).
  `budget.py` now partitions the *timed* window into mutually exclusive states (device kernel running > device copy running > host in an API class with the device idle > host outside traced calls); the states sum exactly to the window and
  so do their deltas between two toolchains. The window is the timed repetitions only: the timer calls an optional hook (`dwarf_probe_mark`, looked up with dlsym; supplied by `libmark.so`), which logs region marks (OpenCL path) and calls
  `cuProfilerStart/Stop` so that Nsight captures only the timed repetitions (CUDA path). Check: window/wall coverage 1.10-1.11 on crc (GPU) and the sum of deltas equals the window gap by construction.
  New view this gives (crc, Clang vs NVHPC, share of the window gap): D2H copy waits +34%, sync +28%, host time outside CUDA calls (libomptarget bookkeeping) +17%, launch +5%, alloc/free/events +12%, kernels only +2%.

- **CPU lowered-code diff (`omp/results/analysis-poc-2026-11-22/lowered-cpu-scoped/FINDINGS.md`):** at one thread the PoCL-vs-Intel gaps fall into four classes: Intel's wider vectorization cuts instruction counts (kmeans 5.5x, lud 19.6x, gem 6.2x, cfd, swat);
  Intel's masked vector code loses IPC on branchy/irregular kernels (nqueens 151 -> 56 ms, tdm 327 -> 221 ms, bfs 8.3 -> 3.7 ms with `CL_CONFIG_USE_VECTORIZER=false`; all three confirmed, instruction count roughly unchanged for nqueens and bfs, 1.7x higher for scalar tdm); crc is NOT vectorizer-related (switching it off: 22.9 -> 22.6 ms, user cycles unchanged) but per-launch runtime overhead in Intel, confirmed by varying the launch count at fixed total work (`analysis/crc_launchsize.py`: Intel - PoCL = +7.7-7.9 us per launch, ~0 intercept, -0.6 ms at a single launch);
  PoCL's hmm gap is kernel-mode time (74% of cycles) with identical user-mode work; nw and srad tie. PoCL-side knobs (`analysis/pocl_knobs.py`) do NOT close the kmeans/lud gaps: kmeans is unchanged by every knob (hot loop not vectorized by PoCL; FP-reassociation allowed: no change);
  on lud the width knobs cut instructions 41-61% but time only 2.5-6.5%, and vectorizer-off is 7% faster, so lud's PoCL time is not limited by vector width (Intel's lud advantage is not explained by vector width). lud scratch-kernel test (`analysis/lud_chain.py`): splitting the accumulation into K independent sums does not help PoCL (K=2,4,8: no change; the serial-FMA-chain hypothesis is falsified), but removing the inner loop entirely (K=16) makes PoCL as fast as Intel (593 vs 602 ms; instructions 11x fewer, stores 55x fewer) -> loop structure, consistent with PoCL vectorizing only innermost loops (not proven).

## Limits (read before citing anything)
- One input per app, a few apps, one GPU; Nsight window = whole process per repetition (not the exact timed region); `ncu` counters come from one profiled launch per kernel.
- A finding is a measured signature, not a proven cause, until its intervention has been run. The "pinned staging" explanation for D1 is a hypothesis.
- Thresholds in the rules are judgement calls and have not been tuned or validated on held-out cases.
- NVIDIA's OpenCL driver is not visible to Nsight (no CUDA API/kernel tables), hence the interposer.
- Prior-art status: no complete equivalent found (about 75% confidence, not a systematic survey); closest work and unverified references are in the literature notes
  (kept outside the repository). Do not claim novelty before a second search pass.


---

## NOTICE (2026-11-25): earlier cross-model results are superseded
An equivalence audit found that the OpenMP and OpenCL versions did not time the same work, and FIX-OMP-14/15/16 and FIX-OCL-35/36 (see omp/FIXES.md, opencl/FIXES.md) changed the sources afterwards:
- every CPU speedup measured against omp-clang (final, Intel and thread-sweep tables) used a clang build without AVX/AVX-512/FMA: **do not quote them**;
- OpenMP-vs-OpenCL rows (GPU and CPU) for lud, spmv, bfs, crc, nw, srad, cfd, gem and tdm included extra OpenMP allocation/mapping inside the timed region, kmeans included a differently executed RMSE step, and CPU nqueens compared two different algorithms;
- the GPU clang-vs-nvhpc analysis of per-launch mapping overhead (crc and others) was measured with the per-repetition mapping in place and must be re-checked.
The PoCL-vs-Intel analysis ran one OpenCL program on both runtimes; its kmeans (RMSE inside the region) and nqueens (it used the FORCE_CPU kernel then) timings refer to the old sources.
No re-measurement has been done yet: it needs a new pre-registered run.
