# OpenDwarfs OpenMP Target Offload

OpenMP-target-offload implementations of 19 Berkeley/OpenDwarfs compute
motifs across 13 dwarf categories. The implementations use OpenMP target
offload and are intended to run on any device with a toolchain and runtime
that support the OpenMP features used here. The supplied CPU, NVIDIA, and
AMD build presets and hardware profiles are reference configurations.

This branch uses plain Makefiles in place of the original autoconf build.

## Layout

```text
common/
  benchmark.mk                 Shared build/run Makefile rules (all devices).
  benchmarks.mk                List of the 19 application directories.
  dwarf_common.hpp             CLI parsing, timing, validation, result-print helpers.
  run_benchmark.sh             Thin wrapper used by each dwarf's run.sh.
  simple_motifs.hpp
scripts/
  install-llvm-offload.sh      Build/install an OpenMP-target-offload Clang.
  detect_devices.py            Auto-detect cache sizes -> devices.json.
  compute_cache_tier_sizes.py  Generic 7-tier cache-bounding problem sizer.
  run_cache_tier_sweep.sh      Generic empirical tier-size sweep runner.
  devices.example.json         Example hardware profile (edit/copy this).
test/
  Makefile                     Smoke tests for all 19 apps plus seven real-input cases.
  test_scripts.py              Script regression tests (Python standard library).
  check_binaries.py            CLI and input checks for built binaries.
  <category>/<dwarf>/          Reference input datasets.

branch-and-bound/nqueens
combinational-logic/crc
dense-linear-algebra/{atax,cg,gemm,kmeans,lud}
dynamic-programming/{nw,swat}
finite-state-machine/tdm
graphical-models/hmm
graph-traversal/bfs
mapreduce/pagerank
n-body-methods/gem
sparse-linear-algebra/spmv
spectral-methods/fft
structured-grids/{fdtd2d,srad}
unstructured-grids/cfd
```

Each dwarf directory has:

- `main.cpp` -- the implementation (CPU path + `#pragma omp target` path).
- `Makefile` -- a few lines setting `BENCHMARK`/defaults and including
  `../../common/benchmark.mk`.
- `run.sh` -- a thin wrapper around `make run DEVICE=<device>`.

## Getting an OpenMP-target-offload Clang

CPU builds need Clang with the OpenMP runtime; GPU builds also need the
offload runtime and the appropriate CUDA or ROCm toolkit. If this machine
does not already have an offload compiler, `scripts/install-llvm-offload.sh`
builds and installs one from upstream LLVM sources -- it auto-detects
whatever GPU vendor toolkit(s) are available (NVIDIA/CUDA, AMD/ROCm,
neither for a host-only build) and only builds support for what it finds,
with environment overrides for local toolkit paths and build settings:

```bash
./scripts/install-llvm-offload.sh

# Force a specific vendor, or point at non-default toolkit locations:
ENABLE_AMDGPU=0 CUDA_ROOT=/opt/cuda-12.9 ./scripts/install-llvm-offload.sh

# Load vendor toolkits via environment-modules first, if this machine uses it:
CUDA_MODULE=cuda/12.6.3 ROCM_MODULE=rocm/6.4.4 ./scripts/install-llvm-offload.sh

# Also build+run a quick multi-device OpenMP probe after installing:
RUN_SMOKE_TEST=1 ./scripts/install-llvm-offload.sh
```

It installs to `~/.local/apps/llvm-offload-<hostname>` by default (override
with `INSTALL_DIR=...`); see the script's own header comment for the full
list of environment-variable overrides.

## Build and run

Select the build target with `DEVICE=cpu`, `DEVICE=a100`, or another
supported device name.

The build uses GNU Make and defaults to `clang++` on `PATH`. Set `CXX` to
use an existing installation. Tests and sizing scripts also need Python 3;
the sweep runner needs Bash 4 or newer and GNU `timeout`.

```bash
export CXX=/path/to/llvm-offload/bin/clang++
# If CUDA is outside Clang's default search paths:
export CUDA_PATH=/path/to/cuda
# Require GPU runs to offload instead of falling back to the host:
export OMP_TARGET_OFFLOAD=MANDATORY
```

From the repository root, across every dwarf:

```bash
make build DEVICE=cpu
make build DEVICE=a100
make run   DEVICE=a100
make all   DEVICES="cpu a100"   # build every dwarf for each listed device
make test  DEVICE=a100         # synthetic and real-input tests
make clean
```

Within an application directory, use `build`, `run`, `all`, or `clean`:

```bash
make build DEVICE=a100
make run   DEVICE=a100

# or the equivalent script:
./run.sh a100
```

Included build presets (see `common/benchmark.mk`): `cpu`, `h100`, `h200`, `a100`,
`mi210`, `mi300a`, `mi100`.

When changing compilers or compiler flags, use `make clean` or a fresh
`BUILD_DIR=build-custom` to avoid reusing binaries built with old settings.

Override problem shape and measurement settings:

```bash
make -C dense-linear-algebra/gemm run DEVICE=a100 \
  RUN_SIZE=2048 RUN_ITERATIONS=5 RUN_REPEATS=3 RUN_WARMUP=1
```

The root `make run` applies overrides to every app; problem sizes have
different meanings for different apps. `RUN_SECONDARY_SIZE` sets the
secondary dimension. Invalid options and numeric values return an error.
Binaries accept the corresponding CLI flags directly:

```bash
./dense-linear-algebra/gemm/build/gemm_a100 --size 2048 --repeats 3 --warmup 1
```

Each run prints a machine-readable line:

```text
DWARF_RESULT,benchmark=gemm,motif=dense-linear-algebra,target=a100,...
```

### Adding a new hardware target

For a device outside the supplied presets, configure its compiler, offload
target, architecture flags, and runtime. Add a name to `KNOWN_DEVICES` and
extend the device rules in `common/benchmark.mk` to use those settings. To use a different architecture with an existing
device name, override its architecture variable:

```bash
make build DEVICE=a100 A100_ARCH=sm_86 CXX=/path/to/clang++
```

When adding a device name, also update the sweep runner's target validation
and GPU environment selection in `scripts/run_cache_tier_sweep.sh`.
An architecture override keeps the existing device label in results.

For analytical sizing, add the device's cache capacities to your hardware
profile; runtime tiering uses measurements from the device itself. Neither
method is limited to the hardware in `scripts/devices.example.json`.

## Tests and reference inputs

Run `make test DEVICE=cpu` or `make test DEVICE=a100` from the repository
root. The suite runs script regression tests, small synthetic cases for all
19 apps, seven real-input cases, and CLI/input regression checks. Each
application must exit successfully and report `valid=1`. Validation varies by app: some compare against a
serial reference, while others only check finiteness or bounds.
Use `OMP_TARGET_OFFLOAD=MANDATORY` to require GPU execution.

The real-input datasets live under `test/<category>/<app>/`. Binaries accept
these flags, with file paths relative to the current directory:

| App | Input flags |
| --- | --- |
| bfs | `--input <graph.txt>` |
| lud | `--input <matrix.dat>` |
| kmeans | `--input <points-file>` |
| swat | `--input <query> --input2 <database>` |
| tdm | `--input <events.csv> --input2 <ivl.txt> --input3 <episodes.txt>` |
| gem | `--input <molecule.pqr> --input2 <molecule.vert>` |
| cfd | `--input <mesh.domn.*>` |

For example, from the root:

```bash
./graph-traversal/bfs/build/bfs_a100 --input test/graph-traversal/bfs/graph4096.txt
```

`test/mapreduce/fsa-blastcl` contains a sample that is not wired into an app.
The upstream CRC and SPMV reference datasets are not included.

## Problem sizing: inversion and runtime tiering

We use two methods to choose problem sizes. For the streaming apps, we
invert an analytical working-set model to obtain sizes at cache-based byte
targets. For the remaining apps, we choose sizes from measured runtimes.
Both methods produce entries in the same tier manifest, which the sweep
runner uses to build and run the apps.

### Analytical inversion

`atax`, `crc`, `fdtd2d`, `fft`, `kmeans`, `pagerank`, `spmv`, and `srad`
use analytical sizing. Given a working-set model `W(N)` and a tier's byte
target `B`, the sizer computes `N = W⁻¹(B)` and rounds to an integer problem
size. Secondary dimensions are fixed by the model or set equal to `N`.
For example, the SpMV model with 16 entries per row is `W(N) = 140N + 4`,
so its size is `floor((B - 4) / 140)`, with a minimum of one row.

The seven byte targets come from a hardware profile:

| Tier | Byte target |
| --- | --- |
| `l1lo`, `l1hi` | Minimum and maximum L1 capacity across the listed devices |
| `l2lo`, `l2hi` | Minimum and maximum L2 capacity |
| `llclo`, `llchi` | Minimum and maximum LLC capacity |
| `stream` | 1 GiB |

Cache sizes are aggregate device capacities. They provide sizing reference
points, not a guarantee of cache residency. The profile needs at least one
non-null value for each cache level. Tier targets can coincide; `stream`
is fixed rather than derived from the profile.

From the repository root:

```bash
./scripts/detect_devices.py --out scripts/devices.json

# Size all eight streaming apps.
./scripts/compute_cache_tier_sizes.py --devices scripts/devices.json

# Size only ATAX and SpMV.
./scripts/compute_cache_tier_sizes.py --devices scripts/devices.json \
  --dwarfs atax spmv
```

The detector reads CPU cache totals from `lscpu` and AMD cache information
from `rocminfo`. For NVIDIA GPUs it identifies the device with `nvidia-smi`
and looks up cache sizes in a compute-capability table. Review those values
for your GPU configuration. To enter hardware details by hand, copy
`scripts/devices.example.json` to `scripts/devices.json` and edit it.
Capacities are in bytes; use `null` for a cache level that does not apply.

The full output is `data/generated/cache_tier_sizes_7tier.json`. A subset
gets a filename such as `cache_tier_sizes_7tier__atax_spmv.json` in the same
directory. Use `--out` to choose another path. Repeating a command overwrites
its output.

### Runtime tiering

`gemm`, `cg`, `lud`, `gem`, `nqueens`, `nw`, `swat`, `hmm`, `bfs`, `tdm`,
and `cfd` use runtime tiering. Run candidate problem sizes on the intended
devices, compare their measured runtimes, and adjust the sizes to select
the tiers. The sweep script handles the measurements; choosing and refining
the candidates is a manual step in the current tooling.

Create `data/generated/` and save candidates in
`data/generated/empirical_seed.json`. For example, this starts a GEMM trial
at size 64; it is not a calibrated tier size:

```json
{"dwarfs": {"gemm": {"l1lo": {"size": 64, "secondary_size": 64}}}}
```

Merge the candidates with the analytical sizes, then measure the trial:

```bash
./scripts/compute_cache_tier_sizes.py --devices scripts/devices.json \
  --seed data/generated/empirical_seed.json \
  --out data/generated/cache_tier_sizes_seeded.json

DWARF_TIERS_MANIFEST=data/generated/cache_tier_sizes_seeded.json \
DWARF_TARGETS="cpu a100" DWARF_FILTER="gemm" DWARF_TIERS="l1lo" \
  ./scripts/run_cache_tier_sweep.sh
```

The runner uses the exported `CXX` and toolkit settings from the build
instructions. It writes timings and validation results to
`results/<run-id>/summary.csv`, with a log for each run. Read the timings,
edit the candidate sizes, regenerate the manifest, and repeat. Add entries
for the other tiers as you select them. The manifest stores `size`,
`secondary_size`, and an optional `iterations` value for each app and tier;
the runner does not assign runtime thresholds or search for sizes itself.

Seed entries replace the entire tier map for each named app. Leave out
`--dwarfs` when merging runtime-tiered apps: that option restricts both
sizing and seed inclusion to the selected streaming apps. Use a separate
`--out` for seed trials, since default filenames record only the streaming
app selection.

### Running a sweep

The same runner accepts analytical sizes, runtime-tiered sizes, or a
manifest containing both. Remove `DWARF_FILTER` and `DWARF_TIERS` from the
example above to run all available entries. Missing tier entries are
skipped. Failed runs return a nonzero script status; a selection with no
runnable entries also fails.

| Setting | Default |
| --- | --- |
| `DWARF_TARGETS` | `cpu` |
| `DWARF_TIERS` | All seven tiers |
| `DWARF_FILTER` | All apps with sizes in the manifest |
| `DWARF_WARMUP` | `1` |
| `DWARF_REPEATS` | `3` |
| `DWARF_TIMEOUT_S` | `600` per build/run command |
| `DWARF_GPU_DEV` | GPU index `0` |
| `DWARF_BUILD_DIR` | `build` in each app directory |

## Generated files

Use `build/` or a `build-*` directory in each app for binaries. Build
outputs, Python caches, temporary files, `data/`, `results/`, and the local
`scripts/devices.json` profile are gitignored. Keep scratch files in `tmp/`
or `temp/`, and save generated manifests and seeds under `data/generated/`.
Reference inputs under `test/` and `scripts/devices.example.json` are kept
in the repository.

`make clean` removes each app's selected `BUILD_DIR` (default `build`);
for a custom directory use `make clean BUILD_DIR=build-custom`. It does not
remove generated manifests, result directories, or Python caches. Custom
output paths outside the ignored locations need their own ignore rules.

## Notes

- `RUN_REPEATS` controls the exact number of measured runs; `RUN_WARMUP`
  controls the number of untimed warmups.
- GPU validation can require actual offload with `OMP_TARGET_OFFLOAD=MANDATORY`.
- Script regression tests can also run alone with
  `python3 -m unittest discover -s test -v`.
- `DWARF_THREAD_LIMIT` (default 256) controls OpenMP target launch
  configuration; pass it as a make/env override to tune per device.

## Citation

If you use this work, please cite:

```bibtex
@inproceedings{chaudhury2026opendwarfs2,
  title     = {OpenDwarfs 2.0: Evolving the Berkeley Dwarfs in OpenCL and OpenMP},
  author    = {Chaudhury, Nabayan and Feng, Wu-chun},
  booktitle = {Proceedings of the IEEE High Performance Extreme Computing Conference (HPEC)},
  year      = {2026},
  note      = {Best Student Paper Award. TBD: volume/pages/DOI once the
               proceedings are indexed on IEEE Xplore.}
}
```
