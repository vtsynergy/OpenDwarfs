# OpenACC (kept as imported from `openacc-beta`)

2014 PGI/CoreTSAR code: PolyBench-GPU ports, `helmholtz`, and `gem`. The sources are **unmodified**.
Build files added: `common.mk`, per-benchmark `Makefile`s, `benchmarks.mk`, and a `libctsar` rule
(`make libctsar`) that builds the CoreTSAR-lite runtime with `nvc++` and adds stand-ins for the removed
PGI OpenCL accessors.

```bash
source ~/.local/apps/openacc-<host>/env.sh        # loads the NVHPC module
make build                                        # GPU (cc70); prints BUILD OK/FAIL per benchmark
make -C polybench-gpu/dense-linear-algebra/2MM run GPU_INDEX=0
make build ACC_TARGET=multicore                   # x86 threads instead of GPU
```
Benchmarks that use OpenACC 1.0 syntax (`acc data region`, `acc for`) fail with a current `nvc`
until the sources are ported; see `build.log` in each directory. Status is recorded in the results
summary produced by `./run_all.sh`.
