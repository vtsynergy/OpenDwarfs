#include "simple_motifs.hpp"
#include <vector>

namespace {
using Value = float;
void Run(std::vector<Value> &img, std::vector<Value> &next, int n, int iterations) {
  Value *pi = img.data(), *pn = next.data();
#if DWARF_USE_TARGET
#pragma omp target data map(tofrom : pi[0 : n * n]) map(alloc : pn[0 : n * n])
  {
#endif
    for (int it = 0; it < iterations; ++it) {
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for collapse(2) thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for collapse(2) schedule(static)
#endif
      for (int i = 1; i < n - 1; ++i) {
        for (int j = 1; j < n - 1; ++j) {
          int p = i * n + j;
          Value lap = pi[p - 1] + pi[p + 1] + pi[p - n] + pi[p + n] - 4.0f * pi[p];
          pn[p] = pi[p] + 0.25f * lap;
        }
      }
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for collapse(2) thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for collapse(2) schedule(static)
#endif
      for (int i = 1; i < n - 1; ++i) for (int j = 1; j < n - 1; ++j) pi[i * n + j] = pn[i * n + j];
    }
#if DWARF_USE_TARGET
  }
#endif
}
}
int main(int argc, char **argv) {
  auto opts = dwarfs_omp::ParseOptions(argc, argv);
  int n = static_cast<int>(opts.size);
  std::vector<Value> img(n * n), next(n * n);
  for (int i = 0; i < n * n; ++i) img[i] = dwarfs_omp::InitValue(i);
  auto times = dwarfs_omp::AdaptiveRun(opts,
      [](){},
      [&]() { Run(img, next, n, opts.iterations); });
  bool valid = std::isfinite(img[(n / 2) * n + n / 2]);
  dwarfs_omp::PrintSimpleResult("srad", "structured-grids", opts, times, valid);
  return valid ? 0 : 2;
}
