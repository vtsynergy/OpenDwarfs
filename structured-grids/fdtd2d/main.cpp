#include "dwarf_common.hpp"

#include <vector>

namespace {
using Value = float;

void Init(std::vector<Value> &ex, std::vector<Value> &ey, std::vector<Value> &hz,
          int n) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      ex[i * n + j] = static_cast<Value>((i * (j + 1)) % 113) / 113.0f;
      ey[i * n + j] = static_cast<Value>(((i + 1) * j) % 109) / 109.0f;
      hz[i * n + j] = static_cast<Value>((i + j + 1) % 107) / 107.0f;
    }
  }
}

void Run(std::vector<Value> &ex, std::vector<Value> &ey, std::vector<Value> &hz,
         int n, int iterations) {
  Value *pex = ex.data();
  Value *pey = ey.data();
  Value *phz = hz.data();
#if DWARF_USE_TARGET
#pragma omp target data map(tofrom : pex[0 : n * n], pey[0 : n * n], phz[0 : n * n])
  {
#endif
    for (int t = 0; t < iterations; ++t) {
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for collapse(2) thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for collapse(2) schedule(static)
#endif
      for (int i = 0; i < n; ++i) {
        for (int j = 1; j < n; ++j) {
          pex[i * n + j] -= 0.5f * (phz[i * n + j] - phz[i * n + j - 1]);
        }
      }
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for collapse(2) thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for collapse(2) schedule(static)
#endif
      for (int i = 1; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
          pey[i * n + j] -= 0.5f * (phz[i * n + j] - phz[(i - 1) * n + j]);
        }
      }
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for collapse(2) thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for collapse(2) schedule(static)
#endif
      for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - 1; ++j) {
          phz[i * n + j] -= 0.7f * (pex[i * n + j + 1] - pex[i * n + j] +
                                    pey[(i + 1) * n + j] - pey[i * n + j]);
        }
      }
    }
#if DWARF_USE_TARGET
  }
#endif
}

bool VerifyFinite(const std::vector<Value> &ex, const std::vector<Value> &ey,
                  const std::vector<Value> &hz) {
  for (Value v : ex) {
    if (!std::isfinite(v)) return false;
  }
  for (Value v : ey) {
    if (!std::isfinite(v)) return false;
  }
  for (Value v : hz) {
    if (!std::isfinite(v)) return false;
  }
  return true;
}
}  // namespace

int main(int argc, char **argv) {
  const auto opts = dwarfs_omp::ParseOptions(argc, argv);
  const int n = static_cast<int>(opts.size);
  std::vector<Value> ex(n * n), ey(n * n), hz(n * n);
  bool valid = true;

  auto times = dwarfs_omp::AdaptiveRun(opts,
      [&]() { Init(ex, ey, hz, n); },
      [&]() { Run(ex, ey, hz, n, opts.iterations); });
  if (opts.verify) {
    valid = VerifyFinite(ex, ey, hz);
  }
  dwarfs_omp::PrintResult("fdtd2d", "structured-grids", opts.size,
                          opts.secondary_size, opts.iterations,
                          static_cast<int>(times.size()), times, valid);
  return valid ? 0 : 2;
}
