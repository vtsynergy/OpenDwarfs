#include "dwarf_common.hpp"

#include <omp.h>

#include <iostream>
#include <vector>

namespace {
using Value = float;

void Init(std::vector<Value> &a, std::vector<Value> &b, std::vector<Value> &c,
          int n) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      a[i * n + j] = static_cast<Value>((i * j + 1) % 97) / 97.0f;
      b[i * n + j] = static_cast<Value>((i + j + 3) % 89) / 89.0f;
      c[i * n + j] = static_cast<Value>((i - j + 11) % 83) / 83.0f;
    }
  }
}

void Run(std::vector<Value> &a, std::vector<Value> &b, std::vector<Value> &c,
         int n) {
  Value *pa = a.data();
  Value *pb = b.data();
  Value *pc = c.data();
  const Value alpha = 1.25f;
  const Value beta = 0.75f;
#if DWARF_USE_TARGET
#pragma omp target data map(to : pa[0 : n * n], pb[0 : n * n]) map(tofrom : pc[0 : n * n])
  {
#pragma omp target teams distribute parallel for collapse(2) thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        Value sum = 0.0f;
        for (int k = 0; k < n; ++k) {
          sum += pa[i * n + k] * pb[k * n + j];
        }
        pc[i * n + j] = alpha * sum + beta * pc[i * n + j];
      }
    }
#if DWARF_USE_TARGET
  }
#endif
}

bool Verify(const std::vector<Value> &a, const std::vector<Value> &b,
            const std::vector<Value> &got, const std::vector<Value> &initial_c,
            int n, int iterations) {
  const Value alpha = 1.25f;
  const Value beta = 0.75f;
  int failures = 0;
  for (int i = 0; i < n; i += std::max(1, n / 16)) {
    for (int j = 0; j < n; j += std::max(1, n / 16)) {
      Value sum = 0.0f;
      for (int k = 0; k < n; ++k) {
        sum += a[i * n + k] * b[k * n + j];
      }
      Value expected = initial_c[i * n + j];
      for (int it = 0; it < iterations; ++it) {
        expected = alpha * sum + beta * expected;
      }
      if (!dwarfs_omp::NearlyEqual(expected, got[i * n + j], 2.0e-3)) {
        ++failures;
      }
    }
  }
  return failures == 0;
}
}  // namespace

int main(int argc, char **argv) {
  const auto opts = dwarfs_omp::ParseOptions(argc, argv);
  const int n = static_cast<int>(opts.size);
  std::vector<Value> a(n * n), b(n * n), c(n * n), initial_c;
  bool valid = true;

  auto times = dwarfs_omp::AdaptiveRun(opts,
      [&]() { Init(a, b, c, n); initial_c = c; },
      [&]() { for (int it = 0; it < opts.iterations; ++it) Run(a, b, c, n); });

  if (opts.verify) {
    valid = Verify(a, b, c, initial_c, n, opts.iterations);
  }
  dwarfs_omp::PrintResult("gemm", "dense-linear-algebra", opts.size,
                          opts.secondary_size, opts.iterations,
                          static_cast<int>(times.size()), times, valid);
  return valid ? 0 : 2;
}
