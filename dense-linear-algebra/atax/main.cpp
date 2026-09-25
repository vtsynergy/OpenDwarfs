#include "dwarf_common.hpp"

#include <vector>

namespace {
using Value = float;

void Init(std::vector<Value> &a, std::vector<Value> &x, int rows, int cols) {
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      a[i * cols + j] = static_cast<Value>((i * j + 7) % 101) / 101.0f;
    }
  }
  for (int j = 0; j < cols; ++j) {
    x[j] = static_cast<Value>((j + 5) % 67) / 67.0f;
  }
}

void Run(std::vector<Value> &a, std::vector<Value> &x, std::vector<Value> &y,
         std::vector<Value> &tmp, int rows, int cols) {
  Value *pa = a.data();
  Value *px = x.data();
  Value *py = y.data();
  Value *pt = tmp.data();
#if DWARF_USE_TARGET
#pragma omp target data map(to : pa[0 : rows * cols], px[0 : cols]) map(alloc : pt[0 : rows]) map(from : py[0 : cols])
  {
#pragma omp target teams distribute parallel for thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < rows; ++i) {
      Value sum = 0.0f;
      for (int j = 0; j < cols; ++j) {
        sum += pa[i * cols + j] * px[j];
      }
      pt[i] = sum;
    }
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for schedule(static)
#endif
    for (int j = 0; j < cols; ++j) {
      Value sum = 0.0f;
      for (int i = 0; i < rows; ++i) {
        sum += pa[i * cols + j] * pt[i];
      }
      py[j] = sum;
    }
#if DWARF_USE_TARGET
  }
#endif
}

bool Verify(const std::vector<Value> &a, const std::vector<Value> &x,
            const std::vector<Value> &y, int rows, int cols) {
  int failures = 0;
  for (int j = 0; j < cols; j += std::max(1, cols / 32)) {
    Value expected = 0.0f;
    for (int i = 0; i < rows; ++i) {
      Value tmp = 0.0f;
      for (int k = 0; k < cols; ++k) {
        tmp += a[i * cols + k] * x[k];
      }
      expected += a[i * cols + j] * tmp;
    }
    if (!dwarfs_omp::NearlyEqual(expected, y[j], 3.0e-3)) {
      ++failures;
    }
  }
  return failures == 0;
}
}  // namespace

int main(int argc, char **argv) {
  const auto opts = dwarfs_omp::ParseOptions(argc, argv);
  const int rows = static_cast<int>(opts.size);
  const int cols = static_cast<int>(opts.secondary_size);
  std::vector<Value> a(rows * cols), x(cols), y(cols), tmp(rows);
  bool valid = true;

  Init(a, x, rows, cols);
  auto times = dwarfs_omp::AdaptiveRun(opts,
      [&]() { std::fill(y.begin(), y.end(), 0.0f); },
      [&]() { for (int it = 0; it < opts.iterations; ++it) Run(a, x, y, tmp, rows, cols); });
  if (opts.verify) {
    valid = Verify(a, x, y, rows, cols);
  }
  dwarfs_omp::PrintResult("atax", "dense-linear-algebra", opts.size,
                          opts.secondary_size, opts.iterations,
                          static_cast<int>(times.size()), times, valid);
  return valid ? 0 : 2;
}
