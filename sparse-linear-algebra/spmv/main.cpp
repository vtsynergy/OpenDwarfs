#include "dwarf_common.hpp"

#include <vector>

namespace {
using Value = float;

void Init(std::vector<int> &row_ptr, std::vector<int> &col_idx,
          std::vector<Value> &vals, std::vector<Value> &x, int rows,
          int nnz_per_row) {
  row_ptr.resize(rows + 1);
  col_idx.resize(static_cast<std::size_t>(rows) * nnz_per_row);
  vals.resize(static_cast<std::size_t>(rows) * nnz_per_row);
  x.resize(rows);
  for (int i = 0; i <= rows; ++i) {
    // cast to size_t before multiply to avoid 32-bit overflow for large rows
    row_ptr[i] = static_cast<int>(static_cast<std::size_t>(i) * nnz_per_row);
  }
  for (int i = 0; i < rows; ++i) {
    x[i] = static_cast<Value>((i + 1) % 71) / 71.0f;
    for (int k = 0; k < nnz_per_row; ++k) {
      const std::size_t idx = static_cast<std::size_t>(i) * nnz_per_row + k;
      col_idx[idx] = (i + k * 17) % rows;
      vals[idx] = static_cast<Value>((idx + 3) % 97) / 97.0f;
    }
  }
}

void Run(std::vector<int> &row_ptr, std::vector<int> &col_idx,
         std::vector<Value> &vals, std::vector<Value> &x,
         std::vector<Value> &y, int rows, std::size_t nnz) {
  int *pr = row_ptr.data();
  int *pc = col_idx.data();
  Value *pv = vals.data();
  Value *px = x.data();
  Value *py = y.data();
#if DWARF_USE_TARGET
#pragma omp target data map(to : pr[0 : rows + 1], pc[0 : nnz], pv[0 : nnz], px[0 : rows]) map(from : py[0 : rows])
  {
#pragma omp target teams distribute parallel for thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < rows; ++i) {
      Value sum = 0.0f;
      for (int jj = pr[i]; jj < pr[i + 1]; ++jj) {
        sum += pv[jj] * px[pc[jj]];
      }
      py[i] = sum;
    }
#if DWARF_USE_TARGET
  }
#endif
}

bool Verify(const std::vector<int> &row_ptr, const std::vector<int> &col_idx,
            const std::vector<Value> &vals, const std::vector<Value> &x,
            const std::vector<Value> &y, int rows) {
  int failures = 0;
  for (int i = 0; i < rows; i += std::max(1, rows / 64)) {
    Value expected = 0.0f;
    for (int jj = row_ptr[i]; jj < row_ptr[i + 1]; ++jj) {
      expected += vals[jj] * x[col_idx[jj]];
    }
    if (!dwarfs_omp::NearlyEqual(expected, y[i], 1.0e-4)) {
      ++failures;
    }
  }
  return failures == 0;
}
}  // namespace

int main(int argc, char **argv) {
  const auto opts = dwarfs_omp::ParseOptions(argc, argv);
  const int rows = static_cast<int>(opts.size);
  const int nnz_per_row = static_cast<int>(opts.secondary_size);
  const std::size_t nnz = static_cast<std::size_t>(rows) * nnz_per_row;
  std::vector<int> row_ptr, col_idx;
  std::vector<Value> vals, x, y(rows);
  bool valid = true;

  Init(row_ptr, col_idx, vals, x, rows, nnz_per_row);
  auto times = dwarfs_omp::AdaptiveRun(opts,
      [](){},
      [&]() { for (int it = 0; it < opts.iterations; ++it) Run(row_ptr, col_idx, vals, x, y, rows, nnz); });
  if (opts.verify) {
    valid = Verify(row_ptr, col_idx, vals, x, y, rows);
  }
  dwarfs_omp::PrintResult("spmv", "sparse-linear-algebra", opts.size,
                          opts.secondary_size, opts.iterations,
                          static_cast<int>(times.size()), times, valid);
  return valid ? 0 : 2;
}
