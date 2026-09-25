#include "simple_motifs.hpp"
#include <cstdint>
#include <vector>

namespace {
using Count = unsigned long long;
Count Solve(int n, int row, uint64_t cols, uint64_t diag1, uint64_t diag2) {
  if (row == n) return 1;
  Count total = 0;
  uint64_t avail = ((1ull << n) - 1ull) & ~(cols | diag1 | diag2);
  while (avail) {
    uint64_t bit = avail & -avail;
    avail -= bit;
    total += Solve(n, row + 1, cols | bit, (diag1 | bit) << 1, (diag2 | bit) >> 1);
  }
  return total;
}
void Run(std::vector<Count> &partial, int n, int tasks) {
  Count *pp = partial.data();
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for map(from : pp[0 : tasks]) firstprivate(n, tasks) thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for schedule(dynamic)
#endif
  for (int t = 0; t < tasks; ++t) {
    int col = t % n;
    uint64_t bit = 1ull << col;
    pp[t] = Solve(n, 1, bit, bit << 1, bit >> 1);
  }
}
}
int main(int argc, char **argv) {
  auto opts = dwarfs_omp::ParseOptions(argc, argv);
  int n = static_cast<int>(std::min<std::size_t>(opts.size, 15));
  int tasks = static_cast<int>(opts.secondary_size);
  std::vector<Count> partial(tasks);
  auto times = dwarfs_omp::AdaptiveRun(opts,
      [](){},
      [&]() { for (int it = 0; it < opts.iterations; ++it) Run(partial, n, tasks); });
  Count total = 0; for (Count v : partial) total += v;
  bool valid = total > 0;
  dwarfs_omp::PrintSimpleResult("nqueens", "branch-and-bound", opts, times, valid);
  return valid ? 0 : 2;
}
