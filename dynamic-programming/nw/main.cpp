#include "simple_motifs.hpp"
#include <algorithm>
#include <vector>

namespace {
void Run(std::vector<int> &score, int n) {
  int *ps = score.data();
  for (int wave = 2; wave <= 2 * (n - 1); ++wave) {
    int lo = std::max(1, wave - (n - 1));
    int hi = std::min(n - 1, wave - 1);
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for map(tofrom : ps[0 : n * n]) firstprivate(n, wave, lo, hi) thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for schedule(static)
#endif
    for (int i = lo; i <= hi; ++i) {
      int j = wave - i;
      int match = ps[(i - 1) * n + j - 1] + ((i + j) % 7 == 0 ? 2 : -1);
      int del = ps[(i - 1) * n + j] - 1;
      int ins = ps[i * n + j - 1] - 1;
      ps[i * n + j] = std::max(match, std::max(del, ins));
    }
  }
}

void RunSerial(std::vector<int> &score, int n) {
  int *ps = score.data();
  for (int i = 1; i < n; ++i) {
    for (int j = 1; j < n; ++j) {
      int match = ps[(i - 1) * n + j - 1] + ((i + j) % 7 == 0 ? 2 : -1);
      int del = ps[(i - 1) * n + j] - 1;
      int ins = ps[i * n + j - 1] - 1;
      ps[i * n + j] = std::max(match, std::max(del, ins));
    }
  }
}

bool Verify(const std::vector<int> &score_gpu, int n) {
  std::vector<int> score_ref(n * n, 0);
  for (int i = 0; i < n; ++i) { score_ref[i * n] = -i; score_ref[i] = -i; }
  RunSerial(score_ref, n);
  for (std::size_t i = 0; i < score_gpu.size(); ++i) {
    if (score_gpu[i] != score_ref[i]) {
      return false;
    }
  }
  return true;
}

}
int main(int argc, char **argv) {
  auto opts = dwarfs_omp::ParseOptions(argc, argv);
  int n = static_cast<int>(opts.size);
  std::vector<int> score(n * n, 0);
  for (int i = 0; i < n; ++i) { score[i * n] = -i; score[i] = -i; }
  std::vector<int> s;
  auto times = dwarfs_omp::AdaptiveRun(opts,
      [&]() { s = score; },
      [&]() { for (int it = 0; it < opts.iterations; ++it) { s = score; Run(s, n); } });
  score.swap(s);

  bool valid = true;
  if (opts.verify) {
    valid = Verify(score, n);
  }
  dwarfs_omp::PrintSimpleResult("nw", "dynamic-programming", opts, times, valid);
  return valid ? 0 : 2;
}
