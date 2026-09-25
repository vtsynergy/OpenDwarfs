#include "simple_motifs.hpp"
#include <fstream>
#include <vector>

namespace {
using Value = float;
bool LoadMatrix(const std::string &path, std::vector<Value> &a, int &n) {
  std::ifstream in(path);
  if (!in) return false;
  in >> n;
  a.resize(static_cast<std::size_t>(n) * n);
  for (int i = 0; i < n * n; ++i) in >> a[i];
  return static_cast<bool>(in);
}
void Run(std::vector<Value> &a, int n) {
  Value *pa = a.data();
  for (int k = 0; k < n; ++k) {
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for map(tofrom : pa[0 : n * n]) firstprivate(k, n) thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for schedule(static)
#endif
    for (int i = k + 1; i < n; ++i) pa[i * n + k] = pa[i * n + k] / pa[k * n + k];
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for collapse(2) map(tofrom : pa[0 : n * n]) firstprivate(k, n) thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int i = k + 1; i < n; ++i) {
      // Use the same rounding on host and device; implicit GPU contraction
      // can otherwise diverge from the serial reference for small pivots.
      for (int j = k + 1; j < n; ++j)
        pa[i * n + j] = std::fma(-pa[i * n + k], pa[k * n + j], pa[i * n + j]);
    }
  }
}
void RunSerial(std::vector<Value> &a, int n) {
  Value *pa = a.data();
  for (int k = 0; k < n; ++k) {
    for (int i = k + 1; i < n; ++i) pa[i * n + k] = pa[i * n + k] / pa[k * n + k];
    for (int i = k + 1; i < n; ++i) {
      for (int j = k + 1; j < n; ++j)
        pa[i * n + j] = std::fma(-pa[i * n + k], pa[k * n + j], pa[i * n + j]);
    }
  }
}

bool Verify(const std::vector<Value> &a_gpu, std::vector<Value> a_ref, int n, int iters) {
  for (int it = 0; it < iters; ++it) {
    RunSerial(a_ref, n);
  }
  for (std::size_t i = 0; i < a_gpu.size(); ++i) {
    if (!dwarfs_omp::NearlyEqual(a_gpu[i], a_ref[i], 1e-2)) {
      return false;
    }
  }
  return true;
}

}
int main(int argc, char **argv) {
  auto opts = dwarfs_omp::ParseOptions(argc, argv);
  int n = static_cast<int>(opts.size);
  std::vector<Value> a(n * n);
  if (!opts.input_file.empty()) {
    if (!LoadMatrix(opts.input_file, a, n)) return 3;
  } else {
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) a[i * n + j] = (i == j ? n : 0.0f) + dwarfs_omp::InitValue(i * n + j);
  }
  auto a_orig = a;
  std::vector<Value> b;
  auto times = dwarfs_omp::AdaptiveRun(opts,
      [&]() { b = a_orig; },
      [&]() { for (int it = 0; it < opts.iterations; ++it) Run(b, n); });
  a.swap(b);

  bool valid = true;
  if (opts.verify) {
    valid = Verify(a, a_orig, n, opts.iterations);
  }
  dwarfs_omp::PrintResult("lud", "dense-linear-algebra",
                          static_cast<std::size_t>(n),
                          static_cast<std::size_t>(n), opts.iterations,
                          static_cast<int>(times.size()), times, valid);
  return valid ? 0 : 2;
}
