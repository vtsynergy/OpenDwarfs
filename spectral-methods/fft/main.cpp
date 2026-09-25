#include "simple_motifs.hpp"
#include <cmath>
#include <vector>

namespace {
using Value = float;
void PrecomputeTwiddles(std::vector<Value> &twiddle_r,
                        std::vector<Value> &twiddle_i, int n, int bins) {
  const float twopi = 6.283185307179586f;
  for (int k = 0; k < bins; ++k) {
    for (int t = 0; t < n; ++t) {
      float angle = -twopi * static_cast<float>(k * t) / static_cast<float>(n);
      twiddle_r[k * n + t] = std::cos(angle);
      twiddle_i[k * n + t] = std::sin(angle);
    }
  }
}

void Run(std::vector<Value> &real, std::vector<Value> &imag,
         std::vector<Value> &twiddle_r, std::vector<Value> &twiddle_i,
         std::vector<Value> &out_r, std::vector<Value> &out_i, int n,
         int bins) {
  Value *pr = real.data(), *pi = imag.data(), *por = out_r.data(), *poi = out_i.data();
  Value *pwr = twiddle_r.data(), *pwi = twiddle_i.data();
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for map(to : pr[0 : n], pi[0 : n], pwr[0 : n * bins], pwi[0 : n * bins]) map(from : por[0 : bins], poi[0 : bins]) thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for schedule(static)
#endif
  for (int k = 0; k < bins; ++k) {
    Value sr = 0.0f, si = 0.0f;
    for (int t = 0; t < n; ++t) {
      Value c = pwr[k * n + t], s = pwi[k * n + t];
      sr += pr[t] * c - pi[t] * s;
      si += pr[t] * s + pi[t] * c;
    }
    por[k] = sr; poi[k] = si;
  }
}
}
int main(int argc, char **argv) {
  auto opts = dwarfs_omp::ParseOptions(argc, argv);
  int n = static_cast<int>(opts.size), bins = static_cast<int>(std::min<std::size_t>(opts.secondary_size, opts.size));
  std::vector<Value> r(n), im(n), orr(bins), oi(bins);
  std::vector<Value> wr(static_cast<std::size_t>(n) * bins), wi(static_cast<std::size_t>(n) * bins);
  for (int i = 0; i < n; ++i) { r[i] = dwarfs_omp::InitValue(i); im[i] = dwarfs_omp::InitValue(i, 113); }
  PrecomputeTwiddles(wr, wi, n, bins);
  auto times = dwarfs_omp::AdaptiveRun(opts,
      [](){},
      [&]() { for (int it = 0; it < opts.iterations; ++it) Run(r, im, wr, wi, orr, oi, n, bins); });
  bool valid = std::isfinite(orr[bins / 2]) && std::isfinite(oi[bins / 2]);
  dwarfs_omp::PrintSimpleResult("fft", "spectral-methods", opts, times, valid);
  return valid ? 0 : 2;
}
