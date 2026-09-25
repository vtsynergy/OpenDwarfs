#ifndef DWARFS_OMP_COMMON_SIMPLE_MOTIFS_HPP
#define DWARFS_OMP_COMMON_SIMPLE_MOTIFS_HPP

#include "dwarf_common.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <numeric>
#include <vector>

namespace dwarfs_omp {

inline float InitValue(int i, int mod = 127) {
  return static_cast<float>((static_cast<int64_t>(i) * 17 + 13) % mod) / static_cast<float>(mod);
}

inline void PrintSimpleResult(const char *benchmark, const char *motif,
                              const Options &opts,
                              const std::vector<double> &times, bool valid) {
  PrintResult(benchmark, motif, opts.size, opts.secondary_size, opts.iterations,
              static_cast<int>(times.size()), times, valid);
}

}  // namespace dwarfs_omp

#endif  // DWARFS_OMP_COMMON_SIMPLE_MOTIFS_HPP
