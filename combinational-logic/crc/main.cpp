#include "simple_motifs.hpp"
#include <cstdint>
#include <vector>

namespace {
uint32_t Update(uint32_t crc, unsigned char byte) {
  crc ^= byte;
  for (int k = 0; k < 8; ++k) {
    crc = (crc >> 1) ^ (0xedb88320u & (0u - (crc & 1u)));
  }
  return crc;
}
void Run(std::vector<unsigned char> &data, std::vector<uint32_t> &partial,
         int n, int chunks) {
  unsigned char *pd = data.data();
  uint32_t *pp = partial.data();
  const int chunk_size = (n + chunks - 1) / chunks;
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for map(to : pd[0 : n]) map(from : pp[0 : chunks]) thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for schedule(static)
#endif
  for (int c = 0; c < chunks; ++c) {
    int begin = c * chunk_size;
    int end = begin + chunk_size < n ? begin + chunk_size : n;
    uint32_t crc = 0xffffffffu;
    for (int i = begin; i < end; ++i) crc = Update(crc, pd[i]);
    pp[c] = ~crc;
  }
}

void RunSerial(const std::vector<unsigned char> &data, std::vector<uint32_t> &partial, int n, int chunks) {
  const unsigned char *pd = data.data();
  uint32_t *pp = partial.data();
  const int chunk_size = (n + chunks - 1) / chunks;
  for (int c = 0; c < chunks; ++c) {
    int begin = c * chunk_size;
    int end = std::min(n, begin + chunk_size);
    uint32_t crc = 0xffffffffu;
    for (int i = begin; i < end; ++i) crc = Update(crc, pd[i]);
    pp[c] = ~crc;
  }
}

bool Verify(const std::vector<uint32_t> &partial_gpu, const std::vector<unsigned char> &data, int n, int chunks) {
  std::vector<uint32_t> partial_ref(chunks);
  RunSerial(data, partial_ref, n, chunks);
  for (std::size_t i = 0; i < partial_gpu.size(); ++i) {
    if (partial_gpu[i] != partial_ref[i]) {
      return false;
    }
  }
  return true;
}

}
int main(int argc, char **argv) {
  auto opts = dwarfs_omp::ParseOptions(argc, argv);
  int n = static_cast<int>(opts.size);
  int chunks = static_cast<int>(std::max<std::size_t>(1, opts.secondary_size));
  std::vector<unsigned char> data(n);
  std::vector<uint32_t> partial(chunks);
  for (int i = 0; i < n; ++i) data[i] = static_cast<unsigned char>((static_cast<unsigned int>(i) * 31u + 7u) & 255u);
  auto times = dwarfs_omp::AdaptiveRun(opts,
      [](){},
      [&]() { for (int it = 0; it < opts.iterations; ++it) Run(data, partial, n, chunks); });
  bool valid = true;
  if (opts.verify) {
    valid = Verify(partial, data, n, chunks);
  }
  dwarfs_omp::PrintSimpleResult("crc", "combinational-logic", opts, times, valid);
  return valid ? 0 : 2;
}
