#include "simple_motifs.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {
int HashToken(const std::string &s) {
  unsigned int h = 0;
  for (unsigned char c : s) h = (h * 131u + c) & 0x7fffffffu;
  return h;
}
bool LoadTdm(const std::string &events_path, const std::string &ivl_path,
             const std::string &episodes_path, std::vector<int> &events,
             std::vector<int> &patterns) {
  std::ifstream events_in(events_path), ivl_in(ivl_path), ep_in(episodes_path);
  if (!events_in || !ivl_in || !ep_in) return false;
  events.clear();
  patterns.clear();
  std::string line;
  while (std::getline(events_in, line)) {
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream ls(line);
    std::string token;
    while (ls >> token) {
      events.push_back(HashToken(token) & 1023);
    }
  }
  while (ivl_in && std::getline(ivl_in, line)) {
    std::istringstream ls(line);
    std::string token;
    while (ls >> token) patterns.push_back(HashToken(token) & 1023);
  }
  while (ep_in && std::getline(ep_in, line)) {
    std::istringstream ls(line);
    std::string token;
    while (ls >> token) patterns.push_back(HashToken(token) & 1023);
  }
  return !events.empty() && !patterns.empty();
}
void Run(std::vector<int> &input, std::vector<int> &patterns,
         std::vector<int> &states, std::vector<int> &accept, int n,
         int machines, int pattern_count) {
  int *pin = input.data();
  int *ppat = patterns.data();
  int *ps = states.data(), *pa = accept.data();
  const int chunk = (n + machines - 1) / machines;
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for map(to : pin[0 : n], ppat[0 : pattern_count]) map(from : ps[0 : machines], pa[0 : machines]) thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for schedule(static)
#endif
  for (int m = 0; m < machines; ++m) {
    int state = 0, count = 0;
    int begin = m * chunk, end = begin + chunk < n ? begin + chunk : n;
    for (int i = begin; i < end; ++i) {
      state = (state * 33 + pin[i]) & 1023;
      for (int p = 0; p < pattern_count; ++p) {
        if (((state ^ ppat[p]) & 63) == 7) ++count;
      }
    }
    ps[m] = state; pa[m] = count;
  }
}
}
int main(int argc, char **argv) {
  auto opts = dwarfs_omp::ParseOptions(argc, argv);
  int n = static_cast<int>(opts.size), machines = static_cast<int>(opts.secondary_size);
  std::vector<int> input(n), patterns(64);
  std::vector<int> states(machines), accept(machines);
  if (!opts.input_file.empty()) {
    if (!LoadTdm(opts.input_file, opts.input_file2, opts.input_file3, input,
                 patterns)) return 3;
    n = static_cast<int>(input.size());
  } else {
    for (int i = 0; i < n; ++i) input[i] = (i * 19 + 5) & 1023;
    for (std::size_t i = 0; i < patterns.size(); ++i) patterns[i] = (i * 17 + 3) & 1023;
  }
  int pattern_count = static_cast<int>(patterns.size());
  auto times = dwarfs_omp::AdaptiveRun(opts,
      [](){},
      [&]() { for (int it = 0; it < opts.iterations; ++it) Run(input, patterns, states, accept, n, machines, pattern_count); });
  bool valid = states[machines / 2] >= 0;
  dwarfs_omp::PrintResult("tdm", "finite-state-machine", input.size(),
                          static_cast<std::size_t>(machines), opts.iterations,
                          static_cast<int>(times.size()), times, valid);
  return valid ? 0 : 2;
}
