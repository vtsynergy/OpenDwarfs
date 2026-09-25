#include "simple_motifs.hpp"
#include <fstream>
#include <sstream>
#include <vector>

namespace {
using Value = float;
bool LoadAsciiPoints(const std::string &path, std::vector<Value> &x,
                     std::vector<Value> &y, int &n) {
  std::ifstream in(path);
  if (!in) return false;
  std::string line;
  std::vector<Value> xs;
  std::vector<Value> ys;
  while (std::getline(in, line)) {
    std::istringstream ls(line);
    std::vector<Value> vals;
    std::string tok;
    while (ls >> tok) {
      try {
        vals.push_back(std::stof(tok));
      } catch (...) {
      }
    }
    if (vals.size() >= 3) {
      xs.push_back(vals[1]);
      ys.push_back(vals[2]);
    } else if (vals.size() == 2) {
      xs.push_back(vals[0]);
      ys.push_back(vals[1]);
    }
  }
  if (xs.empty()) return false;
  x.swap(xs);
  y.swap(ys);
  n = static_cast<int>(x.size());
  return true;
}
void Run(std::vector<Value> &x, std::vector<Value> &y, std::vector<Value> &cx,
         std::vector<Value> &cy, std::vector<int> &labels, int n, int k,
         int iterations) {
  Value *px = x.data(), *py = y.data(), *pcx = cx.data(), *pcy = cy.data();
  int *pl = labels.data();
  for (int it = 0; it < iterations; ++it) {
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for map(to : px[0 : n], py[0 : n], pcx[0 : k], pcy[0 : k]) map(from : pl[0 : n]) thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < n; ++i) {
      int best = 0; Value best_dist = 1.0e30f;
      for (int c = 0; c < k; ++c) {
        Value dx = px[i] - pcx[c], dy = py[i] - pcy[c], d = dx * dx + dy * dy;
        if (d < best_dist) { best_dist = d; best = c; }
      }
      pl[i] = best;
    }
  }
}
}
int main(int argc, char **argv) {
  auto opts = dwarfs_omp::ParseOptions(argc, argv);
  int n = static_cast<int>(opts.size), k = static_cast<int>(opts.secondary_size);
  std::vector<Value> x(n), y(n), cx(k), cy(k);
  std::vector<int> labels(n);
  if (!opts.input_file.empty()) {
    if (!LoadAsciiPoints(opts.input_file, x, y, n)) return 3;
    labels.resize(n);
  } else {
    for (int i = 0; i < n; ++i) { x[i] = dwarfs_omp::InitValue(i); y[i] = dwarfs_omp::InitValue(i, 113); }
  }
  for (int c = 0; c < k; ++c) { cx[c] = dwarfs_omp::InitValue(c, 97); cy[c] = dwarfs_omp::InitValue(c, 89); }
  auto times = dwarfs_omp::AdaptiveRun(opts,
      [](){},
      [&]() { Run(x, y, cx, cy, labels, n, k, opts.iterations); });
  bool valid = labels[n / 2] >= 0 && labels[n / 2] < k;
  dwarfs_omp::PrintResult("kmeans", "dense-linear-algebra",
                          static_cast<std::size_t>(n),
                          static_cast<std::size_t>(k), opts.iterations,
                          static_cast<int>(times.size()), times, valid);
  return valid ? 0 : 2;
}
