#include "simple_motifs.hpp"
#include <fstream>
#include <vector>

namespace {
using Value = float;
bool LoadDomn(const std::string &path, std::vector<int> &nbr,
              std::vector<Value> &q, int &cells, int &degree) {
  std::ifstream in(path);
  if (!in) return false;
  in >> cells;
  degree = 4;
  nbr.assign(static_cast<std::size_t>(cells) * degree, 0);
  q.assign(cells, 0.0f);
  for (int c = 0; c < cells; ++c) {
    Value density = 0.0f;
    in >> density;
    q[c] = density;
    for (int k = 0; k < degree; ++k) {
      int nb = -1;
      Value a = 0.0f, b = 0.0f, d = 0.0f;
      in >> nb >> a >> b >> d;
      nbr[c * degree + k] = (nb >= 0 && nb < cells) ? nb : c;
    }
    if (!in) return false;
  }
  return true;
}
void Run(std::vector<int> &nbr, std::vector<Value> &q, std::vector<Value> &next,
         int cells, int degree, int iterations) {
  int *pn = nbr.data();
  Value *pq = q.data(), *px = next.data();
#if DWARF_USE_TARGET
#pragma omp target data map(to : pn[0 : cells * degree]) map(tofrom : pq[0 : cells]) map(alloc : px[0 : cells])
  {
#endif
    for (int it = 0; it < iterations; ++it) {
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for schedule(static)
#endif
      for (int c = 0; c < cells; ++c) {
        Value flux = 0.0f;
        for (int k = 0; k < degree; ++k) flux += pq[pn[c * degree + k]] - pq[c];
        px[c] = pq[c] + 0.05f * flux;
      }
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for schedule(static)
#endif
      for (int c = 0; c < cells; ++c) pq[c] = px[c];
    }
#if DWARF_USE_TARGET
  }
#endif
}
}
int main(int argc, char **argv) {
  auto opts = dwarfs_omp::ParseOptions(argc, argv);
  int cells = static_cast<int>(opts.size), degree = static_cast<int>(opts.secondary_size);
  std::vector<int> nbr(cells * degree);
  std::vector<Value> q(cells), next(cells);
  if (!opts.input_file.empty()) {
    if (!LoadDomn(opts.input_file, nbr, q, cells, degree)) return 3;
    next.resize(cells);
  } else {
    for (int c = 0; c < cells; ++c) { q[c] = dwarfs_omp::InitValue(c); for (int k = 0; k < degree; ++k) nbr[c * degree + k] = (c + k * 17 + 1) % cells; }
  }
  auto times = dwarfs_omp::AdaptiveRun(opts,
      [](){},
      [&]() { Run(nbr, q, next, cells, degree, opts.iterations); });
  bool valid = std::isfinite(q[cells / 2]);
  dwarfs_omp::PrintResult("cfd", "unstructured-grids", q.size(),
                          static_cast<std::size_t>(degree), opts.iterations,
                          static_cast<int>(times.size()), times, valid);
  return valid ? 0 : 2;
}
