#include "simple_motifs.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {
using Value = float;
void Run(std::vector<Value> &x, std::vector<Value> &y, std::vector<Value> &z,
         std::vector<Value> &q, std::vector<Value> &sx, std::vector<Value> &sy,
         std::vector<Value> &sz, std::vector<Value> &pot, int atoms,
         int samples) {
  Value *px = x.data(), *py = y.data(), *pz = z.data(), *pq = q.data();
  Value *psx = sx.data(), *psy = sy.data(), *psz = sz.data(), *pp = pot.data();
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for map(to : px[0 : atoms], py[0 : atoms], pz[0 : atoms], pq[0 : atoms], psx[0 : samples], psy[0 : samples], psz[0 : samples]) map(from : pp[0 : samples]) thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for schedule(static)
#endif
  for (int s = 0; s < samples; ++s) {
    Value acc = 0.0f;
    for (int a = 0; a < atoms; ++a) {
      Value dx = psx[s] - px[a], dy = psy[s] - py[a], dz = psz[s] - pz[a];
      acc += pq[a] / (0.001f + dx * dx + dy * dy + dz * dz);
    }
    pp[s] = acc;
  }
}
bool LoadPqr(const std::string &path, std::vector<Value> &x,
             std::vector<Value> &y, std::vector<Value> &z,
             std::vector<Value> &q) {
  std::ifstream in(path);
  if (!in) return false;
  std::string line;
  while (std::getline(in, line)) {
    if (line.rfind("ATOM", 0) != 0 && line.rfind("HETATM", 0) != 0) continue;
    std::istringstream ls(line);
    std::vector<std::string> tok;
    std::string t;
    while (ls >> t) tok.push_back(t);
    if (tok.size() < 10) continue;
    x.push_back(std::stof(tok[tok.size() - 5]));
    y.push_back(std::stof(tok[tok.size() - 4]));
    z.push_back(std::stof(tok[tok.size() - 3]));
    q.push_back(std::stof(tok[tok.size() - 2]));
  }
  return !x.empty();
}
bool LoadVertices(const std::string &path, std::vector<Value> &x,
                  std::vector<Value> &y, std::vector<Value> &z) {
  std::ifstream in(path);
  if (!in) return false;
  std::string line;
  bool skipped_header = false;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ls(line);
    Value vx, vy, vz;
    if (!(ls >> vx)) continue;
    if (!skipped_header) {
      skipped_header = true;
      continue;
    }
    if (ls >> vy >> vz) {
      x.push_back(vx);
      y.push_back(vy);
      z.push_back(vz);
    }
  }
  return !x.empty();
}
std::string DefaultVertPath(std::string pqr) {
  const std::string suffix = ".pqr";
  if (pqr.size() >= suffix.size() &&
      pqr.compare(pqr.size() - suffix.size(), suffix.size(), suffix) == 0) {
    pqr.resize(pqr.size() - suffix.size());
  }
  return pqr + ".vert";
}
}
int main(int argc, char **argv) {
  auto opts = dwarfs_omp::ParseOptions(argc, argv);
  int atoms = static_cast<int>(opts.size), samples = static_cast<int>(opts.secondary_size);
  std::vector<Value> x(atoms), y(atoms), z(atoms), q(atoms), sx(samples), sy(samples), sz(samples), pot(samples);
  if (!opts.input_file.empty()) {
    x.clear(); y.clear(); z.clear(); q.clear();
    sx.clear(); sy.clear(); sz.clear();
    if (!LoadPqr(opts.input_file, x, y, z, q)) return 3;
    std::string vert = opts.input_file2.empty() ? DefaultVertPath(opts.input_file) : opts.input_file2;
    if (!LoadVertices(vert, sx, sy, sz)) return 3;
    atoms = static_cast<int>(x.size());
    samples = static_cast<int>(sx.size());
    pot.resize(samples);
  } else {
    for (int i = 0; i < atoms; ++i) { x[i] = dwarfs_omp::InitValue(i); y[i] = dwarfs_omp::InitValue(i, 113); z[i] = dwarfs_omp::InitValue(i, 109); q[i] = 0.5f + dwarfs_omp::InitValue(i, 101); }
    for (int s = 0; s < samples; ++s) { sx[s] = (s % 97) * 0.01f; sy[s] = (s % 89) * 0.01f; sz[s] = (s % 83) * 0.01f; }
  }
  auto times = dwarfs_omp::AdaptiveRun(opts,
      [](){},
      [&]() { for (int it = 0; it < opts.iterations; ++it) Run(x, y, z, q, sx, sy, sz, pot, atoms, samples); });
  bool valid = opts.verify ? std::isfinite(pot[samples / 2]) : true;
  dwarfs_omp::PrintResult("gem", "n-body-methods", x.size(), sx.size(),
                          opts.iterations, static_cast<int>(times.size()), times, valid);
  return valid ? 0 : 2;
}
