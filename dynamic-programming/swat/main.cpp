#include "simple_motifs.hpp"
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

namespace {
bool LoadSequence(const std::string &path, std::vector<int> &seq) {
  std::ifstream in(path);
  if (!in) return false;
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line[0] == '>') continue;
    for (char ch : line) {
      switch (ch) {
        case 'a': case 'A': seq.push_back(0); break;
        case 'c': case 'C': seq.push_back(1); break;
        case 'g': case 'G': seq.push_back(2); break;
        case 't': case 'T': seq.push_back(3); break;
        default: break;
      }
    }
  }
  return !seq.empty();
}
void Run(std::vector<int> &query, std::vector<int> &db, std::vector<int> &score,
         int rows, int cols) {
  int *pq = query.data();
  int *pd = db.data();
  int *ps = score.data();
  for (int wave = 2; wave <= rows + cols - 2; ++wave) {
    int lo = std::max(1, wave - (cols - 1));
    int hi = std::min(rows - 1, wave - 1);
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for map(to : pq[0 : rows - 1], pd[0 : cols - 1]) map(tofrom : ps[0 : rows * cols]) firstprivate(rows, cols, wave, lo, hi) thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for schedule(static)
#endif
    for (int i = lo; i <= hi; ++i) {
      int j = wave - i;
      int match = ps[(i - 1) * cols + j - 1] + (pq[i - 1] == pd[j - 1] ? 3 : -2);
      int del = ps[(i - 1) * cols + j] - 1;
      int ins = ps[i * cols + j - 1] - 1;
      ps[i * cols + j] = std::max(0, std::max(match, std::max(del, ins)));
    }
  }
}

void RunSerial(const std::vector<int> &query, const std::vector<int> &db, std::vector<int> &score, int rows, int cols) {
  const int *pq = query.data();
  const int *pd = db.data();
  int *ps = score.data();
  for (int i = 1; i < rows; ++i) {
    for (int j = 1; j < cols; ++j) {
      int match = ps[(i - 1) * cols + j - 1] + (pq[i - 1] == pd[j - 1] ? 3 : -2);
      int del = ps[(i - 1) * cols + j] - 1;
      int ins = ps[i * cols + j - 1] - 1;
      ps[i * cols + j] = std::max(0, std::max(match, std::max(del, ins)));
    }
  }
}

bool Verify(const std::vector<int> &score_gpu, const std::vector<int> &query, const std::vector<int> &db, int rows, int cols) {
  std::vector<int> score_ref(static_cast<std::size_t>(rows) * cols, 0);
  RunSerial(query, db, score_ref, rows, cols);
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
  std::vector<int> query, db;
  if (!opts.input_file.empty()) {
    if (!LoadSequence(opts.input_file, query)) return 3;
    if (!opts.input_file2.empty()) {
      if (!LoadSequence(opts.input_file2, db)) return 3;
    } else {
      db = query;
    }
  } else {
    int n = static_cast<int>(opts.size);
    query.resize(n);
    db.resize(n);
    for (int i = 0; i < n; ++i) {
      query[i] = (i * 3 + 1) & 3;
      db[i] = (i * 5 + 2) & 3;
    }
  }
  int rows = static_cast<int>(query.size()) + 1;
  int cols = static_cast<int>(db.size()) + 1;
  std::vector<int> score(static_cast<std::size_t>(rows) * cols, 0);
  auto times = dwarfs_omp::AdaptiveRun(opts,
      [](){},
      [&]() { for (int it = 0; it < opts.iterations; ++it) { std::fill(score.begin(), score.end(), 0); Run(query, db, score, rows, cols); } });

  bool valid = true;
  if (opts.verify) {
    valid = Verify(score, query, db, rows, cols);
  }
  dwarfs_omp::PrintResult("swat", "dynamic-programming", query.size(),
                          db.size(), opts.iterations,
                          static_cast<int>(times.size()), times, valid);
  return valid ? 0 : 2;
}
