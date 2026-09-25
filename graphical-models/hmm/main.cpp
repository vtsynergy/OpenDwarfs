#include "simple_motifs.hpp"
#include <algorithm>
#include <vector>

namespace {
using Value = float;
void Run(std::vector<Value> &trans, std::vector<Value> &emit,
         std::vector<Value> &prev, std::vector<Value> &next, int steps,
         int states) {
  Value *pt = trans.data(), *pe = emit.data(), *pp = prev.data(), *pn = next.data();
  for (int step = 0; step < steps; ++step) {
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for map(to : pt[0 : states * states], pe[0 : steps * states], pp[0 : states]) map(from : pn[0 : states]) firstprivate(step, states) thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for schedule(static)
#endif
    for (int s = 0; s < states; ++s) {
      Value best = -1.0e30f;
      for (int p = 0; p < states; ++p) best = std::max(best, pp[p] + pt[p * states + s]);
      pn[s] = best + pe[step * states + s];
    }
    std::swap(prev, next);
    pp = prev.data(); pn = next.data();
  }
}

void RunSerial(const std::vector<Value> &trans, const std::vector<Value> &emit,
               std::vector<Value> &prev, std::vector<Value> &next, int steps,
               int states) {
  const Value *pt = trans.data(), *pe = emit.data();
  Value *pp = prev.data(), *pn = next.data();
  for (int step = 0; step < steps; ++step) {
    for (int s = 0; s < states; ++s) {
      Value best = -1.0e30f;
      for (int p = 0; p < states; ++p) best = std::max(best, pp[p] + pt[p * states + s]);
      pn[s] = best + pe[step * states + s];
    }
    std::swap(prev, next);
    pp = prev.data(); pn = next.data();
  }
}

bool Verify(const std::vector<Value> &prev_gpu, const std::vector<Value> &trans,
            const std::vector<Value> &emit, int steps, int states) {
  std::vector<Value> prev_ref(states, 0.0f), next_ref(states);
  RunSerial(trans, emit, prev_ref, next_ref, steps, states);
  for (std::size_t i = 0; i < prev_gpu.size(); ++i) {
    if (!dwarfs_omp::NearlyEqual(prev_gpu[i], prev_ref[i], 1e-2)) {
      return false;
    }
  }
  return true;
}

}
int main(int argc, char **argv) {
  auto opts = dwarfs_omp::ParseOptions(argc, argv);
  int steps = static_cast<int>(opts.size), states = static_cast<int>(opts.secondary_size);
  std::vector<Value> trans(states * states), emit(steps * states), prev(states), next(states);
  for (int i = 0; i < states * states; ++i) trans[i] = -dwarfs_omp::InitValue(i);
  for (int i = 0; i < steps * states; ++i) emit[i] = -dwarfs_omp::InitValue(i, 113);
  std::fill(prev.begin(), prev.end(), 0.0f);
  auto times = dwarfs_omp::AdaptiveRun(opts,
      [](){},
      [&]() { for (int it = 0; it < opts.iterations; ++it) { std::fill(prev.begin(), prev.end(), 0.0f); Run(trans, emit, prev, next, steps, states); } });

  bool valid = true;
  if (opts.verify) {
    valid = Verify(prev, trans, emit, steps, states);
  }
  dwarfs_omp::PrintSimpleResult("hmm", "graphical-models", opts, times, valid);
  return valid ? 0 : 2;
}
