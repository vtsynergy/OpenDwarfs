#include "dwarf_common.hpp"

#include <cmath>
#include <vector>

namespace {
using Value = float;

// Symmetric positive-definite system A x = b. A has a dominant diagonal (= n) with
// small symmetric off-diagonals, so the system is very well conditioned and CG drives
// the residual down in a handful of iterations regardless of size.
void Init(std::vector<Value> &a, std::vector<Value> &b, int n) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      a[i * n + j] = (i == j) ? static_cast<Value>(n)
                              : 1.0f / (1.0f + static_cast<Value>(std::abs(i - j)));
    }
    b[i] = 1.0f;
  }
}

// Conjugate gradient. One iteration is six canonical loops over n:
//   matvec  q = A p        (S-per-row inner sum)
//   dot     pq = p . q     (reduction)
//   axpy    x  = x + a p
//   axpy    r  = r - a q
//   dot     rs = r . r     (reduction)
//   xpby    p  = r + b p
// The vectors stay resident on the device across all iterations; only the two scalar
// reductions cross back to the host, where alpha and beta are formed.
void Run(std::vector<Value> &a, std::vector<Value> &b, std::vector<Value> &x,
         std::vector<Value> &r, std::vector<Value> &p, std::vector<Value> &q,
         int n, int iterations) {
  Value *pa = a.data();
  Value *px = x.data();
  Value *pr = r.data();
  Value *pp = p.data();
  Value *pq = q.data();

  // x0 = 0, r0 = b - A x0 = b, p0 = r0  (host init).
  for (int i = 0; i < n; ++i) {
    px[i] = 0.0f;
    pr[i] = b[i];
    pp[i] = b[i];
  }
  Value rsold = 0.0f;
  for (int i = 0; i < n; ++i) rsold += pr[i] * pr[i];

#if DWARF_USE_TARGET
#pragma omp target data map(to : pa[0 : n * n]) map(tofrom : px[0 : n]) \
                        map(to : pr[0 : n], pp[0 : n]) map(alloc : pq[0 : n])
  {
#endif
    for (int it = 0; it < iterations; ++it) {
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for schedule(static)
#endif
      for (int i = 0; i < n; ++i) {
        Value s = 0.0f;
        for (int j = 0; j < n; ++j) s += pa[i * n + j] * pp[j];
        pq[i] = s;
      }

      Value pqd = 0.0f;
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for reduction(+ : pqd) thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for reduction(+ : pqd) schedule(static)
#endif
      for (int i = 0; i < n; ++i) pqd += pp[i] * pq[i];

      Value alpha = rsold / (pqd + 1.0e-30f);

#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for schedule(static)
#endif
      for (int i = 0; i < n; ++i) px[i] += alpha * pp[i];

#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for schedule(static)
#endif
      for (int i = 0; i < n; ++i) pr[i] -= alpha * pq[i];

      Value rsnew = 0.0f;
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for reduction(+ : rsnew) thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for reduction(+ : rsnew) schedule(static)
#endif
      for (int i = 0; i < n; ++i) rsnew += pr[i] * pr[i];

      Value beta = rsnew / (rsold + 1.0e-30f);

#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for schedule(static)
#endif
      for (int i = 0; i < n; ++i) pp[i] = pr[i] + beta * pp[i];

      rsold = rsnew;
    }
#if DWARF_USE_TARGET
  }
#endif
}

// Verify that x approximately solves A x = b: the relative infinity-norm residual
// should be small once CG has run. The tolerance is loose because the number of
// iterations is a run parameter, not tuned for exact convergence.
bool Verify(const std::vector<Value> &a, const std::vector<Value> &b,
            const std::vector<Value> &x, int n) {
  Value rmax = 0.0f, bmax = 0.0f;
  for (int i = 0; i < n; ++i) {
    Value ax = 0.0f;
    for (int j = 0; j < n; ++j) ax += a[i * n + j] * x[j];
    rmax = std::max(rmax, std::abs(ax - b[i]));
    bmax = std::max(bmax, std::abs(b[i]));
  }
  return rmax <= 1.0e-2f * (bmax + 1.0e-30f);
}
}  // namespace

int main(int argc, char **argv) {
  const auto opts = dwarfs_omp::ParseOptions(argc, argv);
  const int n = static_cast<int>(opts.size);
  const int iters = std::max(1, opts.iterations);
  std::vector<Value> a(static_cast<size_t>(n) * n), b(n), x(n), r(n), p(n), q(n);
  bool valid = true;

  Init(a, b, n);
  auto times = dwarfs_omp::AdaptiveRun(opts,
      [&]() { std::fill(x.begin(), x.end(), 0.0f); },
      [&]() { Run(a, b, x, r, p, q, n, iters); });
  if (opts.verify) {
    valid = Verify(a, b, x, n);
  }
  dwarfs_omp::PrintResult("cg", "dense-linear-algebra", opts.size,
                          opts.secondary_size, opts.iterations,
                          static_cast<int>(times.size()), times, valid);
  return valid ? 0 : 2;
}
