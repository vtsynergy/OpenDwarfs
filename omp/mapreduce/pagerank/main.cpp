/*
 ** The MIT License (MIT)
 **
 ** Copyright (c) 2014, Erick Lavoie, Faiz Khan, Sujay Kathrotia, Vincent
 ** Foley-Bourgon, Laurie Hendren
 **
 ** Permission is hereby granted, free of charge, to any person obtaining a copy
 **of this software and associated documentation files (the "Software"), to deal
 ** in the Software without restriction, including without limitation the rights
 ** to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 ** copies of the Software, and to permit persons to whom the Software is
 ** furnished to do so, subject to the following conditions:
 **
 ** The above copyright notice and this permission notice shall be included in all
 ** copies or substantial portions of the Software.
 **
 ** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 ** IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 ** FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 ** AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 ** LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 ** OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 ** SOFTWARE.
 **
 **/
#include "dwarf_common.hpp"
#include <omp.h>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iostream>

#define D_FACTOR (0.85f)
#ifndef BLOCK_SIZE
#define BLOCK_SIZE DWARF_THREAD_LIMIT
#endif

// generates an array of random pages and their links
int *random_pages(int n, unsigned int *noutlinks, int divisor){
  int i, j, k;
  int *pages = (int*) calloc((size_t)n * n, sizeof(int)); // matrix 1 means link from j->i, 0 elsewhere

  if (divisor <= 0) {
    fprintf(stderr, "ERROR: Invalid divisor '%d' for random initialization, divisor should be greater or equal to 1\n", divisor);
    exit(1);
  }

  for(i=0; i<n; ++i){
    noutlinks[i] = 0;
    for(j=0; j<n; ++j){
      if(i!=j && (abs(rand())%divisor == 0)){
        pages[(size_t)i*n+j] = 1;
        noutlinks[i] += 1;
      }
    }

    // the case with no outlinks is avoided
    if(noutlinks[i] == 0){
      if (n == 1) k = 0;
      else do { k = abs(rand()) % n; } while ( k == i);
      pages[(size_t)i*n + k] = 1;
      noutlinks[i] = 1;
    }
  }
  return pages;
}

void map_kernel(const int* pages, const float* page_ranks, float* maps, const unsigned int* noutlinks, int n) {
  size_t block_size  = n < BLOCK_SIZE ? n : BLOCK_SIZE;
#if DWARF_USE_TARGET
  #pragma omp target teams distribute parallel for thread_limit(block_size)
#else
  #pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < n; i++) {
    float outbound_rank = page_ranks[i]/(float)noutlinks[i];
    for(int j=0; j<n; ++j) {
      maps[(size_t)i*n+j] = pages[(size_t)i*n+j]*outbound_rank;
    }
  }
}

void reduce_kernel(float* page_ranks, const float* maps, int n) {
  size_t block_size  = n < BLOCK_SIZE ? n : BLOCK_SIZE;
#if DWARF_USE_TARGET
  #pragma omp target teams distribute parallel for thread_limit(block_size)
#else
  #pragma omp parallel for schedule(static)
#endif
  for (int j = 0; j < n; j++) {
    float new_rank = 0.0f;
    for(int i=0; i< n; ++i) {
      new_rank += maps[(size_t)i*n + j];
    }
    page_ranks[j] = ((1.f-D_FACTOR)/n)+(D_FACTOR*new_rank);
  }
}

#if DWARF_USE_TARGET
// NOTE: this loop must live in a plain function, not a lambda. Clang 21.1's
// nvptx64 offload codegen silently fails to copy `target data`-mapped
// results back to the host when the `target data` construct (and the
// target regions nested inside it via map_kernel/reduce_kernel) is lexically
// inside a C++ lambda with reference captures -- e.g. the `[&](){...}`
// RunFn passed to dwarfs_omp::AdaptiveRun below. No diagnostic, no crash:
// the host buffer just silently keeps its pre-kernel value. Verified with a
// minimal ~70-line reproducer outside this codebase: identical map/reduce
// code inside a plain function updates page_ranks correctly; the same code
// inside a lambda called the same way does not. Extracting the device loop
// into this free function (called *by* the lambda, rather than living
// *inside* it) works around it.
void run_device_iterations(const int* pages, float* page_ranks, float* maps,
                            const unsigned int* noutlinks, int n, int pagerank_iters) {
  #pragma omp target data map(to: pages[0:(size_t)n*n], noutlinks[0:n]) \
                          map(tofrom: page_ranks[0:n]) \
                          map(alloc: maps[0:(size_t)n*n])
  {
    for (int it = 0; it < pagerank_iters; ++it) {
      map_kernel(pages, page_ranks, maps, noutlinks, n);
      reduce_kernel(page_ranks, maps, n);
    }
  }
}
#endif

bool verify(const int* pages, const unsigned int* noutlinks, const float* page_ranks_gpu, int n, int iterations) {
  std::vector<float> page_ranks_ref(n, 1.0f / (float) n);
  std::vector<float> maps_ref((size_t)n * n);

  for (int t = 0; t < iterations; ++t) {
    for (int i = 0; i < n; ++i) {
      float outbound_rank = page_ranks_ref[i] / (float)noutlinks[i];
      for (int j = 0; j < n; ++j) {
        maps_ref[(size_t)i*n+j] = pages[(size_t)i*n+j] * outbound_rank;
      }
    }
    for (int j = 0; j < n; ++j) {
      float new_rank = 0.0f;
      for (int i = 0; i < n; ++i) {
        new_rank += maps_ref[(size_t)i*n+j];
      }
      page_ranks_ref[j] = ((1.f-D_FACTOR)/n)+(D_FACTOR*new_rank);
    }
  }

  float max_diff = 0.0f;
  for (int i = 0; i < n; ++i) {
    float diff = std::abs(page_ranks_gpu[i] - page_ranks_ref[i]);
    if (diff > max_diff) {
      max_diff = diff;
    }
  }
  return max_diff < 1e-3f;
}

int main(int argc, char **argv) {
  const auto opts = dwarfs_omp::ParseOptions(argc, argv);
  const int n = static_cast<int>(opts.size);
  const int divisor = static_cast<int>(opts.secondary_size);  // link density
  const int pagerank_iters = opts.iterations;   // inner PR convergence iters

  std::vector<unsigned int> noutlinks(n, 0);
  int* pages = random_pages(n, noutlinks.data(), divisor);
  std::vector<float> page_ranks(n);
  std::vector<float> maps((size_t)n * n);

  bool valid = true;

  auto times = dwarfs_omp::AdaptiveRun(opts,
      [&]() { std::fill(page_ranks.begin(), page_ranks.end(), 1.0f / (float)n); },
      [&]() {
#if DWARF_USE_TARGET
        run_device_iterations(pages, page_ranks.data(), maps.data(), noutlinks.data(), n, pagerank_iters);
#else
        for (int it = 0; it < pagerank_iters; ++it) {
          map_kernel(pages, page_ranks.data(), maps.data(), noutlinks.data(), n);
          reduce_kernel(page_ranks.data(), maps.data(), n);
        }
#endif
      });

  if (opts.verify) {
    valid = verify(pages, noutlinks.data(), page_ranks.data(), n, pagerank_iters);
  }

  dwarfs_omp::PrintResult("pagerank", "map-reduce",
                          opts.size, opts.secondary_size,
                          opts.iterations, static_cast<int>(times.size()), times, valid);

  free(pages);
  return valid ? 0 : 2;
}
