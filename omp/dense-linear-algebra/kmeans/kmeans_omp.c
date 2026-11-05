/* Device layer of K-means -- OpenMP target-offload replacement of kmeans_opencl.cpp.
 * The Rodinia host code (kmeans.c, cluster.c, kmeans_clustering.c, rmse.c) is unchanged apart from the
 * timing hooks; the entry points the host code calls are the same: initCL / allocateMemory /
 * kmeansCuda / deallocateMemory.
 *
 * Kernels (kmeans_opencl_kernel.cl): invert_mapping (row-major -> column-major copy of the features)
 * and kmeansPoint (nearest cluster centre per point; the membership of every point is written).
 * CPU_DELTA_REDUCE / CPU_CENTER_REDUCE are defined in the OpenCL build, so delta and the new centre
 * sums are accumulated on the host exactly as there. The GPU_* reduction variants of the kernel are
 * #ifdef'd out in the OpenCL kernel and are not part of this translation either.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <float.h>
#include <omp.h>
#include "kmeans.h"
#include "ocd_compat.h"

#define THREADS_PER_DIM 16
#define BLOCKS_PER_DIM 16
#define CPU_DELTA_REDUCE
#define CPU_CENTER_REDUCE

int setup(int argc, char **argv);

static unsigned int num_threads_perdim = THREADS_PER_DIM; /* sqrt(256) */
static unsigned int num_blocks_perdim = BLOCKS_PER_DIM;
static unsigned int num_threads = THREADS_PER_DIM * THREADS_PER_DIM;
static unsigned int num_blocks = BLOCKS_PER_DIM * BLOCKS_PER_DIM;

static int *membership_new;       /* newly assigned membership (host) */
static float *feature_d;          /* inverted data array (device) */
static float *feature_flipped_d;  /* original data array (device) */
static int *membership_d;         /* membership (device) */
static float *clusters_d;         /* cluster centres (device) */
static float *block_new_centers;  /* host scratch */
static int m_npoints, m_nfeatures, m_nclusters;

void initCL(void) {
  /* OpenCL: platform/device discovery and kernel build. Nothing to do for OpenMP offload. */
}

/* allocate device memory, calculate number of blocks and threads, and invert the data array */
void allocateMemory(int npoints, int nfeatures, int nclusters, float **features) {
  size_t NP = (size_t)npoints, NF = (size_t)nfeatures;
  m_npoints = npoints; m_nfeatures = nfeatures; m_nclusters = nclusters;
  num_blocks = npoints / num_threads;
  if (npoints % num_threads > 0) num_blocks++; /* defeat truncation */
  num_blocks_perdim = (unsigned int)sqrt((double)num_blocks);
  while (num_blocks_perdim * num_blocks_perdim < num_blocks) num_blocks_perdim++;
  num_blocks = num_blocks_perdim * num_blocks_perdim;

  /* membership_new[] initialised to -1 (host); block_new_centers[] (host) */
  membership_new = (int *)malloc(npoints * sizeof(int));
  for (int i = 0; i < npoints; i++) membership_new[i] = -1;
  block_new_centers = (float *)malloc(nclusters * nfeatures * sizeof(float));

  /* device buffers: "device" pointers are host allocations that carry the OpenMP mapping */
  feature_flipped_d = features[0];                       /* original layout, H2D below */
  feature_d = (float *)malloc(NP * NF * sizeof(float));  /* inverted layout, produced on the device */
  membership_d = membership_new;                         /* host mirror of the device membership */
  clusters_d = (float *)malloc((size_t)nclusters * NF * sizeof(float));

  float *ff = feature_flipped_d, *fd = feature_d;
  int *md = membership_d; float *cd = clusters_d;
#pragma omp target enter data map(to : ff[0 : NP * NF]) map(alloc : fd[0 : NP * NF], md[0 : NP], cd[0 : (size_t)nclusters * NF])

  /* invert_mapping kernel: output[point + npoints*i] = input[point*nfeatures + i]
   * global size = num_blocks*num_threads work-items, guarded by point_id < npoints */
  size_t global = (size_t)num_blocks * num_threads;
#pragma omp target teams distribute parallel for map(present, alloc : ff[0 : NP * NF], fd[0 : NP * NF])
  for (size_t point_id = 0; point_id < global; point_id++) {
    if (point_id < NP)
      for (size_t i = 0; i < NF; i++) fd[point_id + NP * i] = ff[point_id * NF + i];
  }
}

/* free host and device memory */
void deallocateMemory(void) {
  size_t NP = (size_t)m_npoints, NF = (size_t)m_nfeatures;
  float *ff = feature_flipped_d, *fd = feature_d, *cd = clusters_d;
  int *md = membership_d;
#pragma omp target exit data map(delete : ff[0 : NP * NF], fd[0 : NP * NF], md[0 : NP], cd[0 : (size_t)m_nclusters * NF])
  free(membership_new);
  free(block_new_centers);
  free(feature_d);
  free(clusters_d);
}

int kmeansCuda(float **feature,       /* in: [npoints][nfeatures] */
               int nfeatures,         /* number of attributes for each point */
               int npoints,           /* number of data points */
               int nclusters,         /* number of clusters */
               int *membership,       /* which cluster the point belongs to */
               float **clusters,      /* coordinates of cluster centers */
               int *new_centers_len,  /* number of elements in each cluster */
               float **new_centers    /* sum of elements in each cluster */
) {
  int delta = 0; /* if point has moved */
  int i, j;
  size_t NP = (size_t)npoints, NF = (size_t)nfeatures, NC = (size_t)nclusters;
  float *fd = feature_d, *cd = clusters_d;
  int *md = membership_d;

  /* copy membership (host to device) and clusters (host to device) */
  memcpy(cd, clusters[0], NC * NF * sizeof(float));
#pragma omp target update to(md[0 : NP], cd[0 : NC * NF])

  /* kmeansPoint kernel: one work-item per point. The OpenCL global range is
   * (num_blocks_perdim*256) x num_blocks_perdim with local size (256,1); point_id =
   * block_id*256 + local_id covers [0, num_blocks*256) and is guarded by point_id < npoints. */
  size_t total = (size_t)num_blocks * num_threads;
#pragma omp target teams distribute parallel for map(present, alloc : fd[0 : NP * NF], cd[0 : NC * NF], md[0 : NP])
  for (size_t point_id = 0; point_id < total; point_id++) {
    int index = -1;
    if (point_id < NP) {
      float min_dist = FLT_MAX;
      for (int c = 0; c < (int)NC; c++) {
        int cluster_base_index = c * (int)NF; /* base index of cluster centers */
        float ans = 0.0;                      /* Euclidean distance square */
        for (int f = 0; f < (int)NF; f++) {
          int addr = (int)point_id + f * (int)NP; /* appropriate index of data point */
          float diff = fd[addr] - cd[cluster_base_index + f];
          ans += diff * diff;
        }
        float dist = ans;
        if (dist < min_dist) { min_dist = dist; index = c; }
      }
      md[point_id] = index; /* assign the membership to object point_id */
    }
  }

  /* copy back membership (device to host) */
#pragma omp target update from(md[0 : NP])

  /* for each point, sum data points in each cluster and see if membership has changed:
     if so, increase delta and change old membership, and update new_centers;
     otherwise, update new_centers */
  delta = 0;
  for (i = 0; i < npoints; i++) {
    int cluster_id = membership_new[i];
    new_centers_len[cluster_id]++;
    if (membership_new[i] != membership[i]) {
#ifdef CPU_DELTA_REDUCE
      delta++;
#endif
      membership[i] = membership_new[i];
    }
#ifdef CPU_CENTER_REDUCE
    for (j = 0; j < nfeatures; j++) new_centers[cluster_id][j] += feature[i][j];
#endif
  }
  return delta;
}

int main(int argc, char **argv) {
  int rc = setup(argc, argv);
  printf("Setup done\n");
  return rc;
}
