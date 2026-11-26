/* Device layer of K-means -- OpenACC translation of opencl/dense-linear-algebra/kmeans/
 * {kmeans_opencl.cpp, kmeans_opencl_kernel.cl}. The Rodinia host code (kmeans.c, cluster.c,
 * kmeans_clustering.c, rmse.c) calls the same entry points: initCL / allocateMemory / kmeansCuda /
 * deallocateMemory.
 *
 * Kernels: invert_mapping (row-major -> column-major copy of the features) and kmeansPoint (nearest
 * cluster centre per point). As in the OpenCL build (CPU_DELTA_REDUCE / CPU_CENTER_REDUCE), delta
 * and the new centre sums are accumulated on the host; the GPU_* reduction variants of the kernel
 * are #ifdef'd out there and are not translated.
 * Device buffers are created and released inside cluster(), i.e. inside the timed region, as the
 * OpenCL host does.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <float.h>
#include <openacc.h>
#include "kmeans.h"
#include "ocd_compat.h"
#include "dwarf_acc.h" /* [FIX-OCL-34] */

#define CPU_DELTA_REDUCE
#define CPU_CENTER_REDUCE

int setup(int argc, char **argv);

static int *membership_new;       /* newly assigned membership (host) */
static float *feature_d;          /* inverted data array (device) */
static float *feature_flipped_d;  /* original data array (device) */
static int *membership_d;         /* membership (device) */
static float *clusters_d;         /* cluster centres (device) */
static float *block_new_centers;  /* host scratch */
static int m_npoints, m_nfeatures, m_nclusters;

void initCL(void) {
  acc_init(acc_get_device_type()); /* device discovery/initialisation, outside the timed region */
}

/* allocate device memory and invert the data array */
void allocateMemory(int npoints, int nfeatures, int nclusters, float **features) {
  size_t NP = (size_t)npoints, NF = (size_t)nfeatures;
  m_npoints = npoints; m_nfeatures = nfeatures; m_nclusters = nclusters;

  /* membership_new[] initialised to -1 (host); block_new_centers[] (host) */
  membership_new = (int *)malloc(npoints * sizeof(int));
  for (int i = 0; i < npoints; i++) membership_new[i] = -1;
  block_new_centers = (float *)malloc(nclusters * nfeatures * sizeof(float));

  /* device buffers: each "device" pointer is a host allocation that names an OpenACC present-table
   * entry (on multicore it is the storage the kernels use) */
  feature_flipped_d = features[0];                      /* original layout */
  feature_d = (float *)malloc(NP * NF * sizeof(float)); /* inverted layout, produced by the kernel */
  membership_d = membership_new;
  clusters_d = (float *)malloc((size_t)nclusters * NF * sizeof(float));

  float *ff = feature_flipped_d, *fd = feature_d;
  int *md = membership_d;
  float *cd = clusters_d;
#pragma acc enter data create(ff[0 : NP * NF])
#pragma acc update device(ff[0 : NP * NF])
#pragma acc enter data create(fd[0 : NP * NF])

  /* invert_mapping kernel: output[point + npoints*i] = input[point*nfeatures + i];
   * OpenCL global = roundup(npoints, 128), local 128 (the loop bound replaces the point_id < npoints guard) */
#pragma acc parallel loop gang vector ACC_LAUNCH(NP) present(ff[0 : NP * NF], fd[0 : NP * NF])
  for (int point_id = 0; point_id < npoints; point_id++) {
#pragma acc loop seq
    for (int i = 0; i < nfeatures; i++) fd[point_id + NP * i] = ff[point_id * NF + i];
  }

#pragma acc enter data create(md[0 : NP], cd[0 : (size_t)nclusters * NF])
}

/* free host and device memory */
void deallocateMemory(void) {
  size_t NP = (size_t)m_npoints, NF = (size_t)m_nfeatures;
  float *ff = feature_flipped_d, *fd = feature_d, *cd = clusters_d;
  int *md = membership_d;
#pragma acc exit data delete(ff[0 : NP * NF], fd[0 : NP * NF], md[0 : NP], cd[0 : (size_t)m_nclusters * NF])
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
#pragma acc update device(md[0 : NP], cd[0 : NC * NF])

  /* kmeansPoint kernel: one work-item per point, local 128, global roundup(npoints, 128) */
#pragma acc parallel loop gang vector ACC_LAUNCH(NP) present(fd[0 : NP * NF], cd[0 : NC * NF], md[0 : NP])
  for (int point_id = 0; point_id < npoints; point_id++) {
    int index = -1;
    float min_dist = FLT_MAX;
#pragma acc loop seq
    for (int c = 0; c < nclusters; c++) {
      int cluster_base_index = c * nfeatures; /* base index of cluster centers */
      float ans = 0.0f;                       /* Euclidean distance square */
#pragma acc loop seq
      for (int f = 0; f < nfeatures; f++) {
        int addr = point_id + f * npoints; /* appropriate index of data point */
        float diff = fd[addr] - cd[cluster_base_index + f];
        ans += diff * diff;
      }
      float dist = ans;
      if (dist < min_dist) { min_dist = dist; index = c; }
    }
    md[point_id] = index; /* assign the membership to object point_id */
  }

  /* copy back membership (device to host) */
#pragma acc update self(md[0 : NP])

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
