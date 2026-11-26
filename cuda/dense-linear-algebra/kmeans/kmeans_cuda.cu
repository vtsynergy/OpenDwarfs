#include "dwarf_launch.h"  /* [FIX-OCL-34] */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <assert.h>
#include "dwarf_cuda.h"
#include "kmeans_kernel.cuh"

#include <omp.h>

#define THREADS_PER_DIM 16
#define BLOCKS_PER_DIM 16
#define AOCL_ALIGNMENT 64

//#define BLOCK_DELTA_REDUCE
//#define BLOCK_CENTER_REDUCE

#define CPU_DELTA_REDUCE
#define CPU_CENTER_REDUCE

extern "C"
int setup(int argc, char** argv);									/* function prototype */
// GLOBAL!!!!!
unsigned int num_threads_perdim = THREADS_PER_DIM;					/* sqrt(256) -- see references for this choice */
unsigned int num_blocks_perdim = BLOCKS_PER_DIM;					/* temporary */
unsigned int num_threads = num_threads_perdim*num_threads_perdim;	/* number of threads */
unsigned int num_blocks = num_blocks_perdim*num_blocks_perdim;		/* number of blocks */

/* _d denotes it resides on the device */
int    *membership_new;												/* newly assignment membership */
float *feature_d;													/* inverted data array */
float *feature_flipped_d;											/* original (not inverted) data array */
int   *membership_d;												/* membership on the device */
float  *block_new_centers;											/* sum of points in a cluster (per block) */
float *clusters_d;													/* cluster centers on the device */

/* image memory */
/*cl_mem t_features;
  cl_mem t_features_flipped;
  cl_mem t_clusters;*/


/* -------------- initCL() -------------- */
/* kept under this name because kmeans.c calls it; the device was selected in main() and kernels are compiled in */
	extern "C"
void initCL()
{
	num_threads = num_threads_perdim*num_threads_perdim;
}
/* -------------- initCL() end -------------- */

/* -------------- allocateMemory() ------------------- */
/* allocate device memory, calculate number of blocks and threads, and invert the data array */
	extern "C"
void allocateMemory(int npoints, int nfeatures, int nclusters, float **features)
{	

	num_blocks = npoints / num_threads;
	if (npoints % num_threads > 0)		/* defeat truncation */
		num_blocks++;

	num_blocks_perdim = sqrt((double) num_blocks);
	while (num_blocks_perdim * num_blocks_perdim < num_blocks)	// defeat truncation (should run once)
		num_blocks_perdim++;

	num_blocks = num_blocks_perdim*num_blocks_perdim;

	/* allocate memory for memory_new[] and initialize to -1 (host) */
	membership_new = (int*) memalign ( AOCL_ALIGNMENT, npoints * sizeof(int));
	//membership_new = (int*) malloc(npoints * sizeof(int));
	for(int i=0;i<npoints;i++) {
		membership_new[i] = -1;
	}

	/* allocate memory for block_new_centers[] (host) */
	block_new_centers = (float *) memalign ( AOCL_ALIGNMENT, nclusters*nfeatures*sizeof(float));
        //block_new_centers = (float *) malloc(nclusters*nfeatures*sizeof(float));

	/* allocate memory for feature_flipped_d[][], feature_d[][] (device) */
	CUDA_CHECK(cudaMalloc((void**)&feature_flipped_d, npoints*nfeatures*sizeof(float)));
	CUDA_CHECK(cudaMemcpy(feature_flipped_d, features[0], npoints*nfeatures*sizeof(float), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaDeviceSynchronize());
	CUDA_CHECK(cudaMalloc((void**)&feature_d, npoints*nfeatures*sizeof(float)));

	/* invert the data array (kernel execution): one point per thread, DWARF_LS threads per block */
	invert_mapping<<<DWARF_NBLOCKS(npoints), DWARF_LS>>>(feature_flipped_d, feature_d, npoints, nfeatures);
	CUDA_SYNC_LAUNCH();

	/* allocate memory for membership_d[] and clusters_d[][] (device) */
	CUDA_CHECK(cudaMalloc((void**)&membership_d, npoints*sizeof(int)));
	CUDA_CHECK(cudaMalloc((void**)&clusters_d, nclusters*nfeatures*sizeof(float)));
}
/* -------------- allocateMemory() end ------------------- */

/* -------------- deallocateMemory() ------------------- */
/* free host and device memory */
	extern "C"
void deallocateMemory()
{
	free(membership_new);
	free(block_new_centers);
	CUDA_CHECK(cudaFree(feature_d));
	CUDA_CHECK(cudaFree(feature_flipped_d));
	CUDA_CHECK(cudaFree(membership_d));
	CUDA_CHECK(cudaFree(clusters_d));
}
/* -------------- deallocateMemory() end ------------------- */

/* counterpart of releaseCL() [FIX-OCL-33]: nothing to release in CUDA */
	extern "C"
void releaseCL()
{
}



////////////////////////////////////////////////////////////////////////////////
// Program main																  //

	int
main( int argc, char** argv) 
{
	// as done in the CUDA start/help document provided
	dwarf_cuda_init(&argc, &argv);
	setup(argc, argv);   
    printf("Setup done\n");
}

//																			  //
////////////////////////////////////////////////////////////////////////////////


/* ------------------- kmeansCuda() ------------------------ */    
extern "C"
	int	// delta -- had problems when return value was of float type
kmeansCuda(float  **feature,				/* in: [npoints][nfeatures] */
		int      nfeatures,				/* number of attributes for each point */
		int      npoints,				/* number of data points */
		int      nclusters,				/* number of clusters */
		int     *membership,				/* which cluster the point belongs to */
		float  **clusters,				/* coordinates of cluster centers */
		int     *new_centers_len,		/* number of elements in each cluster */
		float  **new_centers				/* sum of elements in each cluster */
	  )
{
	int delta = 0;			/* if point has moved */
	int i,j;				/* counters */


	/* copy membership (host to device) */
	CUDA_CHECK(cudaMemcpy(membership_d, membership_new, npoints*sizeof(int), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaDeviceSynchronize());

	/* copy clusters (host to device) */
	CUDA_CHECK(cudaMemcpy(clusters_d, clusters[0], nclusters*nfeatures*sizeof(float), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaDeviceSynchronize());

	/* one point per thread: DWARF_LS threads per block, ceil(npoints/DWARF_LS) blocks */
	kmeansPoint<<<DWARF_NBLOCKS(npoints), DWARF_LS>>>(feature_d, feature_flipped_d, nfeatures, npoints, nclusters, membership_d, clusters_d);
	CUDA_SYNC_LAUNCH();

	/* copy back membership (device to host) */
	CUDA_CHECK(cudaMemcpy(membership_new, membership_d, npoints*sizeof(int), cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaDeviceSynchronize());


	/* for each point, sum data points in each cluster
	   and see if membership has changed:
	   if so, increase delta and change old membership, and update new_centers;
	   otherwise, update new_centers */
	delta = 0;
	for (i = 0; i < npoints; i++)
	{		
		int cluster_id = membership_new[i];
		new_centers_len[cluster_id]++;
		if (membership_new[i] != membership[i])
		{
#ifdef CPU_DELTA_REDUCE
			delta++;
#endif
			membership[i] = membership_new[i];
		}
#ifdef CPU_CENTER_REDUCE
		for (j = 0; j < nfeatures; j++)
		{			
			new_centers[cluster_id][j] += feature[i][j];
		}
#endif
	}


#ifdef BLOCK_DELTA_REDUCE	
	/*** calculate global sums from per block sums for delta and the new centers ***/    

	//debug
	//printf("\t \t reducing %d block sums to global sum \n",num_blocks_perdim * num_blocks_perdim);
	for(i = 0; i < num_blocks_perdim * num_blocks_perdim; i++) {
		//printf("block %d delta is %d \n",i,block_deltas_h[i]);
		delta += block_deltas_h[i];
	}

#endif
#ifdef BLOCK_CENTER_REDUCE	

	for(int j = 0; j < nclusters;j++) {
		for(int k = 0; k < nfeatures;k++) {
			block_new_centers[j*nfeatures + k] = 0.f;
		}
	}

	for(i = 0; i < num_blocks_perdim * num_blocks_perdim; i++) {
		for(int j = 0; j < nclusters;j++) {
			for(int k = 0; k < nfeatures;k++) {
				block_new_centers[j*nfeatures + k] += block_clusters_h[i * nclusters*nfeatures + j * nfeatures + k];
			}
		}
	}


#ifdef CPU_CENTER_REDUCE
	//debug
	/*for(int j = 0; j < nclusters;j++) {
	  for(int k = 0; k < nfeatures;k++) {
	  if(new_centers[j][k] >	1.001 * block_new_centers[j*nfeatures + k] || new_centers[j][k] <	0.999 * block_new_centers[j*nfeatures + k]) {
	  printf("\t \t for %d:%d, normal value is %e and gpu reduced value id %e \n",j,k,new_centers[j][k],block_new_centers[j*nfeatures + k]);
	  }
	  }
	  }*/
#endif

#ifdef BLOCK_CENTER_REDUCE
	for(int j = 0; j < nclusters;j++) {
		for(int k = 0; k < nfeatures;k++)
			new_centers[j][k]= block_new_centers[j*nfeatures + k];		
	}
#endif

#endif

	return delta;

}
/* ------------------- kmeansCuda() end ------------------------ */    

