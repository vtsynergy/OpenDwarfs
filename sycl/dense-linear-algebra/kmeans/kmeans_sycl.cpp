#include "dwarf_launch.h"  /* [FIX-OCL-34] */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <assert.h>
#include <cfloat>
#include "dwarf_sycl.hpp"

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

static sycl::queue *Q;

/* ---- kernels (translated from kmeans_opencl_kernel.cl) ---- */
static void invert_mapping_kernel(sycl::queue &q, size_t global, float *input, float *output, int npoints, int nfeatures)
{
	q.parallel_for<class km_invert>(sycl::nd_range<1>(global, DWARF_LS), [=](sycl::nd_item<1> it) {
		int point_id = it.get_local_id(0) + it.get_local_range(0)*it.get_group(0);	/* id of thread */
		int i;
		if(point_id < npoints){
			for(i=0;i<nfeatures;i++)
				output[point_id + npoints*i] = input[point_id*nfeatures + i];
		}
	}).wait_and_throw();
}

static void kmeansPoint_kernel(sycl::queue &q, size_t global, float *features, int nfeatures, int npoints, int nclusters,
		int *membership, const float *clusters)
{
	q.parallel_for<class km_point>(sycl::nd_range<1>(global, DWARF_LS), [=](sycl::nd_item<1> it) {
		const unsigned int block_id = it.get_group(0);
		const unsigned int point_id = block_id*it.get_local_range(0) + it.get_local_id(0);
		int  index = -1;
		if (point_id < npoints)
		{
			int i, j;
			float min_dist = FLT_MAX;
			float dist;
			for (i=0; i<nclusters; i++) {
				int cluster_base_index = i*nfeatures;
				float ans=0.0;
				for (j=0; j < nfeatures; j++)
				{
					int addr = point_id + j*npoints;
					float diff = (features[addr] -
							clusters[cluster_base_index + j]);
					ans += diff*diff;
				}
				dist = ans;
				if (dist < min_dist) {
					min_dist = dist;
					index    = i;
				}
			}
		}
		if (point_id < npoints)
			membership[point_id] = index;
	}).wait_and_throw();
}

/* _d denotes it resides on the device */
int    *membership_new;												/* newly assignment membership */
float  *feature_d;													/* inverted data array */
float  *feature_flipped_d;											/* original (not inverted) data array */
int    *membership_d;												/* membership on the device */
float  *block_new_centers;											/* sum of points in a cluster (per block) */
float  *clusters_d;													/* cluster centers on the device */

/* -------------- initCL() -------------- */
	extern "C"
void initCL()
{
	static sycl::queue q = dsycl::make_queue();
	Q = &q;
	num_threads = num_threads_perdim*num_threads_perdim;
	dsycl::prebuild(q);   /* kernel JIT, as clBuildProgram: outside the timed region */
}
/* -------------- initCL() end -------------- */

/* -------------- allocateMemory() ------------------- */
/* allocate device memory, calculate number of blocks and threads, and invert the data array */
	extern "C"
void allocateMemory(int npoints, int nfeatures, int nclusters, float **features)
{
	sycl::queue &q = *Q;

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
	feature_flipped_d = dsycl::alloc<float>(q, (size_t)npoints*nfeatures);
	dsycl::write(q, feature_flipped_d, features[0], (size_t)npoints*nfeatures*sizeof(float));
	feature_d = dsycl::alloc<float>(q, (size_t)npoints*nfeatures);

	/* invert the data array (kernel execution); [FIX-OCL-34] local = DWARF_LS, global = roundup(npoints) */
	invert_mapping_kernel(q, DWARF_ROUNDUP(npoints), feature_flipped_d, feature_d, npoints, nfeatures);

	/* allocate memory for membership_d[] and clusters_d[][] (device) */
	membership_d = dsycl::alloc<int>(q, npoints);
	clusters_d = dsycl::alloc<float>(q, (size_t)nclusters*nfeatures);
}
/* -------------- allocateMemory() end ------------------- */

/* -------------- deallocateMemory() ------------------- */
/* free host and device memory */
	extern "C"
void deallocateMemory()
{
	free(membership_new);
	free(block_new_centers);
	dsycl::release(*Q, feature_d);
	dsycl::release(*Q, feature_flipped_d);
	dsycl::release(*Q, membership_d);
	dsycl::release(*Q, clusters_d);
}
/* -------------- deallocateMemory() end ------------------- */

/* the queue and the built kernels live until process exit (releaseCL() of the OpenCL tree, [FIX-OCL-33]) */
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
	dsycl::parse_args(&argc, &argv);
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
	sycl::queue &q = *Q;

	int delta = 0;			/* if point has moved */
	int i,j;				/* counters */

	/* copy membership (host to device) */
	dsycl::write(q, membership_d, membership_new, npoints*sizeof(int));
	/* copy clusters (host to device) */
	dsycl::write(q, clusters_d, clusters[0], (size_t)nclusters*nfeatures*sizeof(float));

	/* [FIX-OCL-34] pinned launch: 1-D, local = DWARF_LS, global = roundup(npoints); the kernel guards point_id < npoints */
	kmeansPoint_kernel(q, DWARF_ROUNDUP(npoints), feature_d, nfeatures, npoints, nclusters, membership_d, clusters_d);

	/* copy back membership (device to host) */
	dsycl::read(q, membership_new, membership_d, npoints*sizeof(int));

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

