#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "../../include/rdtsc.h"
#include "../../include/common_args.h"

#include <omp.h>

#define THREADS_PER_DIM 16
#define BLOCKS_PER_DIM 16


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

cl_program clProgram;
cl_kernel clKernel_invert_mapping;
cl_kernel clKernel_kmeansPoint;

/* _d denotes it resides on the device */
int    *membership_new;												/* newly assignment membership */
cl_mem feature_d;													/* inverted data array */
cl_mem feature_flipped_d;											/* original (not inverted) data array */
cl_mem membership_d;												/* membership on the device */
float  *block_new_centers;											/* sum of points in a cluster (per block) */
cl_mem clusters_d;													/* cluster centers on the device */
cl_mem block_clusters_d;											/* per block calculation of cluster centers */
cl_mem block_deltas_d;												/* per block calculation of deltas */

/* image memory */
/*cl_mem t_features;
  cl_mem t_features_flipped;
  cl_mem t_clusters;*/


/* -------------- initCL() -------------- */
	extern "C"
void initCL()
{
	FILE *kernelFile;
	char *kernelSource;
	size_t kernelLength;

	cl_int errcode;

	ocd_initCL();

	size_t max_worksize[3];
	errcode = clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_ITEM_SIZES,sizeof(size_t)*3, &max_worksize, NULL);
	CHKERR(errcode, "Failed to get device info!");
	while(num_threads_perdim*num_threads_perdim>max_worksize[0])
		num_threads_perdim = num_threads_perdim/2;

	num_threads = num_threads_perdim*num_threads_perdim;



	kernelFile = fopen("kmeans_opencl_kernel.cl", "r");
	fseek(kernelFile, 0, SEEK_END);
	kernelLength = (size_t) ftell(kernelFile);
	kernelSource = (char *) malloc(sizeof(char)*kernelLength);
	rewind(kernelFile);
	fread((void *) kernelSource, kernelLength, 1, kernelFile);
	fclose(kernelFile);

	clProgram = clCreateProgramWithSource(context, 1, (const char **) &kernelSource, &kernelLength, &errcode);
	CHKERR(errcode, "Failed to create program with source!");

	free(kernelSource);

	errcode = clBuildProgram(clProgram, 1, &device_id, NULL, NULL, NULL);
	if (errcode == CL_BUILD_PROGRAM_FAILURE)
	{
		char *log;
		size_t logLength;
		errcode = clGetProgramBuildInfo(clProgram, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &logLength);
		log = (char *) malloc(sizeof(char)*logLength);
		errcode = clGetProgramBuildInfo(clProgram, device_id, CL_PROGRAM_BUILD_LOG, logLength, (void *) log, NULL);
		fprintf(stderr, "Kernel build error! Log:\n%s", log);
		free(log);
		return;
	}
	CHKERR(errcode, "Failed to build program!");

	clKernel_invert_mapping = clCreateKernel(clProgram, "invert_mapping", &errcode);
	CHKERR(errcode, "Failed to create kernel!");
	clKernel_kmeansPoint = clCreateKernel(clProgram, "kmeansPoint", &errcode);
	CHKERR(errcode, "Failed to create kernel!");
}
/* -------------- initCL() end -------------- */

/* -------------- allocateMemory() ------------------- */
/* allocate device memory, calculate number of blocks and threads, and invert the data array */
	extern "C"
void allocateMemory(int npoints, int nfeatures, int nclusters, float **features)
{	
	cl_int errcode;
	size_t globalWorkSize;
	size_t localWorkSize;

	num_blocks = npoints / num_threads;
	if (npoints % num_threads > 0)		/* defeat truncation */
		num_blocks++;

	num_blocks_perdim = sqrt((double) num_blocks);
	while (num_blocks_perdim * num_blocks_perdim < num_blocks)	// defeat truncation (should run once)
		num_blocks_perdim++;

	num_blocks = num_blocks_perdim*num_blocks_perdim;

	/* allocate memory for memory_new[] and initialize to -1 (host) */
	membership_new = (int*) malloc(npoints * sizeof(int));
	for(int i=0;i<npoints;i++) {
		membership_new[i] = -1;
	}

	/* allocate memory for block_new_centers[] (host) */
	block_new_centers = (float *) malloc(nclusters*nfeatures*sizeof(float));

	/* allocate memory for feature_flipped_d[][], feature_d[][] (device) */
	feature_flipped_d = clCreateBuffer(context, CL_MEM_READ_ONLY, npoints*nfeatures*sizeof(float), NULL, &errcode);
	CHKERR(errcode, "Failed to create buffer!");

	errcode = clEnqueueWriteBuffer(commands, feature_flipped_d, CL_TRUE, 0, npoints*nfeatures*sizeof(float), features[0], 0, NULL, &ocdTempEvent);

	clFinish(commands);
	START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "Point/Feature Copy", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHKERR(errcode, "Failed to create buffer!");
	feature_d = clCreateBuffer(context, CL_MEM_READ_WRITE, npoints*nfeatures*sizeof(float), NULL, &errcode);
	CHKERR(errcode, "Failed to create buffer!");

	/* invert the data array (kernel execution) */	
	unsigned int arg = 0;
	errcode = clSetKernelArg(clKernel_invert_mapping, arg++, sizeof(cl_mem), (void *) &feature_flipped_d);
	errcode |= clSetKernelArg(clKernel_invert_mapping, arg++, sizeof(cl_mem), (void *) &feature_d);
	errcode |= clSetKernelArg(clKernel_invert_mapping, arg++, sizeof(int), (void *) &npoints);
	errcode |= clSetKernelArg(clKernel_invert_mapping, arg++, sizeof(int), (void *) &nfeatures);
	CHKERR(errcode, "Failed to set kernel arg!");
	globalWorkSize = num_blocks*num_threads;
	localWorkSize = num_threads;

	errcode = clEnqueueNDRangeKernel(commands, clKernel_invert_mapping, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &ocdTempEvent);
	clFinish(commands);
	START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "Invert Mapping Kernel", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHKERR(errcode, "Failed to enqueue kernel!");

	/* allocate memory for membership_d[] and clusters_d[][] (device) */
	membership_d = clCreateBuffer(context, CL_MEM_READ_WRITE, npoints*sizeof(int), NULL, &errcode);
	CHKERR(errcode, "Failed to create buffer!");
	clusters_d = clCreateBuffer(context, CL_MEM_READ_ONLY, nclusters*nfeatures*sizeof(float), NULL, &errcode);
	CHKERR(errcode, "Failed to create buffer!");

#ifdef BLOCK_DELTA_REDUCE
	// allocate array to hold the per block deltas on the gpu side

	block_deltas_d = clCreateBuffer(context, CL_MEM_READ_WRITE, num_blocks_perdim*num_blocks_perdim*sizeof(int), NULL, &errcode);
	CHKERR(errcode, "Failed to create buffer!");
	//cudaMemcpy(block_delta_d, &delta_h, sizeof(int), cudaMemcpyHostToDevice);
#endif

#ifdef BLOCK_CENTER_REDUCE
	// allocate memory and copy to card cluster  array in which to accumulate center points for the next iteration
	block_clusters_d = clCreateBuffer(context, CL_MEM_READ_WRITE, num_blocks_perdim*num_blocks_perdim*nclusters*nfeatures*sizeof(float), NULL, &errcode);
	CHKERR(errcode, "Failed to create buffer!");
	//cudaMemcpy(new_clusters_d, new_centers[0], nclusters*nfeatures*sizeof(float), cudaMemcpyHostToDevice);
#endif

}
/* -------------- allocateMemory() end ------------------- */

/* -------------- deallocateMemory() ------------------- */
/* free host and device memory */
	extern "C"
void deallocateMemory()
{
	free(membership_new);
	free(block_new_centers);
	clReleaseMemObject(feature_d);
	clReleaseMemObject(feature_flipped_d);
	clReleaseMemObject(membership_d);

	clReleaseMemObject(clusters_d);
#ifdef BLOCK_CENTER_REDUCE
	clReleaseMemObject(block_clusters_d);
#endif
#ifdef BLOCK_DELTA_REDUCE
	clReleaseMemObject(block_deltas_d);
#endif
	clReleaseKernel(clKernel_invert_mapping);
	clReleaseKernel(clKernel_kmeansPoint);
	clReleaseProgram(clProgram);
	clReleaseCommandQueue(commands);
	clReleaseContext(context);
}
/* -------------- deallocateMemory() end ------------------- */



////////////////////////////////////////////////////////////////////////////////
// Program main																  //

	int
main( int argc, char** argv) 
{
	// as done in the CUDA start/help document provided
	ocd_init(&argc, &argv, NULL);
	setup(argc, argv);   
	ocd_finalize();
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
	cl_int errcode;

	int delta = 0;			/* if point has moved */
	int i,j;				/* counters */


	/* copy membership (host to device) */

	errcode = clEnqueueWriteBuffer(commands, membership_d, CL_TRUE, 0, npoints*sizeof(int), (void *) membership_new, 0, NULL, &ocdTempEvent);
	clFinish(commands);
	START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "Membership Copy", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHKERR(errcode, "Failed to enqueue write buffer!");

	/* copy clusters (host to device) */

	errcode = clEnqueueWriteBuffer(commands, clusters_d, CL_TRUE, 0, nclusters*nfeatures*sizeof(float), (void *) clusters[0], 0, NULL, &ocdTempEvent);
	clFinish(commands);
	START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "Cluster Copy", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHKERR(errcode, "Failed to enqueue write buffer!");

	/* set up texture */
	/*cudaChannelFormatDesc chDesc0 = cudaCreateChannelDesc<float>();
	  t_features.filterMode = cudaFilterModePoint;   
	  t_features.normalized = false;
	  t_features.channelDesc = chDesc0;

	  if(cudaBindTexture(NULL, &t_features, feature_d, &chDesc0, npoints*nfeatures*sizeof(float)) != CUDA_SUCCESS)
	  printf("Couldn't bind features array to texture!\n");

	  cudaChannelFormatDesc chDesc1 = cudaCreateChannelDesc<float>();
	  t_features_flipped.filterMode = cudaFilterModePoint;   
	  t_features_flipped.normalized = false;
	  t_features_flipped.channelDesc = chDesc1;

	  if(cudaBindTexture(NULL, &t_features_flipped, feature_flipped_d, &chDesc1, npoints*nfeatures*sizeof(float)) != CUDA_SUCCESS)
	  printf("Couldn't bind features_flipped array to texture!\n");

	  cudaChannelFormatDesc chDesc2 = cudaCreateChannelDesc<float>();
	  t_clusters.filterMode = cudaFilterModePoint;   
	  t_clusters.normalized = false;
	  t_clusters.channelDesc = chDesc2;

	  if(cudaBindTexture(NULL, &t_clusters, clusters_d, &chDesc2, nclusters*nfeatures*sizeof(float)) != CUDA_SUCCESS)
	  printf("Couldn't bind clusters array to texture!\n");*/

	/* copy clusters to constant memory */
	//cudaMemcpyToSymbol("c_clusters",clusters[0],nclusters*nfeatures*sizeof(float),0,cudaMemcpyHostToDevice);


	/* setup execution parameters.
	   changed to 2d (source code on NVIDIA CUDA Programming Guide) */
	size_t localWorkSize[2] = {num_threads_perdim*num_threads_perdim, 1};
	size_t globalWorkSize[2] = {num_blocks_perdim*localWorkSize[0], num_blocks_perdim*localWorkSize[1]};

	unsigned int arg = 0;
	errcode = clSetKernelArg(clKernel_kmeansPoint, arg++, sizeof(cl_mem), (void *) &feature_d);
	errcode |= clSetKernelArg(clKernel_kmeansPoint, arg++, sizeof(cl_mem), (void *) &feature_flipped_d);
	errcode |= clSetKernelArg(clKernel_kmeansPoint, arg++, sizeof(int), (void *) &nfeatures);
	errcode |= clSetKernelArg(clKernel_kmeansPoint, arg++, sizeof(int), (void *) &npoints);
	errcode |= clSetKernelArg(clKernel_kmeansPoint, arg++, sizeof(int), (void *) &nclusters);
	errcode |= clSetKernelArg(clKernel_kmeansPoint, arg++, sizeof(cl_mem), (void *) &membership_d);
	errcode |= clSetKernelArg(clKernel_kmeansPoint, arg++, sizeof(cl_mem), (void *) &clusters_d);
#ifdef BLOCK_DELTA_REDUCE
	errcode |= clSetKernelArg(clKernel_kmeansPoint, arg++, sizeof(cl_mem), (void *) &block_clusters_d);
#endif
#ifdef BLOCK_CENTER_REDUCE
	errcode |= clSetKernelArg(clKernel_kmeansPoint, arg++, sizeof(cl_mem), (void *) &block_deltas_d);
#endif
	CHKERR(errcode, "Failed to set kernel arg!");

	/* execute the kernel */

	errcode = clEnqueueNDRangeKernel(commands, clKernel_kmeansPoint, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
	CHKERR(errcode, "Failed to enqueue kernel!");
	errcode = clFinish(commands);
	START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "Point Kernel", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHKERR(errcode, "Failed to clFinish!");
	/* copy back membership (device to host) */

	errcode = clEnqueueReadBuffer(commands, membership_d, CL_TRUE, 0, npoints*sizeof(int), (void *) membership_new, 0, NULL, &ocdTempEvent);
	clFinish(commands);
	START_TIMER(ocdTempEvent, OCD_TIMER_D2H, "Membership Copy", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHKERR(errcode, "Failed to enqueue read buffer!");

#ifdef BLOCK_CENTER_REDUCE
	/*** Copy back arrays of per block sums ***/
	float * block_clusters_h = (float *) malloc(
			num_blocks_perdim * num_blocks_perdim * 
			nclusters * nfeatures * sizeof(float));

	errcode = clEnqueueReadBuffer(commands, block_clusters_d, CL_TRUE, 0, num_blocks_perdim*num_blocks_perdim*nclusters*nfeatures*sizeof(float), (void *) block_clusters_h, 0, NULL, &ocdTempEvent);
	clFinish(commands);
	START_TIMER(ocdTempEvent, OCD_TIMER_D2H, "Block_clusters Copy", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHKERR(errcode, "Failed to enqueue read buffer!");
#endif
#ifdef BLOCK_DELTA_REDUCE
	int * block_deltas_h = (int *) malloc(
			num_blocks_perdim * num_blocks_perdim * sizeof(int));

	errcode = clEnqueueReadBuffer(commands, block_deltas_d, CL_TRUE, 0, num_blocks_perdim*num_blocks_perdim*sizeof(int), (void *) block_deltas_h, 0, NULL, &ocdTempEvent);
	clFinish(commands);
	START_TIMER(ocdTempEvent, OCD_TIMER_D2H, "Block_deltas Copy", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHKERR(errcode, "Failed to enqueue read buffer!");
#endif

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

