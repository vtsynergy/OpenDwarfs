/* CUDA translation of kmeans_opencl_kernel.cl (GPU_DELTA_REDUCTION / GPU_NEW_CENTER_REDUCTION are undefined there, so only the
 * code actually compiled is translated). One point per thread, DWARF_LS threads/block, guards point_id < npoints. The OpenCL
 * __constant buffer `clusters` is ordinary device memory here. */
#include <float.h>

__global__ void invert_mapping(float *input, float *output, int npoints, int nfeatures)
{
	int point_id = threadIdx.x + blockDim.x*blockIdx.x;	/* id of thread */
	int i;

	if(point_id < npoints){
		for(i=0;i<nfeatures;i++)
			output[point_id + npoints*i] = input[point_id*nfeatures + i];
	}
	return;
}

__global__ void kmeansPoint(float *features, float *features_flipped, int nfeatures, int npoints, int nclusters,
		int *membership, const float *clusters)
{
	const unsigned int block_id = gridDim.x*blockIdx.y+blockIdx.x;
	const unsigned int point_id = block_id*blockDim.x*blockDim.y + threadIdx.x;

	int  index = -1;

	if (point_id < npoints)
	{
		int i, j;
		float min_dist = FLT_MAX;
		float dist;													/* distance square between a point to cluster center */

		/* find the cluster center id with min distance to pt */
		for (i=0; i<nclusters; i++) {
			int cluster_base_index = i*nfeatures;					/* base index of cluster centers for inverted array */
			float ans=0.0f;												/* Euclidean distance sqaure */

			for (j=0; j < nfeatures; j++)
			{
				int addr = point_id + j*npoints;					/* appropriate index of data point */
				float diff = (features[addr] -
						clusters[cluster_base_index + j]);	/* distance between a data point to cluster centers */
				ans += diff*diff;									/* sum of squares */
			}
			dist = ans;

			if (dist < min_dist) {
				min_dist = dist;
				index    = i;
			}
		}
	}

	if (point_id < npoints)
	{
		/* assign the membership to object point_id */
		membership[point_id] = index;
	}
}
