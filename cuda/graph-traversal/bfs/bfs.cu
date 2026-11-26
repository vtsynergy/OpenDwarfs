/***********************************************************************
 * BFS.c
 *
 * CUDA translation of the OpenCL port (opencl/graph-traversal/bfs/bfs.cpp) of the Rodinia BFS test.
 * 
 * Written by: Kenneth Lee
 * Last Modified: 2/4/2010
 **********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

//#include "/home/pallavid/altera/14.0/hld/host/include/CL/cl_ext.h" //for memory banks

#include "dwarf_cuda.h"
#include "dwarf_timer.h"
static int bfs_valid = 1;

#include <ctime>

#define AOCL_ALIGNMENT 64

unsigned int no_of_nodes;
unsigned int edge_list_size;
FILE *fp;

//Structure for Nodes in the graph
//__attribute__((packed))
//__attribute__((aligned(16)))
struct Node
{
	int starting;     //Index where the edges of the node start
	int no_of_edges;  //The degree of the node
};

#include "bfs_kernel.cuh"

void BFSGraph(int argc, char** argv);

/******************************************************************************
 * MAIN METHOD
 *****************************************************************************/
int main(int argc, char** argv)
{
	dwarf_cuda_init(&argc, &argv);
	no_of_nodes = 0;
	edge_list_size = 0;
	BFSGraph(argc, argv);
	return dwarf_report(bfs_valid);
}

/******************************************************************************
 * Apply BFS on a Graph using OpenCL
 *****************************************************************************/
void BFSGraph(int argc, char ** argv)
{
	//ocd_options opts = ocd_get_options();
	//n_platform = opts.platform_id;
	//n_device = opts.device_id;
	//_deviceType = opts.device_type;

	if(argc < 2)
	{
		printf("Usage: <filename>\n");
		exit(1);
	} 
	printf("Reading File\n");
	//Read in Graph from a file
	fp = fopen(argv[1], "r");
	if(!fp)
	{
		printf("Error Reading graph file\n");
		return;
	}

	int source = 0;
	fscanf(fp, "%d", &no_of_nodes);
         
         Node* h_graph_nodes = (Node*) memalign ( AOCL_ALIGNMENT, sizeof(Node) * no_of_nodes);
         int* h_graph_mask= (int*) memalign ( AOCL_ALIGNMENT, sizeof(int) * no_of_nodes);
         int* h_updating_graph_mask= (int*) memalign (AOCL_ALIGNMENT, sizeof(int) * no_of_nodes);
         int* h_graph_visited= (int*) memalign (AOCL_ALIGNMENT, sizeof(int) * no_of_nodes);


        //allocate host memory
       
       /* Node* h_graph_nodes = (Node*) malloc(sizeof(Node) * no_of_nodes);
	int* h_graph_mask = (int*) malloc(sizeof(int) * no_of_nodes);
	int* h_updating_graph_mask = (int*) malloc(sizeof(int) * no_of_nodes);
	int* h_graph_visited = (int*) malloc(sizeof(int) * no_of_nodes);
*/

	int start, edgeno;
	//initalize the memory
	for(unsigned int i = 0; i < no_of_nodes; i++)
	{
		fscanf(fp, "%d %d", &start, &edgeno);
		h_graph_nodes[i].starting = start; 
		h_graph_nodes[i].no_of_edges = edgeno;
		h_graph_mask[i] = 0;
		h_updating_graph_mask[i] = 0;
		h_graph_visited[i] = 0;
	}

	//read the source node from the file
	fscanf(fp, "%d", &source);
	source = 0;					     //Hmmmm.... Seems that the fscanf is not really used....

	//set the source node as true in the masks
	h_graph_mask[source] = 1;
	h_graph_visited[source] = 1;

	fscanf(fp, "%d", &edge_list_size);

	int id, cost;
	//int* h_graph_edges = (int*) memalign(AOCL_ALIGNMENT, sizeof(int) * edge_list_size);
        int* h_graph_edges = (int*) malloc(sizeof(int) * edge_list_size);

	for(unsigned int i = 0; i < edge_list_size; i++)
	{
		fscanf(fp, "%d", &id);
		fscanf(fp, "%d", &cost);
		h_graph_edges[i] = id;
	}

	if(fp)
		fclose(fp);

	printf("Read File\n");

	//Allocate the graph arrays in device memory and upload them once (as the OpenCL host does before its timed loop)
	Node *d_graph_nodes; int *d_graph_edges, *d_graph_mask, *d_updating_graph_mask, *d_graph_visited, *d_cost, *d_over;
	CUDA_CHECK(cudaMalloc((void**)&d_graph_nodes, sizeof(Node) * no_of_nodes));
	CUDA_CHECK(cudaMalloc((void**)&d_graph_edges, sizeof(int) * edge_list_size));
	CUDA_CHECK(cudaMalloc((void**)&d_graph_mask, sizeof(int) * no_of_nodes));
	CUDA_CHECK(cudaMalloc((void**)&d_updating_graph_mask, sizeof(int) * no_of_nodes));
	CUDA_CHECK(cudaMalloc((void**)&d_graph_visited, sizeof(int) * no_of_nodes));
#define H2D(dst, src, bytes) CUDA_CHECK(cudaMemcpy(dst, src, bytes, cudaMemcpyHostToDevice))
	H2D(d_graph_nodes, h_graph_nodes, sizeof(Node) * no_of_nodes);
	H2D(d_graph_edges, h_graph_edges, sizeof(int) * edge_list_size);
	H2D(d_graph_mask, h_graph_mask, sizeof(int) * no_of_nodes);
	H2D(d_updating_graph_mask, h_updating_graph_mask, sizeof(int) * no_of_nodes);
	H2D(d_graph_visited, h_graph_visited, sizeof(int) * no_of_nodes);
	CUDA_CHECK(cudaDeviceSynchronize());


	int* h_cost = (int*) memalign(AOCL_ALIGNMENT, sizeof(int) * no_of_nodes);
	//int* h_cost = (int*) malloc(sizeof(int) * no_of_nodes);

	for(unsigned int i = 0; i < no_of_nodes; i++)
		h_cost[i] = -1;
	h_cost[source] = 0;
	//Allocate device memory for result
	CUDA_CHECK(cudaMalloc((void**)&d_cost, sizeof(int) * no_of_nodes));
	//Make a bool to check if the execution is over
	CUDA_CHECK(cudaMalloc((void**)&d_over, sizeof(int)));
	H2D(d_cost, h_cost, sizeof(int) * no_of_nodes);
	CUDA_CHECK(cudaDeviceSynchronize());


	printf("Copied Everything to device memory\n");

	int k = 0;
	int stop;

	/* one node per thread: DWARF_LS threads per block, ceil(nodes/DWARF_LS) blocks; both kernels guard tid < no_of_nodes */
	unsigned int nblocks = DWARF_NBLOCKS(no_of_nodes);
	printf("maxThreads[0]=%d WorkSize[0]=%d localWorkSize[0]=%d\n", DWARF_LS, (int)DWARF_ROUNDUP(no_of_nodes), DWARF_LS);
	/* common timer: the timed region is the six H2D uploads + the level loop + the final cost D2H (JIT/build above is excluded).
	 * Host arrays are never modified by the device except h_cost (overwritten by the read-back), so a pristine copy is kept. */
	dwarf_init("bfs", "graph-traversal", "cuda");
	dwarf_set_size(no_of_nodes, edge_list_size);
	int *h_cost_init = (int*) malloc(sizeof(int) * no_of_nodes);
	memcpy(h_cost_init, h_cost, sizeof(int) * no_of_nodes);
	int rep;
	for(rep = dwarf_first_rep(); rep < dwarf_reps(); rep++)
	{
	k = 0;
	dwarf_region_begin();
	H2D(d_graph_nodes, h_graph_nodes, sizeof(Node) * no_of_nodes);
	H2D(d_graph_edges, h_graph_edges, sizeof(int) * edge_list_size);
	H2D(d_graph_mask, h_graph_mask, sizeof(int) * no_of_nodes);
	H2D(d_updating_graph_mask, h_updating_graph_mask, sizeof(int) * no_of_nodes);
	H2D(d_graph_visited, h_graph_visited, sizeof(int) * no_of_nodes);
	H2D(d_cost, h_cost_init, sizeof(int) * no_of_nodes);
	do
	{
		stop = 0;
		//Copy stop to device
		H2D(d_over, &stop, sizeof(int));
		CUDA_CHECK(cudaDeviceSynchronize());
		//Run Kernel1 and Kernel2
		kernel1<<<nblocks, DWARF_LS>>>(d_graph_nodes, d_graph_edges, d_graph_mask, d_updating_graph_mask, d_graph_visited, d_cost, no_of_nodes);
		CUDA_SYNC_LAUNCH();
		kernel2<<<nblocks, DWARF_LS>>>(d_graph_mask, d_updating_graph_mask, d_graph_visited, d_over, no_of_nodes);
		CUDA_SYNC_LAUNCH();

		//Copy stop from device
		CUDA_CHECK(cudaMemcpy(&stop, d_over, sizeof(int), cudaMemcpyDeviceToHost));
		CUDA_CHECK(cudaDeviceSynchronize());
		k++;
	}while(stop == 1);

	printf("Kernel Executed %d times\n", k);

	//copy result form device to host

	CUDA_CHECK(cudaMemcpy(h_cost, d_cost, sizeof(int)*no_of_nodes, cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaDeviceSynchronize());
	dwarf_region_end(rep);
	}

	//Store the result into a file
	FILE* fpo = fopen("bfs_result.txt", "w");
	for(unsigned int i = 0; i < no_of_nodes; i++)
		fprintf(fpo, "%d) cost:%d\n", i, h_cost[i]);
	fclose(fpo);
	printf("Result stored in bfs_result.txt\n");

	/* [FIX-OCL-05] the app had no correctness check. Added: serial level-synchronous BFS from the same source,
	 * compared with the device cost array (the kernels are deterministic per level, so the costs must be identical). */
	{
		int *ref = (int*) malloc(sizeof(int) * no_of_nodes), *q = (int*) malloc(sizeof(int) * no_of_nodes);
		unsigned int qh = 0, qt = 0, bad = 0;
		for(unsigned int i = 0; i < no_of_nodes; i++) ref[i] = -1;
		ref[source] = 0; q[qt++] = source;
		while(qh < qt)
		{
			int u = q[qh++];
			for(int e = h_graph_nodes[u].starting; e < h_graph_nodes[u].starting + h_graph_nodes[u].no_of_edges; e++)
				if(ref[h_graph_edges[e]] < 0) { ref[h_graph_edges[e]] = ref[u] + 1; q[qt++] = h_graph_edges[e]; }
		}
		for(unsigned int i = 0; i < no_of_nodes; i++) if(ref[i] != h_cost[i]) bad++;
		printf("%s: %u cost mismatches vs serial BFS\n", bad ? "FAILURE" : "SUCCESS", bad);
		bfs_valid = (bad == 0);
		free(ref); free(q);
	}

	//cleanup memory
	//Free Host memory
	free(h_graph_nodes);
	free(h_graph_edges);
	free(h_graph_mask);
	free(h_updating_graph_mask);
	free(h_graph_visited);
	free(h_cost);
	//Free device memory
	CUDA_CHECK(cudaFree(d_graph_nodes));
	CUDA_CHECK(cudaFree(d_graph_edges));
	CUDA_CHECK(cudaFree(d_graph_mask));
	CUDA_CHECK(cudaFree(d_updating_graph_mask));
	CUDA_CHECK(cudaFree(d_graph_visited));
	CUDA_CHECK(cudaFree(d_cost));
	CUDA_CHECK(cudaFree(d_over));
}
