/***********************************************************************
 * BFS.c
 *
 * A port of the Rodinia BFS test from CUDA to OpenCL.
 * 
 * Written by: Kenneth Lee
 * Last Modified: 2/4/2010
 **********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include "dwarf_sycl.hpp"
#include "dwarf_timer.h"
static int bfs_valid = 1;

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


/* kernels (translated from bfs_kernel.cl) */
static void kernel1(sycl::queue &q, size_t global, const Node* g_graph_nodes, int* g_graph_edges, int* g_graph_mask,
		int* g_updating_graph_mask, int* g_graph_visited, int* g_cost, int no_of_nodes)
{
	q.parallel_for<class bfs_kernel1>(sycl::nd_range<1>(global, DWARF_LS), [=](sycl::nd_item<1> it) {
		unsigned int tid = it.get_global_id(0);
		if(tid < no_of_nodes && g_graph_mask[tid] != 0)
		{
			g_graph_mask[tid] = 0;
			int max = (g_graph_nodes[tid].no_of_edges + g_graph_nodes[tid].starting);
			for(int i = g_graph_nodes[tid].starting; i < max; i++)
			{
				int id = g_graph_edges[i];
				if(!g_graph_visited[id])
				{
					g_cost[id] = g_cost[tid] + 1;
					g_updating_graph_mask[id] = 1;
				}
			}
		}
	}).wait_and_throw();
}
static void kernel2(sycl::queue &q, size_t global, int* g_graph_mask, int* g_updating_graph_mask, int* g_graph_visited, int* g_over, int no_of_nodes)
{
	q.parallel_for<class bfs_kernel2>(sycl::nd_range<1>(global, DWARF_LS), [=](sycl::nd_item<1> it) {
		unsigned int tid = it.get_global_id(0);
		if(tid < no_of_nodes && g_updating_graph_mask[tid] == 1)
		{
			g_graph_mask[tid] = 1;
			g_graph_visited[tid] = 1;
			*g_over = 1;
			g_updating_graph_mask[tid] = 0;
		}
	}).wait_and_throw();
}

void BFSGraph(int argc, char** argv);

/******************************************************************************
 * MAIN METHOD
 *****************************************************************************/
int main(int argc, char** argv)
{
	dsycl::parse_args(&argc, &argv);
	no_of_nodes = 0;
	edge_list_size = 0;
	BFSGraph(argc, argv);
	return dwarf_report(bfs_valid);
}

/******************************************************************************
 * Apply BFS on a Graph using SYCL
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

	//Copy the Node list to device memory
	sycl::queue q = dsycl::make_queue();
	Node *d_graph_nodes = dsycl::alloc<Node>(q, no_of_nodes);
	int *d_graph_edges = dsycl::alloc<int>(q, edge_list_size);
	int *d_graph_mask = dsycl::alloc<int>(q, no_of_nodes);
	int *d_updating_graph_mask = dsycl::alloc<int>(q, no_of_nodes);
	int *d_graph_visited = dsycl::alloc<int>(q, no_of_nodes);

	int* h_cost = (int*) memalign(AOCL_ALIGNMENT, sizeof(int) * no_of_nodes);
	for(unsigned int i = 0; i < no_of_nodes; i++)
		h_cost[i] = -1;
	h_cost[source] = 0;
	//Allocate device memory for result
	int *d_cost = dsycl::alloc<int>(q, no_of_nodes);
	//Make a bool to check if the execution is over
	int *d_over = dsycl::alloc<int>(q, 1);

	/* untimed first uploads, as in the OpenCL tree (bfs.cpp:180-213): the device pages are touched before the timed region */
	dsycl::write(q, d_graph_nodes, h_graph_nodes, sizeof(Node) * no_of_nodes);
	dsycl::write(q, d_graph_edges, h_graph_edges, sizeof(int) * edge_list_size);
	dsycl::write(q, d_graph_mask, h_graph_mask, sizeof(int) * no_of_nodes);
	dsycl::write(q, d_updating_graph_mask, h_updating_graph_mask, sizeof(int) * no_of_nodes);
	dsycl::write(q, d_graph_visited, h_graph_visited, sizeof(int) * no_of_nodes);
	dsycl::write(q, d_cost, h_cost, sizeof(int) * no_of_nodes);

	dsycl::prebuild(q);   /* kernel JIT, as clBuildProgram: outside the timed region */

	int k = 0;
	int stop;
	size_t WorkSize = DWARF_ROUNDUP(no_of_nodes);
	printf("maxThreads[0]=%d WorkSize[0]=%d localWorkSize[0]=%d\n", DWARF_LS, (int)WorkSize, DWARF_LS);

	/* common timer: the timed region is the six H2D uploads + the level loop + the final cost D2H (JIT/build above is excluded).
	 * Host arrays are never modified by the device except h_cost (overwritten by the read-back), so a pristine copy is kept. */
	dwarf_init("bfs", "graph-traversal", "sycl");
	dwarf_set_size(no_of_nodes, edge_list_size);
	int *h_cost_init = (int*) malloc(sizeof(int) * no_of_nodes);
	memcpy(h_cost_init, h_cost, sizeof(int) * no_of_nodes);
	int rep;
	for(rep = dwarf_first_rep(); rep < dwarf_reps(); rep++)
	{
	k = 0;
	dwarf_region_begin();
	dsycl::write(q, d_graph_nodes, h_graph_nodes, sizeof(Node) * no_of_nodes);
	dsycl::write(q, d_graph_edges, h_graph_edges, sizeof(int) * edge_list_size);
	dsycl::write(q, d_graph_mask, h_graph_mask, sizeof(int) * no_of_nodes);
	dsycl::write(q, d_updating_graph_mask, h_updating_graph_mask, sizeof(int) * no_of_nodes);
	dsycl::write(q, d_graph_visited, h_graph_visited, sizeof(int) * no_of_nodes);
	dsycl::write(q, d_cost, h_cost_init, sizeof(int) * no_of_nodes);
	do
	{
		stop = 0;
		//Copy stop to device
		dsycl::write(q, d_over, &stop, sizeof(int));
		//Run Kernel1 and Kernel2
		kernel1(q, WorkSize, d_graph_nodes, d_graph_edges, d_graph_mask, d_updating_graph_mask, d_graph_visited, d_cost, no_of_nodes);
		kernel2(q, WorkSize, d_graph_mask, d_updating_graph_mask, d_graph_visited, d_over, no_of_nodes);
		//Copy stop from device
		dsycl::read(q, &stop, d_over, sizeof(int));
		k++;
	}while(stop == 1);

	printf("Kernel Executed %d times\n", k);

	//copy result form device to host
	dsycl::read(q, h_cost, d_cost, sizeof(int)*no_of_nodes);
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
	dsycl::release(q, d_graph_nodes); dsycl::release(q, d_graph_edges); dsycl::release(q, d_graph_mask);
	dsycl::release(q, d_updating_graph_mask); dsycl::release(q, d_graph_visited); dsycl::release(q, d_cost); dsycl::release(q, d_over);
}
