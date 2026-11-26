/* CUDA translation of bfs_kernel.cl: one node per thread, DWARF_LS threads/block; both kernels guard tid < no_of_nodes.
 * struct Node is declared in bfs.cu (same layout as the OpenCL typedef). */
__global__ void kernel1(const Node* g_graph_nodes,
		int* g_graph_edges,
		int* g_graph_mask,
		int* g_updating_graph_mask,
		int* g_graph_visited,
		int* g_cost,
		int no_of_nodes)
{
	unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;

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
}

__global__ void kernel2(int* g_graph_mask,
		int* g_updating_graph_mask,
		int* g_graph_visited,
		int* g_over,
		int no_of_nodes)
{
	unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if(tid < no_of_nodes && g_updating_graph_mask[tid] == 1)
	{
		g_graph_mask[tid] = 1;
		g_graph_visited[tid] = 1;
		*g_over = 1;
		g_updating_graph_mask[tid] = 0;
	}
}
