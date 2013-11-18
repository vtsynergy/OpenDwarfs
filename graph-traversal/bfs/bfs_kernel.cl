typedef struct
{
	int starting;
	int no_of_edges;
}Node;

__kernel void kernel1(__global const Node* g_graph_nodes,
		__global int* g_graph_edges,
		__global int* g_graph_mask,
		__global int* g_updating_graph_mask,
		__global int* g_graph_visited,
		__global int* g_cost,
		int no_of_nodes) 
{
	unsigned int tid = get_global_id(0);

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

__kernel void kernel2(__global int* g_graph_mask,
		__global int* g_updating_graph_mask,
		__global int* g_graph_visited,
		__global int* g_over,
		int no_of_nodes)
{
	unsigned int tid = get_global_id(0);
	if(tid < no_of_nodes && g_updating_graph_mask[tid] == 1)
	{
		g_graph_mask[tid] = 1;
		g_graph_visited[tid] = 1;
		*g_over = 1;
		g_updating_graph_mask[tid] = 0;
	}	
}
