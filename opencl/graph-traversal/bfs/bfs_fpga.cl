
#define UNROLL_FACTOR 8	

typedef struct
{
	int starting;
	int no_of_edges;
}Node;

__kernel void kernel1(__global const Node* restrict  g_graph_nodes,
		__global int* restrict  g_graph_edges,
		__global int* restrict  g_graph_mask,
		__global int* restrict  g_updating_graph_mask,
		__global int* restrict  g_graph_visited,
		__global int* restrict  g_cost,
		int no_of_nodes) 
{
	
	__private size_t tid= get_global_id(0);

	if(tid < no_of_nodes && g_graph_mask[tid] != 0)
	{
		g_graph_mask[tid] = 0;
		__private const unsigned int start = g_graph_nodes[tid].starting;
		__private const unsigned int end = g_graph_nodes[tid].no_of_edges + g_graph_nodes[tid].starting;
		__private unsigned int id[UNROLL_FACTOR];
		__private int Mask[UNROLL_FACTOR];
		
                __private const int cost= g_cost[tid]+1;
		__private int id_,Mask_;
		
              
		for (int i = start; i < end; i=i+UNROLL_FACTOR)
		{

		   if(i+UNROLL_FACTOR-1 < end)
	           {

		   #pragma unroll UNROLL_FACTOR
			for(int ii=0; ii<UNROLL_FACTOR; ii++) 
			{
				id[ii] = g_graph_edges[i+ii];				
               			Mask[ii] = -(!(g_graph_visited[id[ii]])); // (!0)=1  -0=0 -1=all ones 
               			g_cost[id[ii]] =  ((cost)&(Mask[ii]))|((g_cost[id[ii]])&~(Mask[ii]));
				g_updating_graph_mask[id[ii]] = 1 & !(g_graph_visited[id[ii]]);
				
			 }	
		 }
		 else
		 {
		 	for(int ii=i; ii<end; ii++) 
			{
				 id_ = g_graph_edges[ii];
				 Mask_ = -(!(g_graph_visited[id_]));
				 g_cost[id_] =  ((cost)&Mask_)|((g_cost[id_])&~Mask_);
				 g_updating_graph_mask[id_] = 1 & !(g_graph_visited[id_]);
			}
		 }
	    }
	}
}

__kernel void kernel2(__global int* restrict  g_graph_mask,
		__global int* restrict  g_updating_graph_mask,
		__global int* restrict  g_graph_visited,
		__global int* restrict  g_over,
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
