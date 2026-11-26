/* Breadth-first search -- OpenMP target-offload translation of opencl/graph-traversal/bfs/bfs.cpp
 * (Rodinia BFS: level-synchronous, two kernels per level with mask / updating-mask / visited arrays).
 *
 * Same graph file format, same data structures (Node{starting,no_of_edges}), same host loop:
 *   do { stop=0; copy stop to device; kernel1; kernel2; copy stop back; } while(stop)
 * Same quirks: the source node is forced to 0 and the cost field of each edge is ignored.
 * Usage: bfs <graph-file>      (writes bfs_result.txt like the OpenCL app)
 * Added (not in the OpenCL app): validation of the cost array against a serial BFS.
 * Timed region: all H2D copies, the level loop, the final cost D2H.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "dwarf_timer.h"
#include "dwarf_launch.h" /* [FIX-OMP-13] pinned launch configuration */

struct Node {
  int starting;     // index where the edges of the node start
  int no_of_edges;  // the degree of the node
};

static unsigned int no_of_nodes, edge_list_size;

/* kernel1 of bfs_kernel.cl: expand the frontier */
static void kernel1(const Node *nodes, const int *edges, int *mask, int *updating, const int *visited, int *cost,
                    int n) {
/* [FIX-OMP-13] was thread_limit(256) only; pinned: num_teams=ceil(n/256), thread_limit=num_threads=256 (same as the OpenCL NDRange) */
#pragma omp target teams distribute parallel for DWARF_LAUNCH(n)
  for (int tid = 0; tid < n; tid++) {
    if (mask[tid] != 0) {
      mask[tid] = 0;
      int max = nodes[tid].no_of_edges + nodes[tid].starting;
      for (int i = nodes[tid].starting; i < max; i++) {
        int id = edges[i];
        if (!visited[id]) {
          cost[id] = cost[tid] + 1;
          updating[id] = 1;
        }
      }
    }
  }
}

/* kernel2 of bfs_kernel.cl: promote updated nodes to the next frontier, raise the "over" flag */
static void kernel2(int *mask, int *updating, int *visited, int *over, int n) {
/* [FIX-OMP-13] was thread_limit(256) only; pinned: num_teams=ceil(n/256), thread_limit=num_threads=256 (same as the OpenCL NDRange) */
#pragma omp target teams distribute parallel for DWARF_LAUNCH(n)
  for (int tid = 0; tid < n; tid++) {
    if (updating[tid] == 1) {
      mask[tid] = 1;
      visited[tid] = 1;
#pragma omp atomic write
      *over = 1;
      updating[tid] = 0;
    }
  }
}

int main(int argc, char **argv) {
  if (argc < 2) { printf("Usage: <filename>\n"); exit(1); }
  printf("Reading File\n");
  FILE *fp = fopen(argv[1], "r");
  if (!fp) { printf("Error Reading graph file\n"); return 1; }

  int source = 0;
  if (fscanf(fp, "%u", &no_of_nodes) != 1) return 1;
  Node *h_graph_nodes = (Node *)malloc(sizeof(Node) * no_of_nodes);
  int *h_graph_mask = (int *)malloc(sizeof(int) * no_of_nodes);
  int *h_updating_graph_mask = (int *)malloc(sizeof(int) * no_of_nodes);
  int *h_graph_visited = (int *)malloc(sizeof(int) * no_of_nodes);
  int *h_cost = (int *)malloc(sizeof(int) * no_of_nodes);
  int start, edgeno;
  for (unsigned int i = 0; i < no_of_nodes; i++) {
    if (fscanf(fp, "%d %d", &start, &edgeno) != 2) return 1;
    h_graph_nodes[i].starting = start;
    h_graph_nodes[i].no_of_edges = edgeno;
  }
  if (fscanf(fp, "%d", &source) != 1) return 1;
  source = 0;  // as in the OpenCL app: the value read from the file is not used
  if (fscanf(fp, "%u", &edge_list_size) != 1) return 1;
  int *h_graph_edges = (int *)malloc(sizeof(int) * edge_list_size);
  int id, cost_unused;
  for (unsigned int i = 0; i < edge_list_size; i++) {
    if (fscanf(fp, "%d", &id) != 1 || fscanf(fp, "%d", &cost_unused) != 1) return 1;
    h_graph_edges[i] = id;
  }
  fclose(fp);
  printf("Read File\n");

  const int n = (int)no_of_nodes;
  const unsigned int E = edge_list_size;
  int stop = 0, levels = 0;
  Node *nodes = h_graph_nodes;
  int *edges = h_graph_edges, *mask = h_graph_mask, *updating = h_updating_graph_mask, *visited = h_graph_visited,
      *cost = h_cost;

  dwarf_init("bfs", "graph-traversal", "omp");
  dwarf_set_size(n, E);
  /* [FIX-OMP-14] OpenCL creates its buffers before the timed region and only transfers inside it: allocate once outside the region, `target update` inside */
#pragma omp target enter data map(alloc : nodes[0 : n], edges[0 : E], mask[0 : n], updating[0 : n], visited[0 : n], cost[0 : n], stop)
  for (int rep = dwarf_first_rep(); rep < dwarf_reps(); ++rep) {
    /* untimed host initialisation (identical to the OpenCL host code) */
    for (int i = 0; i < n; i++) { mask[i] = 0; updating[i] = 0; visited[i] = 0; cost[i] = -1; }
    mask[source] = 1; visited[source] = 1; cost[source] = 0;
    levels = 0;

    dwarf_region_begin();
#pragma omp target update to(nodes[0 : n], edges[0 : E], mask[0 : n], updating[0 : n], visited[0 : n], cost[0 : n])
    do {
      stop = 0;
#pragma omp target update to(stop)  /* copy stop to device */
      kernel1(nodes, edges, mask, updating, visited, cost, n);
      kernel2(mask, updating, visited, &stop, n);
#pragma omp target update from(stop) /* copy stop from device */
      levels++;
    } while (stop == 1);
#pragma omp target update from(cost[0 : n])
    dwarf_region_end(rep);
  }
#pragma omp target exit data map(delete : nodes[0 : n], edges[0 : E], mask[0 : n], updating[0 : n], visited[0 : n], cost[0 : n], stop)
  printf("Kernel Executed %d times\n", levels);

  FILE *fpo = fopen("bfs_result.txt", "w");
  for (int i = 0; i < n; i++) fprintf(fpo, "%d) cost:%d\n", i, cost[i]);
  fclose(fpo);
  printf("Result stored in bfs_result.txt\n");

  /* validation: serial level-synchronous BFS from the same source */
  int *ref = (int *)malloc(sizeof(int) * n), *q = (int *)malloc(sizeof(int) * n);
  for (int i = 0; i < n; i++) ref[i] = -1;
  int qh = 0, qt = 0; ref[source] = 0; q[qt++] = source;
  while (qh < qt) {
    int u = q[qh++];
    for (int e = nodes[u].starting; e < nodes[u].starting + nodes[u].no_of_edges; e++)
      if (ref[edges[e]] < 0) { ref[edges[e]] = ref[u] + 1; q[qt++] = edges[e]; }
  }
  int bad = 0;
  for (int i = 0; i < n; i++) if (ref[i] != cost[i]) bad++;
  printf("%s: %d cost mismatches vs serial BFS\n", bad ? "FAILURE" : "SUCCESS", bad);
  return dwarf_report(bad == 0);
}
