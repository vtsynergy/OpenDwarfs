/* Breadth-first search -- OpenACC translation of opencl/graph-traversal/bfs/{bfs.cpp,bfs_kernel.cl}
 * (Rodinia BFS: level-synchronous, two kernels per level with mask / updating-mask / visited arrays).
 *
 * Same graph file format, same data structures (Node{starting,no_of_edges}), same host loop:
 *   do { stop=0; write stop; kernel1; kernel2; read stop; } while(stop)
 * Same quirks: the source node is forced to 0 and the cost field of each edge is ignored.
 * Usage: bfs <graph-file>      (writes bfs_result.txt like the OpenCL app)
 * Validation as in the OpenCL app ([FIX-OCL-05]): cost array against a serial BFS.
 * Device buffers are created outside the timed region; the timed region is the six H2D copies,
 * the level loop and the final cost D2H (as in the OpenCL host).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <openacc.h>
#include "dwarf_timer.h"
#include "dwarf_acc.h" /* [FIX-OCL-34] pinned launch configuration */

struct Node {
  int starting;     // index where the edges of the node start
  int no_of_edges;  // the degree of the node
};

static unsigned int no_of_nodes, edge_list_size;

/* kernel1 of bfs_kernel.cl: expand the frontier */
static void kernel1(const Node *nodes, const int *edges, int *mask, int *updating, const int *visited, int *cost,
                    int n, int E) {
#pragma acc parallel loop gang vector ACC_LAUNCH(n) present(nodes[0 : n], edges[0 : E], mask[0 : n], updating[0 : n], visited[0 : n], cost[0 : n])
  for (int tid = 0; tid < n; tid++) {
    if (mask[tid] != 0) {
      mask[tid] = 0;
      int max = nodes[tid].no_of_edges + nodes[tid].starting;
#pragma acc loop seq
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
#pragma acc parallel loop gang vector ACC_LAUNCH(n) present(mask[0 : n], updating[0 : n], visited[0 : n], over[0 : 1])
  for (int tid = 0; tid < n; tid++) {
    if (updating[tid] == 1) {
      mask[tid] = 1;
      visited[tid] = 1;
#pragma acc atomic write
      over[0] = 1;
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
  const int E = (int)edge_list_size;
  int levels = 0;
  Node *nodes = h_graph_nodes;
  int *edges = h_graph_edges, *mask = h_graph_mask, *updating = h_updating_graph_mask, *visited = h_graph_visited,
      *cost = h_cost;
  static int over[1];

  acc_init(acc_get_device_type()); /* device initialisation, outside the timed region */
#pragma acc enter data create(nodes[0 : n], edges[0 : E], mask[0 : n], updating[0 : n], visited[0 : n], cost[0 : n], over[0 : 1])

  dwarf_init("bfs", "graph-traversal", "openacc");
  dwarf_set_size(n, E);
  for (int rep = dwarf_first_rep(); rep < dwarf_reps(); ++rep) {
    /* untimed host initialisation (identical to the OpenCL host code) */
    for (int i = 0; i < n; i++) { mask[i] = 0; updating[i] = 0; visited[i] = 0; cost[i] = -1; }
    mask[source] = 1; visited[source] = 1; cost[source] = 0;
    levels = 0;

    dwarf_region_begin();
#pragma acc update device(nodes[0 : n], edges[0 : E], mask[0 : n], updating[0 : n], visited[0 : n], cost[0 : n])
    do {
      over[0] = 0;
#pragma acc update device(over[0 : 1])
      kernel1(nodes, edges, mask, updating, visited, cost, n, E);
      kernel2(mask, updating, visited, over, n);
#pragma acc update self(over[0 : 1])
      levels++;
    } while (over[0] == 1);
#pragma acc update self(cost[0 : n])
    dwarf_region_end(rep);
    printf("Kernel Executed %d times\n", levels);
  }
#pragma acc exit data delete(nodes[0 : n], edges[0 : E], mask[0 : n], updating[0 : n], visited[0 : n], cost[0 : n], over[0 : 1])

  FILE *fpo = fopen("bfs_result.txt", "w");
  for (int i = 0; i < n; i++) fprintf(fpo, "%d) cost:%d\n", i, cost[i]);
  fclose(fpo);
  printf("Result stored in bfs_result.txt\n");

  /* validation: serial level-synchronous BFS from the same source */
  int *ref = (int *)malloc(sizeof(int) * n), *q = (int *)malloc(sizeof(int) * n);
  for (int i = 0; i < n; i++) ref[i] = -1;
  int qh = 0, qt = 0; ref[source] = 0; q[qt++] = source;
  /* on multicore the "device" arrays are the host arrays; nodes/edges are never modified by the kernels */
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
