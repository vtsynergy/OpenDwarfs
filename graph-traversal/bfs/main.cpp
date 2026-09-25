#include "dwarf_common.hpp"

#include <fstream>
#include <vector>

namespace {

bool LoadGraph(const std::string &path, std::vector<int> &row_ptr,
               std::vector<int> &col_idx, int &nodes, int &edges, int &source) {
  std::ifstream in(path);
  if (!in) {
    return false;
  }
  if (!(in >> nodes) || nodes < 1) return false;
  std::vector<int> starts(nodes), counts(nodes);
  for (int i = 0; i < nodes; ++i) {
    if (!(in >> starts[i] >> counts[i]) || starts[i] < 0 || counts[i] < 0)
      return false;
  }
  if (!(in >> source >> edges) || source < 0 || source >= nodes || edges < 0)
    return false;
  long long expected = 0;
  for (int i = 0; i < nodes; ++i) {
    if (starts[i] != expected) return false;
    expected += counts[i];
  }
  if (expected != edges) return false;
  col_idx.resize(edges);
  for (int e = 0; e < edges; ++e) {
    int dst = 0;
    int cost = 0;
    if (!(in >> dst >> cost) || dst < 0 || dst >= nodes) return false;
    (void)cost;
    col_idx[e] = dst;
  }
  row_ptr.assign(nodes + 1, 0);
  for (int i = 0; i < nodes; ++i) {
    row_ptr[i] = starts[i];
  }
  row_ptr[nodes] = edges;
  return true;
}

void InitGraph(std::vector<int> &row_ptr, std::vector<int> &col_idx, int nodes,
               int degree) {
  row_ptr.resize(nodes + 1);
  col_idx.resize(static_cast<std::size_t>(nodes) * degree);
  for (int i = 0; i <= nodes; ++i) {
    row_ptr[i] = i * degree;
  }
  for (int i = 0; i < nodes; ++i) {
    for (int k = 0; k < degree; ++k) {
      const int edge = i * degree + k;
      const int stride = 1 + k * 13;
      col_idx[edge] = (i + stride) % nodes;
    }
  }
}

void Run(std::vector<int> &row_ptr, std::vector<int> &col_idx,
         std::vector<unsigned char> &frontier,
         std::vector<unsigned char> &next_frontier,
         std::vector<int> &visited, std::vector<int> &cost, int nodes,
         int edges, int max_levels) {
  int *pr = row_ptr.data();
  int *pc = col_idx.data();
  unsigned char *pf = frontier.data();
  unsigned char *pn = next_frontier.data();
  int *pv = visited.data();
  int *pco = cost.data();
#if DWARF_USE_TARGET
#pragma omp target data map(to : pr[0 : nodes + 1], pc[0 : edges]) map(tofrom : pf[0 : nodes], pn[0 : nodes], pv[0 : nodes], pco[0 : nodes])
  {
#endif
    for (int level = 0; level < max_levels; ++level) {
      int changed = 0;
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for reduction(+ : changed) thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for reduction(+ : changed) schedule(dynamic, 256)
#endif
      for (int i = 0; i < nodes; ++i) {
        if (!pf[i]) {
          continue;
        }
        pf[i] = 0;
        for (int e = pr[i]; e < pr[i + 1]; ++e) {
          const int dst = pc[e];
          int already_seen;
#pragma omp atomic capture
          { already_seen = pv[dst]; pv[dst] = 1; }
          if (!already_seen) {
            pco[dst] = level + 1;
            pn[dst] = 1;
            changed += 1;
          }
        }
      }
#if DWARF_USE_TARGET
#pragma omp target teams distribute parallel for thread_limit(DWARF_THREAD_LIMIT)
#else
#pragma omp parallel for schedule(static)
#endif
      for (int i = 0; i < nodes; ++i) {
        pf[i] = pn[i];
        pn[i] = 0;
      }
      if (changed == 0) {
        break;
      }
    }
#if DWARF_USE_TARGET
  }
#endif
}

bool Verify(const std::vector<int> &row_ptr, const std::vector<int> &col_idx,
            const std::vector<int> &cost, int source, int max_levels) {
  std::vector<int> expected(cost.size(), -1), queue;
  expected[source] = 0;
  queue.push_back(source);
  for (std::size_t head = 0; head < queue.size(); ++head) {
    const int node = queue[head];
    if (expected[node] >= max_levels) continue;
    for (int e = row_ptr[node]; e < row_ptr[node + 1]; ++e) {
      const int dst = col_idx[e];
      if (expected[dst] == -1) {
        expected[dst] = expected[node] + 1;
        queue.push_back(dst);
      }
    }
  }
  return cost == expected;
}

}  // namespace

int main(int argc, char **argv) {
  const auto opts = dwarfs_omp::ParseOptions(argc, argv);
  int nodes = static_cast<int>(opts.size);
  const int degree = static_cast<int>(opts.secondary_size);
  int edges = nodes * degree;
  int source = 0;
  std::vector<int> row_ptr, col_idx, cost(nodes);
  std::vector<unsigned char> frontier(nodes), next_frontier(nodes);
  std::vector<int> visited(nodes);
  bool valid = true;

  if (!opts.input_file.empty()) {
    if (!LoadGraph(opts.input_file, row_ptr, col_idx, nodes, edges, source)) {
      return 3;
    }
  } else {
    InitGraph(row_ptr, col_idx, nodes, degree);
  }
  cost.resize(nodes);
  frontier.resize(nodes);
  next_frontier.resize(nodes);
  visited.resize(nodes);
  auto times = dwarfs_omp::AdaptiveRun(opts,
      [&]() {
        std::fill(frontier.begin(), frontier.end(), 0);
        std::fill(next_frontier.begin(), next_frontier.end(), 0);
        std::fill(visited.begin(), visited.end(), 0);
        std::fill(cost.begin(), cost.end(), -1);
        frontier[source] = 1; visited[source] = 1; cost[source] = 0;
      },
      [&]() { Run(row_ptr, col_idx, frontier, next_frontier, visited, cost, nodes, edges, opts.iterations); });
  if (opts.verify) {
    valid = Verify(row_ptr, col_idx, cost, source, opts.iterations);
  }
  dwarfs_omp::PrintResult("bfs", "graph-traversal",
                          static_cast<std::size_t>(nodes),
                          opts.input_file.empty() ? opts.secondary_size
                                                  : static_cast<std::size_t>(edges),
                          opts.iterations, static_cast<int>(times.size()),
                          times, valid);
  return valid ? 0 : 2;
}
