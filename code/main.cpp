#include <iostream>
#include <map>
#include <vector>

#include "algorithms/solver.hpp"

using namespace std;

int main(int argc, char* argv[]) {
  MatchingSolver::initialize();

  int n, m;
  cin >> n >> m;

  vector<vector<int>> graph(n, vector<int>(n));
  vector<pair<int, int>> edges(m);
  for (int i = 0; i < m; i++) {
    int u, v;
    cin >> u >> v;
    edges[i] = {u, v};
    graph[u][v] = graph[v][u] = 1;
  }

  assert(argc >= 2);

  MatchingSolver solver(graph, edges);
  solver.set_strategy(stoi(argv[1]));
  solver.PerfectMatching();

  return 0;
}
