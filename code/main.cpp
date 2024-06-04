#include <iostream>

#include "simple-algorithm.hpp"

using namespace std;

int main() {
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

  SimpleAlgorithm::solve(graph, edges, true);

  return 0;
}
