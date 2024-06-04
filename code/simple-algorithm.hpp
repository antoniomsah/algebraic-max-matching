#pragma once

#include <map>
#include <vector>

#include "config.hpp"
#include "tutte-matrix.hpp"

namespace SimpleAlgorithm {

using namespace std;

/**
 * @brief Finds a perfect matching. Complexity: O(n^{\omega+2}).
 *
 * @param n number of vertices
 * @param edges std::vector of edges
 */
vector<pair<int, int>> solve(const vector<vector<int>>& graph,
                             const vector<pair<int, int>>& edges) {
  const int n = graph.size(), m = edges.size();
  TutteMatrix<MOD> T(n, edges);
  vector<bool> is_matching_edge(m);
  for (int i = 0; i < m; i++) {
    T.remove_edge(i);
    if (not T.has_perfect_matching()) {
      is_matching_edge[i] = true;
      T.add_edge(i);
    }
  }

  vector<pair<int, int>> matching;
  for (int i = 0; i < m; i++) {
    if (is_matching_edge[i]) matching.push_back(edges[i]);
  }
  return matching;
}

}  // namespace SimpleAlgorithm