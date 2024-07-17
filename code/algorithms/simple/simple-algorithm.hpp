#pragma once

#include <vector>

#include "../../config.hpp"
#include "../algorithm-strategy-interface.hpp"

using namespace std;

class SimpleAlgorithmStrategy : public IAlgorithmStrategy {
 public:
  /**
   * @brief Finds a perfect matching.
   * Complexity: O(n^{\omega+2}).
   */
  vector<pair<int, int>> solve(
      const vector<vector<int>>& graph,
      const vector<pair<int, int>>& edges) const override {
    const int n = graph.size(), m = edges.size();
    TutteMatrix<MOD> T(n, edges);

    if (not T.has_perfect_matching()) {
      return {};
    }

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
};