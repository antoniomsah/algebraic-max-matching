#pragma once

#include <vector>

#include "../../config.hpp"
#include "../algorithm-strategy-interface.hpp"
#include "../../matrix/tutte-matrix.hpp"

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
    const size_t n = graph.size(), m = edges.size();

    TutteMatrix T = build_tutte_matrix<MOD>(n, edges);

    if (T.is_singular()) {
      return {};
    }

    vector<bool> is_matching_edge(m);
    for (size_t i = 0; i < m; i++) {
      const auto& [u, v] = edges[i];
      const modular_int<MOD> val = T(u, v);

      // remove edge uv 
      T(u, v) = 0;
      T(v, u) = 0;

      if (T.is_singular()) {
        is_matching_edge[i] = true;

        // add edge uv
        T(u, v) = val;
        T(v, u) = val * (-1);
      }
    }

    vector<pair<int, int>> matching;
    for (size_t i = 0; i < m; i++) {
      if (is_matching_edge[i]) matching.push_back(edges[i]);
    }
    return matching;
  }
};