#pragma once

#include <cassert>
#include <map>
#include <vector>

#include "algorithm-strategy-interface.hpp"
#include "../utils.hpp"

using std::pair;
using std::vector;

/**
 * Implements the RankTwoAlgorithm.
 * Time complexity: O(n^4).
 */
class RankTwoAlgorithmStrategy : public IAlgorithmStrategy {
 public:
  /**
   * @brief Finds a perfect matching using rank two updates.
   * Complexity: O(n^4).
   */
  vector<pair<int, int>> solve(const Graph& G) const override {
    auto T = GetTutteMatrix(G);
    if (T.isSingular()) {
      return {};
    } 

    auto N = T.inverse();
    for (const auto& [u, v] : E(G)) {
      if (N(u, v) != -T(u, v).inv()) { // This edge can be removed.
        vector<int> S = {u, v};
        N = rankTwoUpdate(S, T, N);
        T.removeEdge(u, v);
      } 
    }
    return E(T);
  }
};
