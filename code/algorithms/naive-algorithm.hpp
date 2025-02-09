#pragma once

#include <vector>

#include "algorithm-strategy-interface.hpp"

using namespace std;

/**
 * Implements the NaiveAlgorithm. 
 * Time complexity: O(n^2 O(multiply)), where O(multiply) is the time complexity of the matrix multiplication algorithm used.
 */
class NaiveAlgorithmStrategy : public IAlgorithmStrategy {
 public:
  /**
   * Finds a perfect matching.
   * Complexity: O(n^{\omega+2}).
   */
  vector<pair<int, int>> solve(const Graph& G) const override {
    auto T = GetTutteMatrix(G);
    if (T.isSingular()) {
      return {};
    }

    vector<pair<int, int>> M;
    for (const auto& [u, v] : E(G)) {
      T.removeEdge(u, v);
      if (T.isSingular()) { // This edge is essential.
        T.addEdge(u, v);
      }
    }
    return E(T);
  }

  string name() const override {
    return "Simple";
  }
};
