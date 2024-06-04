#pragma once

#include <map>
#include <vector>

#include "config.hpp"
#include "tutte-matrix.hpp"

/**
 * @brief Implements the Rabin and Vazirani's improvement
 * to the "simple" algorithm.
 */
namespace RabinVaziraniAlgorithm {

using namespace std;

vector<pair<int, int>> solve(const vector<vector<int>> graph,
                             const vector<pair<int, int>>& edges) {
  const int n = graph.size();
  TutteMatrix<MOD> T(n, edges);

  bool found;
  vector<pair<int, int>> matching;
  do {
    found = false;
    auto N = T.inverse();
    for (int u = 0; not found and u < n; u++) {
      for (int v = 0; not found and v < n; v++) {
        if (u == v or N[u][v] == 0) continue;
        found = true;
        matching.emplace_back(u, v);
      }
    }

    if (found) {
      auto [u, v] = matching.back();
      T.remove_vertice(u);
      T.remove_vertice(v);
    }
  } while (found);
  return matching;
}

}  // namespace RabinVaziraniAlgorithm