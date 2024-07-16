#pragma once

#include <map>
#include <vector>

#include "../../config.hpp"
#include "../../matrix/tutte-matrix.hpp"

namespace Rank2Algorithm {

using std::pair;
using std::vector;

template <class F>
Matrix<F> rank2_update(const vector<int>& S, const TutteMatrix<MOD>& T,
                       const Matrix<F>& N) {
  assert(S.size() == 2);

  Matrix<F> I(2, 2);
  I(0, 0) = I(1, 1) = 1;
  Matrix<F> TN = (I - T(S, S) * N(S, S)).inverse();

  return N + N('*', S) * TN * T(S, S) * N(S, '*');
}

/**
 * @brief Finds a perfect matching using rank2-updates.
 * Complexity: O(n^4).
 */
vector<pair<int, int>> solve(const vector<vector<int>>& graph,
                             const vector<pair<int, int>>& edges) {
  const int n = graph.size(), m = edges.size();
  TutteMatrix<MOD> T(n, edges);
  auto N = T.inverse();

  vector<bool> is_matching_edge(m);
  for (int i = 0; i < m; i++) {
    const auto& [u, v] = edges[i];

    try {
      vector<int> S = {u, v};
      N = rank2_update(S, T, N);
      T.remove_edge(i);
    } catch (std::exception e) {
      is_matching_edge[i] = true;
    }
  }

  vector<pair<int, int>> matching;
  for (int i = 0; i < m; i++) {
    if (is_matching_edge[i]) matching.push_back(edges[i]);
  }
  return matching;
}

}  // namespace Rank2Algorithm