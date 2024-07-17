#pragma once

#include <cassert>
#include <map>
#include <vector>

#include "../../config.hpp"
#include "../../matrix/tutte-matrix.hpp"

using std::pair;
using std::vector;
class Rank2AlgorithmStrategy : public IAlgorithmStrategy {
  template <class F>
  Matrix<F> rank2_update(const vector<int>& S, const TutteMatrix<MOD>& T,
                         const Matrix<F>& N) const {
    assert(S.size() == 2);

    const int u = S[0], v = S[1];
    Matrix<F> TN(2, 2);
    TN(0, 0) = (T(u, v) * N(u, v) + 1).inv();
    TN(1, 1) = (T(u, v) * N(u, v) + 1).inv();

    return N + N('*', S) * TN * T(S, S) * N(S, '*');
  }

 public:
  /**
   * @brief Finds a perfect matching using rank2-updates.
   * Complexity: O(n^4).
   */
  vector<pair<int, int>> solve(
      const vector<vector<int>>& graph,
      const vector<pair<int, int>>& edges) const override {
    const int n = graph.size(), m = edges.size();
    TutteMatrix<MOD> T(n, edges);

    if (not T.has_perfect_matching()) {
      return {};
    }

    auto N = T.inverse();
    vector<bool> is_matching_edge(m);
    for (int i = 0; i < m; i++) {
      const auto& [u, v] = edges[i];

      if (N(u, v) != T(u, v).inv() * (-1)) {
        vector<int> S = {u, v};
        N = rank2_update(S, T, N);
      } else {
        is_matching_edge[i] = true;
      }
    }

    vector<pair<int, int>> matching;
    for (int i = 0; i < m; i++) {
      if (is_matching_edge[i]) matching.push_back(edges[i]);
    }
    return matching;
  }
};
