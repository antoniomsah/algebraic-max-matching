#pragma once

#include <cassert>
#include <map>
#include <vector>

#include "algorithm-strategy-interface.hpp"

using std::pair;
using std::vector;

/**
 * Implements the RankTwoAlgorithm.
 * Time complexity: O(n^4).
 */
class RankTwoAlgorithmStrategy : public IAlgorithmStrategy {

  template <const int P>
  TutteMatrix<P> rankTwoUpdate(const vector<int>& S, 
                                const TutteMatrix<P>& T,
                                const TutteMatrix<P>& N) const {
    assert(S.size() == 2);

    const int u = S[0], v = S[1];
    TutteMatrix<P> TN(2);
    TN(0, 0) = (T(u, v) * N(u, v) + 1).inv();
    TN(1, 1) = (T(u, v) * N(u, v) + 1).inv();

    return N + N('*', S) * TN * T(S, S) * N(S, '*');
  }

 public:
  /**
   * @brief Finds a perfect matching using rank two updates.
   * Complexity: O(n^4).
   */
  vector<pair<int, int>> solve(const Graph& G) const override {
    const int n = V(G).size(), m = E(G).size();

    auto T = GetTutteMatrix(G);
    if (T.isSingular()) {
      return {};
    } 

    auto N = T.getInverse();
    for (const auto& [u, v] : E(G)) {
      if (N(u, v) != T(u, v).inv() * (-1)) {
        vector<int> S = {u, v};
        N = rankTwoUpdate(S, T, N);
      } 
    }
    return E(T);
  }
};
