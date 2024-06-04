#pragma once

#include <chrono>
#include <random>

#include "matrix.hpp"

std::mt19937 rng(
    (int)std::chrono::steady_clock::now().time_since_epoch().count());

/**
 * @brief Implements a Tutte's matrix in the field Z_{P}.
 *
 * @tparam an integer P such that Z_{P} is a field.
 */
template <const int P = 998244353>
class TutteMatrix : Matrix<int> {
 public:
  /**
   * @brief Construct a new Tutte Matrix. Complexity: O(n^2).
   *
   * @param edges std::vector<pair<int,int>> such that each entry is a
   * pair {u,v} that represents and undirected edge between u and v.
   */
  TutteMatrix(const int& num_vertices,
              const std::vector<std::pair<int, int>>& edges)
      : Matrix<int>(num_vertices, num_vertices) {
    for (long unsigned int e = 0; e < edges.size(); e++) {
      auto [i, j] = edges[e];
      assert(i != j);

      int val = rng() % P;
      (*this)[i][j] = val;
      (*this)[j][i] = P - val;
    }
  }
};