#pragma once

#include <chrono>
#include <random>

#include "matrix.hpp"

std::mt19937 rng(
    (int)std::chrono::steady_clock::now().time_since_epoch().count());

/**
 * @brief Implements a Tutte's matrix in the field Z_{P}.
 *
 * @tparam P an integer such that Z_{P} is a field.
 */
template <class F, const int P>
class TutteMatrix {
 private:
  Matrix<F> T;
  std::vector<std::pair<int, int>> edges;
  std::vector<F> val;

 public:
  /**
   * @brief Construct a new Tutte Matrix. Complexity: O(n^2).
   *
   * @param edges std::vector<pair<int,int>> such that each entry is a
   * pair {u,v} that represents and undirected edge between u and v.
   */
  TutteMatrix(const int& num_vertices,
              const std::vector<std::pair<int, int>>& edges)
      : T(num_vertices, num_vertices), edges(edges), val(edges.size()) {
    for (unsigned int i = 0; i < edges.size(); i++) {
      while (val[i] == 0) {
        val[i] = rng() % P;
      }
      add_edge(i);
    }
  }

  void add_edge(int i) {
    auto [u, v] = edges[i];
    T[u][v] = val[i];
    T[v][u] = val[i] * -1;
  }

  void remove_edge(int i) {
    auto [u, v] = edges[i];
    T[u][v] = 0;
    T[v][u] = 0;
  }

  bool has_perfect_matching() { return T.is_nonsingular(); }
};