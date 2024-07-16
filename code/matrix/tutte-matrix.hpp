#pragma once

#include <chrono>
#include <random>

#include "matrix.hpp"
#include "modular-integer.hpp"

std::mt19937 rng(
    (int)std::chrono::steady_clock::now().time_since_epoch().count());

/**
 * @brief Implements a Tutte's matrix in the field Z_{P}.
 *
 * @tparam P an integer such that Z_{P} is a field.
 */
template <const int P>
class TutteMatrix {
  using F = modular_int<P>;

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
        val[i] = rng();
      }
      add_edge(i);
    }
  }

  Matrix<F> get_matrix() { return T; }

  F& operator()(const int& i, const int& j) { return T(i, j); }
  F operator()(const int& i, const int& j) const { return T(i, j); }

  void add_edge(int i) {
    auto [u, v] = edges[i];
    T(u, v) = val[i];
    T(v, u) = val[i] * -1;
  }

  void remove_edge(int i) {
    auto [u, v] = edges[i];
    T(u, v) = F();
    T(v, u) = F();
  }

  /**
   * @brief Removes a vertex.
   * CAUTION: At the moment cannot add a vertex back!
   *
   * @param u
   */
  void remove_vertice(int u) {
    for (int i = 0; i < T.num_rows(); i++) {
      T(u, i) = T(i, u) = F();
    }
    T(u, u) = 1;
  }

  // Matrices operations
  Matrix<F> inverse() { return T.inverse(); }

  Matrix<F> operator()(const vector<int>& R, const vector<int>& C) const {
    return T(R, C);
  }

  bool has_perfect_matching() { return T.is_nonsingular(); }
};