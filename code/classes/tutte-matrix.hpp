#pragma once

#include <chrono>
#include <random>
#include <iostream>

#include "matrix.hpp"
#include "modular-integer.hpp"
#include "../config.hpp"

std::mt19937 rng(
    SEED ? SEED : (int)std::chrono::steady_clock::now().time_since_epoch().count());

/**
 * Tutte matrix class.
 * Builds a Tutte matrix over field \Integers_P.
 * @tparam P an integer such that \Integers_P is a field.
 */
template <int P>
class TutteMatrix : public Matrix<modular_int<P>> {
  /**
   * Indeterminates assigned to each edge.
   * If no indeterminate was assigned yet, it is zero.
   */
  Matrix<modular_int<P>> t; 

public:
  TutteMatrix() : Matrix<modular_int<P>>() {}

  TutteMatrix(int n) : Matrix<modular_int<P>>(n, n) {
    t = Matrix<modular_int<P>>(n, n);
  }

  TutteMatrix(Matrix<modular_int<P>> M) 
  : Matrix<modular_int<P>>(M.numRows(), M.numColumns()) {
    for (size_t i = 0; i < M.numRows(); i++) {
      for (size_t j = 0; j < M.numColumns(); j++) (*this)(i, j) = M(i, j);
    }
  }

  void addEdge(int u, int v) {
    if (t(u, v) == modular_int<MOD>()) {
      t(u, v) = rng();
      t(v, u) = -t(u, v);
    }
    (*this)(u, v) = t(u, v);
    (*this)(v, u) = t(v, u);
  }

  void removeEdge(int u, int v) {
    (*this)(u, v) = 0;
    (*this)(v, u) = 0;
  }

  // E(T) returns the edges remaining in T, i.e. entries
  // that are not zero.
  friend std::vector<std::pair<int, int>> E(const TutteMatrix& T) {
    int n = T.numRows();
    std::vector<std::pair<int, int>> edges;
    for (int u = 0; u < n; u++) {
      for (int v = u+1; v < n; v++) {
        if (T(u, v).x) edges.emplace_back(u, v);
      }
    }
    return edges;
  }

  TutteMatrix getInverse() { return TutteMatrix((*this).inverse()); }
};
