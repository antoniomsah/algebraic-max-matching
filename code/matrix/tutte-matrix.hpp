#pragma once

#include <chrono>
#include <random>

#include "matrix.hpp"
#include "modular-integer.hpp"

std::mt19937 rng(
    SEED ? SEED
         : (int)std::chrono::steady_clock::now().time_since_epoch().count());

template <const int P>
using TutteMatrix = Matrix<modular_int<P>>;

template <const int T>
TutteMatrix<T> build_tutte_matrix(const size_t& n,
                                  const std::vector<pair<int, int>>& edges) {
  TutteMatrix<T> M(n, n);
  for (size_t i = 0; i < edges.size(); i++) {
    modular_int<T> val = rng();
    const auto& [u, v] = edges[i];
    M(u, v) = val;
    M(v, u) = val * (-1);
  }
  return M;
}
