#pragma once

#include <chrono>
#include <random>
#include <iostream>

#include "matrix.hpp"
#include "modular-integer.hpp"
#include "../config.hpp"

std::mt19937 rng(
    SEED ? SEED : (int)std::chrono::steady_clock::now().time_since_epoch().count());

template <int P>
class TutteMatrix : public Matrix<modular_int<P>> {
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
      t(v, u) = t(u, v) * (-1);
    }
    (*this)(u, v) = t(u, v);
    (*this)(v, u) = t(v, u);
  }

  void removeEdge(int u, int v) {
    (*this)(u, v) = 0;
    (*this)(v, u) = 0;
  }

  TutteMatrix getInverse() { return TutteMatrix((*this).inverse()); }
};
