#pragma once

#include "./matrix-strategy-interface.hpp"

template <typename T>
class DefaultMultiplicationStrategy : public IMatrixMultiplicationStrategy<T> {
 public:
  Matrix<T> multiply(const Matrix<T>& A, const Matrix<T>& B) const override {
    assert(A.numColumns() == B.numRows());
    Matrix<T> C(A.numRows(), B.numColumns());
    for (size_t i = 0; i < A.numRows(); i++) {
      for (size_t j = 0; j < B.numColumns(); j++) {
        for (size_t k = 0; k < A.numColumns(); k++) {
          C(i, j) += A(i, k) * B(k, j);
        }
      }
    }
    return C;
  }
};