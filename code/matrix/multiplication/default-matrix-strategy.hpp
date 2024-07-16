#pragma once

#include "./matrix-strategy-interface.hpp"

template <typename T>
class DefaultMultiplicationStrategy : public IMatrixMultiplicationStrategy<T> {
 public:
  Matrix<T> multiply(const Matrix<T>& A, const Matrix<T>& B) const override {
    assert(A.num_columns() == B.num_rows());
    Matrix<T> C(A.num_rows(), B.num_columns());
    for (int i = 0; i < A.num_rows(); i++) {
      for (int j = 0; j < B.num_columns(); j++) {
        for (int k = 0; k < A.num_columns(); k++) {
          C(i, j) += A(i, k) * B(k, j);
        }
      }
    }
    return C;
  }
};