#pragma once

#include <cmath>
#include <iostream>

#include "./matrix-strategy-interface.hpp"

// StrassenMultiplicationStrategy implements matrix multiplication in O(n^{lg7}).
template <typename T>
class StrassenMultiplicationStrategy : public IMatrixMultiplicationStrategy<T> {
  // extendMatrix extends a matrix to size x size matrix.
  // The new entries are set as 0.
  Matrix<T> extendMatrix(const Matrix<T>& M, const int size) const {
    // First power of two NOT smaller than numRows().
    Matrix<T> ret(size, size);
    for (size_t i = 0; i < M.numRows(); i++) {
      for (size_t j = 0; j < M.numColumns(); j++) {
        ret(i, j) = M(i, j);
      }
    }
    return ret;
  }

  Matrix<T> strassenMultiply(const Matrix<T> A, const Matrix<T> &B) const {
    const size_t n = A.numRows();
    if (n == 1) {
      return Matrix<T>(1, 1, A(0, 0) * B(0, 0));
    }

    const size_t new_size = n/2;

    // Get matrix quadrants.
    auto [A11, A12,
          A21, A22] = split(A);

    auto [B11, B12,
          B21, B22] = split(B);

    // Compute the 7 products
    auto M1 = strassenMultiply(A11+A22, B11+B22);
    auto M2 = strassenMultiply(A21+A22, B11);
    auto M3 = strassenMultiply(A11, B12-B22);
    auto M4 = strassenMultiply(A22, B21-B11);
    auto M5 = strassenMultiply(A11+A12, B22);
    auto M6 = strassenMultiply(A21-A11, B11+B12);
    auto M7 = strassenMultiply(A12-A22, B21+B22);

    // Compute the quadrants of the result
    auto C11 = M1 + M4 - M5 + M7;
    auto C12 = M3 + M5;
    auto C21 = M2 + M4;
    auto C22 = M1 + M3 - M2 + M6;

    return unite(C11, C12, 
                 C21, C22);
  }

 public:
  Matrix<T> multiply(const Matrix<T>& A, const Matrix<T>& B) const override {
    assert(A.numColumns() == B.numRows());

    // Next power of two not smaller than any of the dimensions.
    const size_t maxDim = std::max({A.numRows(), A.numColumns(), B.numRows(), B.numColumns()});
    const size_t new_size = std::__bit_ceil(maxDim);

    // Pad the matrices.
    auto newA = extendMatrix(A, new_size), newB = extendMatrix(B, new_size);

    // Multiply the padded matrices.
    auto C = strassenMultiply(newA, newB);

    // Extract the result from C.
    Matrix<T> result(A.numRows(), B.numColumns());
    for (size_t i = 0; i < A.numRows(); i++) {
      for (size_t j = 0; j < B.numColumns(); j++) {
        result(i, j) = C(i, j);
      }
    }
    return result;
  }
};

