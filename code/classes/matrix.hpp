#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <exception>
#include <iomanip>
#include <tuple>
#include <vector>

#include "memory"
#include "multiplication/default-matrix-strategy.hpp"
#include "multiplication/matrix-strategy-interface.hpp"
#include "singular-exception.hpp"

using std::vector;

/**
 * Matrix's class.
 * @tparam T Type of the matrices entries (MUST have addition and multiplication).
 **/
template <class T>
class Matrix {
protected:
  vector<T> M;

private:
  int n, m;
  std::shared_ptr<IMatrixMultiplicationStrategy<T>> multiplicationStrategy;

  /**
   * invert computes a matrix inverse using matrix multiplication.
   * Caution: matrix A must be positive semi-definite.
   * Complexity: O(O(multiply)).
   *
   * @return the inverse of matrix A.
   */
  Matrix<T> invert(const Matrix<T>& A);

 public:
  Matrix() : n(0), m(0), M(0) {}

  /**
   *	@brief Builds an n x m matrix, is the zero-matrix by default.
   *	Complexity: O(nm)
   **/
  Matrix(int n, int m) : n(n), m(m), M(n * m) {
    setStrategy(std::make_shared<DefaultMultiplicationStrategy<T>>());
  }

  /**
   * Given a n x m std::vector builds a n x m matrix.
   * Complexity: O(nm)
   **/
  Matrix(const vector<vector<T>>& A)
      : n(A.size()), m(A.size() ? A[0].size() : 0), M(n * m) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) M[i * m + j] = A[i][j];
    }
  }

  // isSquare returns true if the matrix is a square matrix.
  bool isSquare() const { return (n == m); };

  T& operator()(const int& i, const int& j) { return M[i * m + j]; }
  const T& operator()(const int& i, const int& j) const { return M[i * m + j]; }

  size_t numRows() const { return n; }
  size_t numColumns() const { return m; }

  void setStrategy(
      std::shared_ptr<IMatrixMultiplicationStrategy<T>> strategy) {
    multiplicationStrategy = strategy;
  }

  friend Matrix<T> operator+(const Matrix<T>& A, const Matrix<T>& B) {
    assert(A.n == B.n and A.m == B.m);
    Matrix<T> C(A.n, A.m);
    for (int i = 0; i < A.n; i++) {
      for (int j = 0; j < A.m; j++) {
        C(i, j) = A(i, j) + B(i, j);
      }
    }
    return C;
  }

  friend Matrix<T> operator-(const Matrix<T>& A, const Matrix<T>& B) {
    assert(A.n == B.n and A.m == B.m);
    Matrix<T> C(A.n, A.m);
    for (int i = 0; i < A.n; i++) {
      for (int j = 0; j < A.m; j++) {
        C(i, j) = A(i, j) - B(i, j);
      }
    }
    return C;
  }

  Matrix<T> operator-() { return (*this) * (-1); }

  Matrix<T> operator*(const T& t) {
    Matrix<T> C(n, m);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        C(i, j) = (*this)(i, j) * t;
      }
    }
    return C;
  }

  Matrix<T> operator/(const T& t) {
    Matrix<T> C(n, m);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        C(i, j) = (*this)(i, j) / t;
      }
    }
    return C;
  }

  Matrix<T> operator*(const Matrix<T>& A) {
    if (!multiplicationStrategy) {
      throw std::runtime_error("Multiplication strategy not implemented.");
    }
    return multiplicationStrategy->multiply((*this), A);
  }

  Matrix<T> operator+=(const Matrix<T>& B) { return (*this) + B; }
  Matrix<T> operator-=(const Matrix<T>& B) { return (*this) - B; }
  Matrix<T> operator*=(const Matrix<T>& B) { return (*this) * B; }
  Matrix<T> operator*=(const T& t) { return (*this) * t; }
  Matrix<T> operator/=(const T& t) { return (*this) / t; }

  bool operator==(const Matrix<T>& B) {
    if ((*this).num_rows() != B.num_rows() or
        (*this).num_columns() != B.num_columns())
      return false;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        if ((*this)(i, j) != B(i, j)) return false;
      }
    }
    return true;
  }

  bool operator!=(const Matrix<T>& B) { return not((*this) == B); }

  /**
   * Computes matrix M[R, C].
   * Complexity: O(|R|*|C|).
   */
  Matrix<T> operator()(const vector<int>& R, const vector<int>& C) const {
    Matrix<T> M(R.size(), C.size());
    for (size_t i = 0; i < R.size(); i++) {
      for (size_t j = 0; j < C.size(); j++) {
        M(i, j) = (*this)(R[i], C[j]);
      }
    }
    return M;
  }

  /**
   * Computes matrix M[*, S].
   * Important: M[*, S] column entries will be in S order.
   * Complexity: O(n*|S|).
   */
  Matrix<T> operator()(char c, const vector<int>& S) const {
    assert(c == '*');
    Matrix<T> M(numRows(), S.size());
    for (size_t i = 0; i < numRows(); i++) {
      for (size_t j = 0; j < S.size(); j++) {
        M(i, j) = (*this)(i, S[j]);
      }
    }
    return M;
  }

  /**
   * Computes matrix M[S, *].
   * Important: M[S, *] row entries will be in S order.
   * Complexity: O(|S|*m).
   */
  Matrix<T> operator()(const vector<int>& S, char c) const {
    assert(c == '*');
    Matrix<T> M(S.size(), numColumns());
    for (size_t i = 0; i < S.size(); i++) {
      for (size_t j = 0; j < numColumns(); j++) {
        M(i, j) = (*this)(S[i], j);
      }
    }
    return M;
  }

  // transpose computes the matrix's transpose. Complexity: O(n^2)
  Matrix<T> transpose() {
    Matrix<T> result(m, n);
    for (size_t i = 0; i < numRows(); i++) {
      for (size_t j = 0; j < numColumns(); j++) {
        result(j, i) = (*this)(i, j);
      }
    }
    return result;
  }

  /**
   * inverse calculates a matrix's inverse.
   * Time complexity: O(O(multiply)), where O(multiply) is the time complexity for the matrix multiplication algorithm used.
   **/
  Matrix<T> inverse();

  // rank returns the rank of the matrix.
  int rank();

  /**
   * isSingular returns true if the matrix is singular.
   * Reminder: a singular matrix is a matrix that does not have an inverse.
   */
  bool isSingular() {
    try {
      (*this).inverse();
    } catch (const SingularMatrixError& e) {
      return true;
    }
    return false;
  }

  // isNonSingular returns true if the matrix is not singular.
  bool isNonSingular() { return not isSingular(); }
};

// Implementations of inverse, rank and invert.

template <typename T>
Matrix<T> Matrix<T>::inverse() {
  assert(isSquare());

  // First power of two NOT smaller than numRows().
  const size_t n = std::__bit_ceil(numRows());
  Matrix<T> A(n, n);
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      if (i < numRows() and j < numColumns()) {
        A(i, j) = (*this)(i, j);
      } else {
        A(i, j) = (i == j);
      }
    }
  }

  Matrix<T> B = A.transpose() * A;

  Matrix<T> NA = invert(B) * A.transpose();

  Matrix<T> N(numRows(), numColumns());
  for (size_t i = 0; i < numRows(); i++) {
    for (size_t j = 0; j < numColumns(); j++) {
      N(i, j) = NA(i, j);
    }
  }
  return N;
}

template <typename T>
Matrix<T> Matrix<T>::invert(const Matrix<T>& A) {
  const int n = A.numRows(), m = A.numColumns();

  Matrix<T> NA(n, m);
  if (n == 1) {
    if (A(0, 0) == T()) {
      throw SingularMatrixError();
    }
    NA(0, 0) = T(1) / A(0, 0);
    return NA;
  } 

  auto [B, CT, C, D] = split(A);
  Matrix<T> NB = invert(B), 
            S = D - C * NB * CT, // Schur complement.
            NS = invert(S);

  return unite(
    NB + NB * CT * NS * C * NB, -NB * CT * NS,
    -NS * C * NB              , NS
  );
}

template <typename T>
int Matrix<T>::rank() {
  Matrix<T> A = (*this);
  const int n = A.numRows(), m = A.numColumns();
  int rnk = 0;
  vector<bool> selected(n, false);
  for (int i = 0; i < m; i++) {
    int j;
    for (j = 0; j < n; j++) {
      if (!selected[j] and A(i, j) != T()) break;
    }

    if (j == n) {
      continue;
    }
    rnk++;
    selected[j] = true;
    for (int p = i+1; p < m; p++) {
      A(j, p) /= A(j, i);
    }
    for (int k = 0; k < n; k++) {
      if (k == j or A(k, i) == T()) continue;
      for (int p = i+1; p < m; p++) {
        A(k, p) -= A(j, p) * A(k, i);
      }
    }
  }
  return rnk;
}