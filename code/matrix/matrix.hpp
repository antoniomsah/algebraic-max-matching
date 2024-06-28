#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <exception>
#include <tuple>
#include <vector>

#include "memory"
#include "multiplication/default-matrix-strategy.hpp"
#include "multiplication/matrix-strategy-interface.hpp"
#include "singular-exception.hpp"

using std::vector;

/**
 * @brief Matrix's class.
 **/
template <class T>
class Matrix {
  int n, m;
  vector<T> M;
  std::shared_ptr<IMatrixMultiplicationStrategy<T>> multiplicationStrategy;

 public:
  /**
   *	@brief Builds an n x m matrix, is the zero-matrix by default.
   *	Complexity: O(nm)
   **/
  Matrix(int n, int m) : n(n), m(m), M(n * m) {
    setStrategy(std::make_shared<DefaultMultiplicationStrategy<T>>());
  }

  /**
   * @brief Given a n x m std::vector builds a n x m matrix.
   * Complexity: O(nm)
   **/
  Matrix(const vector<vector<T>>& A)
      : n(A.size()), m(A.size() ? A[0].size() : 0), M(n * m) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) M[i * m + j] = A[i][j];
    }
  }

  /**
   * @brief Checks if a matrix is a square matrix.
   * A matrix is said to be a square matrix if the number of columns is equal
   * to the number of rows.
   *
   * @return true, if it is; False, otherwise
   *
   **/
  bool is_square() { return (n == m); };

  T& operator()(const int& i, const int& j) { return M[i * m + j]; }
  const T& operator()(const int& i, const int& j) const { return M[i * m + j]; }

  int num_rows() const { return n; }
  int num_columns() const { return m; }

  void setStrategy(std::shared_ptr<IMatrixMultiplicationStrategy<T>> strategy) {
    multiplicationStrategy = strategy;
  }

  /**
   * @brief Matrix addition.
   * Complexity: O(n^2)
   *
   * @return A+B
   **/
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

  /**
   * @brief Matrix subtraction.
   * Complexity: O(n^2)
   * @return A-B
   **/
  friend Matrix<T> operator-(const Matrix<T>& A, const Matrix<T>& B) {
    return A + (B * (-1));
  }

  /**
   * @brief Matrix scalar multiplication.
   * Complexity: O(n^2)
   **/
  Matrix<T> operator*(const T& t) {
    Matrix<T> C(n, m);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        C(i, j) = M(i, j) * t;
      }
    }
    return C;
  }

  /**
   * @brief Matrix scalar division.
   * Complexity: O(n^2)
   **/
  Matrix<T> operator/(const T& t) {
    Matrix<T> C(n, m);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        C(i, j) = M(i, j) / t;
      }
    }
    return C;
  }

  /**
   * @brief Matrix multiplication.
   * Complexity: O(n^3)
   **/
  Matrix<T> operator*(const Matrix<T>& A) {
    if (!multiplicationStrategy) {
      throw new std::runtime_error("Multiplication strategy not implemented.");
    }
    return multiplicationStrategy.multiply((*this), A);
  }

  Matrix<T> operator+=(const Matrix<T>& B) { return (*this) + B; }
  Matrix<T> operator-=(const Matrix<T>& B) { return (*this) - B; }
  Matrix<T> operator*=(const Matrix<T>& B) { return (*this) * B; }
  Matrix<T> operator*=(const T& t) { return (*this) * t; }
  Matrix<T> operator/=(const T& t) { return (*this) / t; }

  /**
   * @brief Computes the matrix's transpose. Complexity: O(n^2)
   *
   * @return Returns the transpose of the matrix.
   **/
  Matrix<T> transpose() {
    Matrix<T> result(m, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        result(j, i) = M(i, j);
      }
    }
    return result;
  }

  /**
   * @brief
   * This function calculates a matrix's determinant.
   *
   * @return The matrix's determinant.
   * @throws std::invalid_argument if the matrix is not a square matrix
   *
   **/
  T determinant();

  /**
   * @brief
   * This function calculates a matrix's inverse using LUP-decomposition.
   * Complexity: O(n^3).
   *
   * @return The matrix's inverse
   *
   **/
  Matrix<T> inverse();

  /**
   * @brief Checks if a matrix is singular (i.e., has inverse).
   *
   * @return True, if it has an inverse; False, otherwise.
   */
  bool is_singular() { return true; };

  /**
   * @brief Checks if a matrix is nonsingular (i.e., does not have an inverse).
   *
   * @return True, if it does not have an inverse; False, otherwise.
   */
  bool is_nonsingular() { return false; };
};