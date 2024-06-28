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

  Matrix<T> invert(const Matrix<T>& A) {
    const int n = A.num_rows(), m = A.num_columns();

    Matrix<T> NA(n, m);
    if (n == 1) {
      if (A(0, 0) == T()) {
        throw SingularMatrixError();
      }
      NA(0, 0) = T(1) / A(0, 0);
      return NA;
    }

    Matrix<T> B(n / 2, m / 2), C(n / 2, m / 2), D(n / 2, m / 2);
    for (int i = 0; i < n / 2; i++) {
      for (int j = 0; j < m / 2; j++) {
        B(i, j) = A(i, j);
      }
    }

    for (int i = n / 2; i < n; i++) {
      for (int j = 0; j < m / 2; j++) {
        C(i - n / 2, j) = A(i, j);
      }
    }

    for (int i = n / 2; i < n; i++) {
      for (int j = m / 2; j < m; j++) {
        D(i - n / 2, j - n / 2) = A(i, j);
      }
    }

    Matrix<T> NB = invert(B), CT = C.transpose(), S = D - C * NB * CT,
              NS = invert(S);

    Matrix<T> top_left = NB + NB * CT * NS * C * NB,
              top_right = NB * CT * NS * (-1), bot_left = NS * C * NB * (-1),
              bot_right = NS;

    for (int i = 0; i < n / 2; i++) {
      for (int j = 0; j < m / 2; j++) {
        NA(i, j) = top_left(i, j);
      }
    }

    for (int i = n / 2; i < n; i++) {
      for (int j = 0; j < m / 2; j++) {
        NA(i, j) = top_right(i - n / 2, j);
      }
    }

    for (int i = 0; i < n / 2; i++) {
      for (int j = m / 2; j < m; j++) {
        NA(i, j) = bot_left(i, j - m / 2);
      }
    }

    for (int i = n / 2; i < n; i++) {
      for (int j = m / 2; j < m; j++) {
        NA(i, j) = bot_right(i - n / 2, j - m / 2);
      }
    }

    return NA;
  }

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
    assert(A.n == B.n and A.m == B.m);
    Matrix<T> C(A.n, A.m);
    for (int i = 0; i < A.n; i++) {
      for (int j = 0; j < A.m; j++) {
        C(i, j) = A(i, j) - B(i, j);
      }
    }
    return C;
  }

  /**
   * @brief Matrix scalar multiplication.
   * Complexity: O(n^2)
   **/
  Matrix<T> operator*(const T& t) {
    Matrix<T> C(n, m);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        C(i, j) = (*this)(i, j) * t;
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
        C(i, j) = (*this)(i, j) / t;
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
    if (n != B.num_rows() or m != B.num_columns()) return false;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        if ((*this)(i, j) != B(i, j)) return false;
      }
    }
    return true;
  }

  bool operator!=(const Matrix<T>& B) { return not((*this) == B); }

  /**
   * @brief Computes the matrix's transpose. Complexity: O(n^2)
   *
   * @return Returns the transpose of the matrix.
   **/
  Matrix<T> transpose() {
    Matrix<T> result(m, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        result(j, i) = (*this)(i, j);
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
   * This function calculates a matrix's inverse using divider and conquer with
   * matrix multiplication. Complexity: O(n^3).
   *
   * @return The matrix's inverse
   *
   **/
  Matrix<T> inverse() {
    assert(num_rows() == num_columns());

    // First power of two NOT smaller than num_rows()
    int n = std::__bit_ceil(num_rows());
    Matrix<T> A(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (i < num_rows() and j < num_columns()) {
          A(i, j) = (*this)(i, j);
        } else {
          A(i, j) = (i == j);
        }
      }
    }

    Matrix<T> B = A.transpose() * A;

    Matrix<T> NA = invert(B) * A.transpose();

    Matrix<T> N(num_rows(), num_columns()), I(num_rows(), num_columns());
    for (int i = 0; i < num_rows(); i++) {
      for (int j = 0; j < num_columns(); j++) {
        N(i, j) = NA(i, j);
        I(i, j) = T(1);
      }
    }
    return N;
  };

  /**
   * @brief Checks if a matrix is singular (i.e., has inverse).
   *
   * @return True, if it has an inverse; False, otherwise.
   */
  bool is_singular() {
    try {
      (*this).inverse();
    } catch (const SingularMatrixError& e) {
      return true;
    }
    return false;
  };

  /**
   * @brief Checks if a matrix is nonsingular (i.e., does not have an inverse).
   *
   * @return True, if it does not have an inverse; False, otherwise.
   */
  bool is_nonsingular() { return not is_singular(); };
};