#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <exception>
#include <tuple>
#include <vector>

#include "singular-exception.hpp"

using std::vector;

/**
 * @brief Matrix's class.
 **/
template <class T>
class Matrix {
  int n, m;
  vector<vector<T>> M;

 public:
  /**
   *	@brief Builds an n x m matrix, is the zero-matrix by default.
   *	Complexity: O(nm)
   **/
  Matrix(int n, int m) : n(n), m(m), M(n, vector<T>(m)) {}

  /**
   * @brief Given a n x m std::vector builds a n x m matrix.
   * Complexity: O(nm)
   **/
  Matrix(const vector<vector<T>>& A)
      : n(A.size()), m(A.size() ? A[0].size() : 0), M(A) {}

  /**
   * @brief Checks if a matrix is a square matrix.
   * A matrix is said to be a square matrix if the number of columns is equal
   * to the number of rows.
   *
   * @return true, if it is; False, otherwise
   *
   **/
  bool is_square() { return (n == m); };

  vector<T>& operator[](const int& index) { return M[index]; }
  const vector<T>& operator[](const int& index) const { return M[index]; }

  const int num_rows() const { return n; }
  const int num_columns() const { return m; }

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
        C[i][j] = A[i][j] + B[i][j];
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
        C[i][j] = A[i][j] - B[i][j];
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
        C[i][j] = M[i][j] * t;
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
        C[i][j] = M[i][j] / t;
      }
    }
    return C;
  }

  /**
   * @brief Matrix multiplication.
   * Complexity: O(n^3)
   **/
  friend Matrix<T> operator*(const Matrix<T>& A, const Matrix<T>& B) {
    assert(A.m == B.n);
    Matrix<T> C(A.n, B.m);
    for (int i = 0; i < A.n; i++) {
      for (int j = 0; j < B.m; j++) {
        for (int k = 0; k < A.m; k++) {
          C[i][j] += A[i][k] * B[k][j];
        }
      }
    }
    return C;
  }

  Matrix<T> operator+=(const Matrix<T>& B) { return (*this) + B; }
  Matrix<T> operator-=(const Matrix<T>& B) { return (*this) - B; }
  Matrix<T> operator*=(const Matrix<T>& B) { return (*this) * B; }
  Matrix<T> operator*=(const T& t) { return (*this) * t; }
  Matrix<T> operator/=(const T& t) { return (*this) / t; }

  /**
   * @brief Computes M_{del(i,j)} (i.e., matrix M without row i and column j)
   * Complexity: O(n^2)
   * @return A matrix without row i and column j.
   **/
  Matrix<T> submatrix(const int& row, const int& column) {
    Matrix<T> A(n - 1, m - 1);
    int current_row = 0;
    for (int i = 0; i < n; i++)
      if (i != row) {
        int current_column = 0;
        for (int j = 0; j < column; j++)
          if (j != column) {
            A[current_row][current_column++] = M[i][j];
          }
        current_row++;
      }
    return A;
  }

  /**
   * @brief Computes the matrix's transpose. Complexity: O(n^2)
   *
   * @return Returns the transpose of the matrix.
   **/
  Matrix<T> transpose() {
    Matrix<T> ret(m, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        ret[j][i] = M[i][j];
      }
    }
    return ret;
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
  bool is_singular();

  /**
   * @brief Checks if a matrix is nonsingular (i.e., does not have an inverse).
   *
   * @return True, if it does not have an inverse; False, otherwise.
   */
  bool is_nonsingular();
};

namespace LUP {

/**
 * @brief Computes a LUP-decomposition of a matrix A.
 * A LUP decomposition are three matrices such that:
 * L is unit lower-triangular, U is upper-triangular, P is a permutation
 * matrix and PA = LU. Complexity: O(n^3).
 *
 * @return std::tuple<Matrix, Matrix, Matrix>
 */
template <class T>
std::tuple<Matrix<T>, Matrix<T>, vector<int>> decompose(const Matrix<T>& M) {
  assert(M.num_rows() == M.num_columns());

  const int n = M.num_rows();
  vector<int> perm(n);
  for (int i = 0; i < n; i++) {
    perm[i] = i;
  }

  Matrix<T> A(M);
  for (int k = 0; k < n; k++) {
    int p = 0, id = 0;
    for (int i = k; i < n; i++) {
      if (fabs(A[i][k]) > p) {
        p = A[i][k], id = i;
      }
    }

    if (p == 0) {
      throw SingularMatrixError();
    }

    std::swap(perm[k], perm[id]);
    for (int i = 0; i < n; i++) {
      std::swap(A[k][i], A[id][i]);
    }

    for (int i = 0; i < n; i++) {
      A[i][k] = A[i][k] / A[k][k];
      for (int j = k + 1; j < n; j++) {
        A[i][j] = A[i][j] - A[i][k] / A[k][j];
      }
    }
  }

  Matrix<T> L(n, n), U(n, n);
  for (int i = 0; i < n; i++) {
    L[i][i] = 1;
    for (int j = 0; j < i; j++) {
      L[i][j] = A[i][j];
    }

    U[i][i] = A[i][i];
    for (int j = i + 1; j < n; j++) {
      U[i][j] = A[i][j];
    }
  }

  return std::tuple(L, U, perm);
}

template <class T>
std::vector<T> solve(const Matrix<T>& A, const vector<T>& b) {
  auto [L, U, perm] = decompose(A);
  return solve(L, U, perm, b);
}

/**
 * @brief Given L, U, perm such that LUP is a LUP-decomposition of A,
 * finds a vector x such that PAx = Pb
 * Complexity: O(n^3);
 *
 * @return std::vector<T>
 */
template <class T>
std::vector<T> solve(const Matrix<T>& L, const Matrix<T>& U,
                     const vector<T>& perm, const vector<T>& b) {
  const int n = L.num_rows();

  // Ly = PA
  vector<T> y(n);
  for (int i = 0; i < n; i++) {
    y[i] = b[perm[i]];
    for (int j = 0; j < i; j++) {
      y[i] -= L[i][j] * y[j];
    }
  }

  // LUx = PA
  vector<T> x(n);
  for (int i = n - 1; i >= 0; i--) {
    x[i] = y[i];
    for (int j = n - 1; j > i; j--) {
      x[i] -= U[i][j] * x[j];
    }
    x[i] /= U[i][i];
  }
  return x;
}

};  // namespace LUP

template <class T>
T Matrix<T>::determinant() {
  assert(is_nonsingular());
  auto [L, U, _] = LUP::decompose(*this);
  T det = 1;
  for (int i = 0; i < num_rows(); i++) {
    det *= L[i][i];
    det *= U[i][i];
  }
  return det;
}

template <class T>
Matrix<T> Matrix<T>::inverse() {
  auto [L, U, perm] = LUP::decompose(*this);

  const int n = this->num_rows();
  Matrix<T> X(n, n);
  vector<T> b(n);
  for (int i = 0; i < n; i++) {
    if (i > 0) b[i - 1] = 0;
    b[i] = 1;
    vector<T> x = LUP::solve(L, U, perm, b);
    for (int j = 0; j < n; j++) {
      X[i][j] = x[j];
    }
  }
  return X;
}

template <class T>
bool Matrix<T>::is_singular() {
  try {
    LUP::decompose((*this));
  } catch (SingularMatrixError e) {
    return true;
  }
  return false;
}

template <class T>
bool Matrix<T>::is_nonsingular() {
  return not is_singular();
}

template <class T>
T det(const Matrix<T>& M) {
  return M.determinant();
}
