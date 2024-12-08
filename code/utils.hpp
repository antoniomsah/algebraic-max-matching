#pragma once

#include "vector"
#include "array"

#include "classes/tutte-matrix.hpp"

using std::array;
using std::vector;

array<vector<int>, 2> DivideInTwo(const vector<int>& A) {
    const int n = A.size();
    array<vector<int>, 2> S;
    S[0].reserve(n / 2);
    S[1].reserve(n - n / 2);
    for (int i = 0; i < n; i++) {
        S[i >= n / 2].push_back(A[i]);
    }
    return S;
}

// Unite unites two vectors A and B, i.e. A \cup B.
vector<int> Unite(const vector<int>& A, const vector<int> &B) {
    vector<int> AB(A.size() + B.size());
    for (int i = 0; i < AB.size(); i++) {
        if (i < A.size()) {
            AB[i] = A[i];
        } else {
            AB[i] = B[i - A.size()];
        }
    }
    return AB;
}

// Unite unite four matrices.
template <typename T>
Matrix<T> unite(
  Matrix<T> top_left, Matrix<T> top_right, 
  Matrix<T> bot_left, Matrix<T> bot_right
) {
  int n = top_left.numRows() + bot_left.numRows(),
      m = top_left.numColumns() + top_right.numColumns();

  Matrix<T> A(n, m);
  for (int i = 0; i < n / 2; i++) {
    for (int j = 0; j < m / 2; j++) {
      A(i, j) = top_left(i, j);
    }
  }

  for (int i = 0; i < n / 2; i++) {
    for (int j = m / 2; j < m; j++) {
      A(i, j) = top_right(i, j - m / 2);
    }
  }

  for (int i = n / 2; i < n; i++) {
    for (int j = 0; j < m / 2; j++) {
      A(i, j) = bot_left(i - n / 2, j);
    }
  }

  for (int i = n / 2; i < n; i++) {
    for (int j = m / 2; j < m; j++) {
      A(i, j) = bot_right(i - n / 2, j - m / 2);
    }
  }
  return A;
}

template <typename T>
array<Matrix<T>, 4> split(const Matrix<T> &M) {
    assert(M.isSquare());
    size_t n = M.numRows();
    Matrix<T> A(n/2, n/2), B(n/2, n/2), C(n/2, n/2), D(n/2, n/2);
    for (int i = 0; i < n/2; i++) {
        for (int j = 0; j < n/2; j++) {
            A(i, j) = M(i, j);
            B(i, j) = M(i, n/2 + j);
            C(i, j) = M(n/2 + i, j);
            D(i, j) = M(n/2 + i, n/2 + j);
        }
    }
    return {A, B, C, D};
}

template <const int P>
TutteMatrix<P> Identity(size_t n) {
    TutteMatrix<P> I(n);
    for (int i = 0; i < n; i++) {
        I(i, i) = 1;
    }
    return I;
}

// Performs a RankTwoUpdate on TutteMatrices.
// S must have size two.
template <const int P>
TutteMatrix<P> rankTwoUpdate(const vector<int>& S, 
                            const TutteMatrix<P>& T,
                            const TutteMatrix<P>& N) {
    assert(S.size() == 2);

    const int u = S[0], v = S[1];
    TutteMatrix<P> TN(2);
    TN(0, 0) = (T(u, v) * N(u, v) + 1).inv();
    TN(1, 1) = (T(u, v) * N(u, v) + 1).inv();

    return N + N('*', S) * TN * T(S, S) * N(S, '*');
}