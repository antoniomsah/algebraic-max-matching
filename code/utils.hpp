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