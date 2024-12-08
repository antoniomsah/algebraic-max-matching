#pragma once

#include "vector"
#include "array"

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