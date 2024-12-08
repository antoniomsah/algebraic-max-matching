#include <bits/stdc++.h>

#include "../../classes/matrix.hpp"
#include "../../classes/modular-integer.hpp"
#include "../../utils.hpp"

using namespace std;

const int P = 998244353;

signed main() {
  int n;
  cin >> n;
  Matrix<modular_int<P>> M(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      int x;
      cin >> x;
      M(i, j) = x;
    }
  }

  if (M.isNonSingular()) {
    cout << "YES\n";
    Matrix<modular_int<P>> N = M.inverse();
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        cout << N(i, j) << ' ';
      }
      cout << '\n';
    }
  } else {
    cout << "NO\n";
  }
  return 0;
}
