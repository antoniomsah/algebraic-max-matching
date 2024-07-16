#include <bits/stdc++.h>

#include "../../matrix/matrix.hpp"
#include "../../matrix/modular-integer.hpp"

using namespace std;

const int MOD = 998244353;

signed main() {
  int n;
  cin >> n;
  Matrix<modular_int<MOD>> M(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      int x;
      cin >> x;
      M(i, j) = x;
    }
  }

  if (M.is_nonsingular()) {
    cout << "YES\n";
    Matrix<modular_int<MOD>> N = M.inverse();
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
