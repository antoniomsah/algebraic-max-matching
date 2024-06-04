#include <iostream>
#include <map>
#include <vector>

#include "matrix.hpp"
#include "tutte-matrix.hpp"

using namespace std;

template <const int MOD>
struct Mint {
  long long x;

  Mint() : x(0) {}
  Mint(int _x) : x(_x % MOD < 0 ? _x % MOD + MOD : _x % MOD) {}

  void operator+=(Mint rhs) {
    x += rhs.x;
    if (x >= MOD) x -= MOD;
  }

  void operator-=(Mint rhs) {
    x -= rhs.x;
    if (x < 0) x += MOD;
  }

  void operator*=(Mint rhs) {
    x *= rhs.x;
    x %= MOD;
  }

  void operator/=(Mint rhs) { *this *= rhs.inv(); }

  Mint operator+(Mint rhs) {
    Mint res = *this;
    res += rhs;
    return res;
  }

  Mint operator-(Mint rhs) {
    Mint res = *this;
    res -= rhs;
    return res;
  }

  Mint operator*(Mint rhs) {
    Mint res = *this;
    res *= rhs;
    return res;
  }

  Mint operator/(Mint rhs) {
    Mint res = *this;
    res /= rhs;
    return res;
  }

  Mint inv() { return this->pow(MOD - 2); }

  Mint pow(int e) {
    Mint res(1);
    for (Mint p = *this; e > 0; e /= 2, p *= p)
      if (e % 2) res *= p;
    return res;
  }

  bool operator<(int y) { return x < y; }
  bool operator>(int y) { return x > y; }
  bool operator==(int y) { return x == y; }

  bool operator<(Mint y) { return x < y.x; }
  bool operator>(Mint y) { return x > y.x; }
  bool operator==(Mint y) { return x == y.x; }

  friend ostream& operator<<(ostream& os, Mint x) { return (os << x.x); }
};

const int MOD = 998244353;

const int IT_SIZE = 100;

int main() {
  int n, m;
  cin >> n >> m;
  vector<pair<int, int>> edges(m);
  for (int i = 0; i < m; i++) {
    int u, v;
    cin >> u >> v;
    edges[i] = {u - 1, v - 1};
  }

  map<vector<pair<int, int>>, int> res;
  for (int it = 0; it < IT_SIZE; it++) {
    TutteMatrix<Mint<MOD>, MOD> T(n, edges);
    vector<bool> is_essential(m);
    for (int e = 0; e < m; e++) {
      T.remove_edge(e);
      if (not T.has_perfect_matching()) {
        is_essential[e] = true;
        T.add_edge(e);
      }
    }

    vector<pair<int, int>> matching;
    for (int e = 0; e < edges.size(); e++) {
      if (is_essential[e]) matching.push_back(edges[e]);
    }
    res[matching] += 1;
  }

  int id = 0;
  for (auto [matching, ocr] : res) {
    cout << id++ << " (" << ocr << ")\n";
    for (auto [u, v] : matching) cout << u << ' ' << v << '\n';
  }

  return 0;
}
