#include <iostream>

#include "simple-algorithm.hpp"

using namespace std;

const int MOD = 998244353;
const int IT_SIZE = 1000;

int main() {
  int n, m;
  cin >> n >> m;
  vector<pair<int, int>> edges(m);
  for (int i = 0; i < m; i++) {
    int u, v;
    cin >> u >> v;
    edges[i] = {u, v};
  }

  SimpleAlgorithm::solve(n, edges);

  return 0;
}
