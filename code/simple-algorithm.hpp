#pragma once

#include <map>
#include <vector>

#include "tutte-matrix.hpp"

namespace SimpleAlgorithm {

const int MAX_IT = 100;
const int MOD = 998244353;

using namespace std;

void solve(const int& n, const vector<pair<int, int>>& edges) {
  const int m = edges.size();
  map<vector<pair<int, int>>, int> res;
  for (int it = 0; it < MAX_IT; it++) {
    vector<bool> is_matching_edge(m);
    TutteMatrix<MOD> T(n, edges);
    for (int i = 0; i < m; i++) {
      T.remove_edge(i);
      if (not T.has_perfect_matching()) {
        is_matching_edge[i] = true;
        T.add_edge(i);
      }
    }

    vector<pair<int, int>> matching;
    for (int i = 0; i < m; i++) {
      if (is_matching_edge[i]) matching.push_back(edges[i]);
    }
    res[matching] += 1;
  }

  for (auto [matching, ocr] : res) {
    cout << "ocurred " << ocr << " times: ";
    for (auto [u, v] : matching) cout << "(" << u << ", " << v << ") ";
  }
  cout << '\n';
}
}  // namespace SimpleAlgorithm