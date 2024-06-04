#pragma once

#include <map>
#include <vector>

#include "matching-verifier.hpp"
#include "tutte-matrix.hpp"

namespace SimpleAlgorithm {

const int MAX_IT = 100;
const int MOD = 998244353;

using namespace std;

/**
 * @brief Finds a perfect matching. Complexity: O(n^{\omega+2}).
 *
 * @param n number of vertices
 * @param edges std::vector of edges
 */
vector<pair<int, int>> solve(const vector<vector<int>> graph,
                             const vector<pair<int, int>>& edges,
                             bool test = false) {
  const int n = graph.size(), m = edges.size();

  map<vector<pair<int, int>>, int> found;
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
    found[matching] += 1;
  }

  vector<pair<int, int>> res;
  for (auto [matching, ocr] : found) {
    if (MatchingValidator::validate(graph, matching, true)) res = matching;
  }

  if (test) {  // for testing purposes
    int correct_outputs = 0;
    for (auto [matching, ocr] : found) {
      if (MatchingValidator::validate(graph, matching, true))
        correct_outputs += ocr;
    }

    cout << "found " << found.size() << " different outputs and generated "
         << correct_outputs << " correct outputs out of " << MAX_IT
         << " iterations (prob = " << double(correct_outputs) / MAX_IT << ")\n";
  }
  return res;
}

}  // namespace SimpleAlgorithm