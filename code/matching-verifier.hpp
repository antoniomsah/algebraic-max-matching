#pragma once

#include <set>
#include <vector>

namespace MatchingValidator {

using std::pair;
using std::set;
using std::vector;

/**
 * @brief Given a adjacency matrix (i.e., 1 if edge[u,v] exists and 0,
 * otherwise), checks if a matching is valid (and if wants perfect matching,
 * checks that as well).
 */
bool validate(const vector<vector<int>> g,
              const vector<pair<int, int>>& matching,
              bool want_perfect = false) {
  set<int> matched_vertices;
  for (auto [u, v] : matching) {
    if (g[u][v]) {
      if (matched_vertices.find(u) != matched_vertices.end()) return false;
      if (matched_vertices.find(v) != matched_vertices.end()) return false;
      matched_vertices.insert(u);
      matched_vertices.insert(v);
    } else {
      return false;
    }
  }
  if (want_perfect and matched_vertices.size() != g.size()) return false;
  return true;
}

}  // namespace MatchingValidator