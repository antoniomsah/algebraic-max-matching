#pragma once

#include <vector>

using std::pair;
using std::vector;

class IAlgorithmStrategy {
 public:
  virtual ~IAlgorithmStrategy() = default;
  virtual vector<pair<int, int>> solve(
      const vector<vector<int>> &graph,
      const vector<pair<int, int>> &edges) const = 0;
};