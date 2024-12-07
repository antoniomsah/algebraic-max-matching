#pragma once

#include <vector>
#include "../classes/graph.hpp"

using std::pair;
using std::vector;

class IAlgorithmStrategy {
 public:
  virtual ~IAlgorithmStrategy() = default;
  virtual vector<pair<int, int>> solve(const Graph& graph) const = 0;
};