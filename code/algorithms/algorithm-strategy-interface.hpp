#pragma once

#include <vector>
#include "../classes/graph.hpp"

class IAlgorithmStrategy {
 public:
  virtual ~IAlgorithmStrategy() = default;
  virtual std::vector<std::pair<int, int>> solve(const Graph& graph) const = 0;
  virtual std::string name() const = 0;
};
