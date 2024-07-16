#pragma once

#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

#include "algorithm-strategy-interface.hpp"
#include "rank2/rank2-algorithm.hpp"
#include "simple/simple-algorithm.hpp"

using std::cout;
using std::pair;
using std::vector;

class MatchingSolver {
  std::shared_ptr<IAlgorithmStrategy> algorithmStrategy;

 public:
  // Algorithm strategies
  static array<std::shared_ptr<IAlgorithmStrategy>, NUM_STRATEGIES> strategies;

  static void initialize() {
    strategies[0] = std::make_shared<SimpleAlgorithmStrategy>();
    strategies[1] = std::make_shared<Rank2AlgorithmStrategy>();
  }

  void set_strategy(const size_t& id) {
    if (id < NUM_STRATEGIES) {
      algorithmStrategy = strategies[id];
    } else {
      std::cout << "Invalid strategy index, number of strategies is: "
                << NUM_STRATEGIES << "\n";
    }
  }

  MatchingSolver() { set_strategy(0); }

  void solve(const vector<vector<int>>& graph,
             const vector<pair<int, int>>& edges) {
    vector<pair<int, int>> matching = algorithmStrategy->solve(graph, edges);
    cout << matching.size() << '\n';
    if (matching.size() > 1 or matching[0] != pair(-1, -1)) {
      cout << "YES\n";
      for (const auto& [u, v] : matching) {
        cout << u << ' ' << v << '\n';
      }
    } else {
      cout << "NO\n";
    }
  }
};