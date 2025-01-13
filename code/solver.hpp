#pragma once

#include <algorithm>
#include <iostream>
#include <memory>
#include <set>
#include <vector>

#include "algorithms/algorithm-strategy-interface.hpp"
#include "algorithms/rank2-algorithm.hpp"
#include "algorithms/naive-algorithm.hpp"
#include "algorithms/harvey-algorithm.hpp"
#include "algorithms/blossom-algorithm.hpp"
#include "classes/graph.hpp"
#include "config.hpp"

using std::cout;
using std::pair;
using std::vector;

/**
 * Class that solves perfect/maximum matching in general graphs.
 */
class MatchingSolver {
  Graph G;

  std::shared_ptr<IAlgorithmStrategy> algorithmStrategy;

 public:
  // Algorithm strategies
  static array<std::shared_ptr<IAlgorithmStrategy>, NUM_STRATEGIES> strategies;

  static void initialize() {
    strategies[0] = std::make_shared<NaiveAlgorithmStrategy>();
    strategies[1] = std::make_shared<RankTwoAlgorithmStrategy>();
    strategies[2] = std::make_shared<HarveyAlgorithmStrategy>();
    strategies[3] = std::make_shared<BlossomAlgorithmStrategy>();
  }

  // setStrategy changes the algorithm strategy.
  // The default values are 
  //  - 0: NaiveAlgorithm // O(n^{\omega + 2})
  //  - 1: RankTwoAlgorithm // O(n^4)
  //  - 2: HarveyAlgorithm // O(n^\omega)
  void setStrategy(const int& id) {
    if (id >= NUM_STRATEGIES) {
      std::cout << "Invalid strategy index, number of strategies is: "
                << NUM_STRATEGIES << "\n";
      return;
    } 
    algorithmStrategy = strategies[id];
  }

  string algorithmName() {
    return algorithmStrategy->name();
  }

  /**
   * Default strategy is the naive one.
   * @param graph a Graph class;
   * @param strategy_id id of the selected strategy (default is 0).
   */
  MatchingSolver(Graph graph, int strategy_id = 0) : G(graph) {
    setStrategy(strategy_id);
  }

  /**
   * PerfectMatching finds a perfect matching in the underlying graph of MatchingSolver.
   * @returns A perfect matching, if one exists. Else, returns an empty set.
   */
  vector<pair<int, int>> PerfectMatching();

  /**
   * MaximumMatching finds a maximum matching in the underlying graph of MatchingSolver.
   * It prints any maximum matching from the graph.
   * @returns A maximum matching.
   */
  vector<pair<int, int>> MaximumMatching();
};

array<std::shared_ptr<IAlgorithmStrategy>, NUM_STRATEGIES>
    MatchingSolver::strategies;

vector<pair<int, int>> MatchingSolver::PerfectMatching() {
    vector<pair<int, int>> matching;

    size_t matching_size = 0;
    for (int it = 0; it < MAX_IT; it++) {
      matching = algorithmStrategy->solve(G);
      if (G.hasMatching(matching)) {
        matching_size = max(matching_size, matching.size());
        if (2 * matching_size == V(G).size()) break;
      }
    }

    if (2 * matching_size != V(G).size()) {
      return {};
    }
    return matching;
}

vector<pair<int, int>> MatchingSolver::MaximumMatching() {
  const int n = V(G).size();

  int unmatched = n - 2 * G.MatchingNumber();
  Graph aG(n + unmatched); // Augmented Graph.

  for (const auto& [u, v] : E(G)) {
    aG.addEdge(u, v);
  }

  for (int u = n; u < n + unmatched; u++) {
    for (int v = 0; v < n; v++) {
      aG.addEdge(u, v);
    }
  }

  // Augmented graph aG has been created.
  vector<pair<int, int>> aM; // Matching in the Augmented Graph.
  for (int it = 0; it < MAX_IT; it++) {
    aM = algorithmStrategy->solve(aG);
    if (aG.hasMatching(aM)) {
      break;
    }
  }

  vector<pair<int, int>> M;
  for (const auto &[u, v] : aM) {
    if (G.isAdj(u, v)) { // This edge exists in the original graph.
      M.emplace_back(u, v);
    }
  }
  return M;
}
