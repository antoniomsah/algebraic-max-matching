#pragma once

#include <algorithm>
#include <iostream>
#include <memory>
#include <set>
#include <vector>

#include "algorithm-strategy-interface.hpp"
#include "rank2/rank2-algorithm.hpp"
#include "naive/naive-algorithm.hpp"
#include "harvey/harvey-algorithm.hpp"
#include "../classes/graph.hpp"
#include "../config.hpp"

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
    strategies[1] = std::make_shared<Rank2AlgorithmStrategy>();
    strategies[2] = std::make_shared<HarveyAlgorithmStrategy>();
  }

  // setStrategy changes the algorithm strategy.
  // The default values are 
  //  * 0 := NaiveAlgorithm // O(n^{\omega + 2})
  //  * 1 := RankTwoAlgorithm // O(n^4)
  //  * 2 := HarveyAlgorithm // O(n^\omega)
  void setStrategy(const size_t& id) {
    if (id >= NUM_STRATEGIES) {
      std::cout << "Invalid strategy index, number of strategies is: "
                << NUM_STRATEGIES << "\n";
      return;
    } 
    algorithmStrategy = strategies[id];
  }

  /**
   * @brief MatchingSolver's constructor.
   * Default strategy is the naive one.
   * @param graph vector<vector<int>> that represents a graph
   * @param edges vector<pair<int,int>> of edges
   * @param strategy_id id of the selected strategy (default is 0).
   */
  MatchingSolver(Graph graph, int strategy_id = 0) : G(graph) {
    setStrategy(strategy_id);
  }

  /**
   * Finds a perfect matching in the underlying graph of MatchingSolver.
   * It prints any perfect matching from the graph.
   */
  void PerfectMatching();

  /**
   * Finds a maximum matching in the underlying graph of MatchingSolver.
   * It prints any maximum matching from the graph.
   */
  void MaximumMatching();

  int get_maximum_matching_size() {
    TutteMatrix<MOD> T = GetTutteMatrix(G);
    return T.rank() / 2;
  }

  /**
   * Validates checks if a matching is a valid in the underlying graph of MatchingSolver.
   */
  bool validate(const vector<pair<int, int>>& matching) {
    set<int> matched_vertices;
    for (auto [u, v] : matching) {
      if (!G.isAdj(u, v)) { // This edge DOES not exist in the original graph.
        return false;
      } 
      if (matched_vertices.find(u) != matched_vertices.end()) return false;
      if (matched_vertices.find(v) != matched_vertices.end()) return false;
      matched_vertices.insert(u);
      matched_vertices.insert(v);
    }
    return true;
  }
};

array<std::shared_ptr<IAlgorithmStrategy>, NUM_STRATEGIES>
    MatchingSolver::strategies;

void MatchingSolver::PerfectMatching() {
    vector<pair<int, int>> matching;

    size_t matching_size = 0;
    for (int it = 0; it < MAX_IT; it++) {
      matching = algorithmStrategy->solve(G);
      if (validate(matching)) {
        matching_size = max(matching_size, matching.size());
        if (2 * matching_size == G.size()) break;
      }
    }

    if (2 * matching_size != G.size()) {
      cout<<"NO\n";
      return;
    }

    cout << "YES\n";
    //for (const auto& [u, v] : matching) {
      //cout << u << ' ' << v << '\n';
    //}
}

void MatchingSolver::MaximumMatching() {
  const int n = G.size();

  int sz = get_maximum_matching_size();
  int unmatched_vertices = n - 2 * sz;
  // vector<vector<int>> G(n + unmatched_vertices, vector<int>(n + unmatched_vertices));
  // vector<pair<int, int>> E = E(G);
  // for (int i = 0; i < unmatched_vertices; i++) {
    // int u = n + i; // Index of the new vertex.
    // for (int v = 0; v < n; v++) {
      // E.emplace_back(u, v);
    // }
  // }

  // for (auto [u, v] : E) {
    // G[u][v] = G[v][u] = 1;
  // }

  // // Augmented graph G has been created.
  // vector<pair<int, int>> matching;
  // for (int it = 0; it < MAX_IT; it++) {
    // matching = algorithmStrategy->solve(G, E);
    // if (validate(matching) and matching.size() == sz) {
      // break;
    // }
  // }

  // cout << sz << '\n';
  // for (auto [u, v] : matching) {
    // if (u >= n or v >= n) continue; // Not an original vertex.
    // cout << u << ' ' << v << '\n';
  // }
}