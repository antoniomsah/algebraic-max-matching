#pragma once

#include <algorithm>
#include <iostream>
#include <memory>
#include <set>
#include <vector>

#include "algorithm-strategy-interface.hpp"
#include "rank2/rank2-algorithm.hpp"
#include "simple/simple-algorithm.hpp"

using std::cout;
using std::pair;
using std::vector;

/**
 * @brief Class that solves perfect/maximum matching in general graphs.
 */
class MatchingSolver {
  vector<vector<int>> graph;
  vector<pair<int, int>> edges;

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

  /**
   * @brief MatchingSolver's constructor.
   * Default strategy is the naive one.
   * @param graph vector<vector<int>> that represents a graph
   * @param edges vector<pair<int,int>> of edges
   * @param strategy_id id of the selected strategy (default is 0).
   */
  MatchingSolver(vector<vector<int>> graph, vector<pair<int, int>> edges,
                 int strategy_id = 0)
      : graph(graph), edges(edges) {
    set_strategy(strategy_id);
  }

  /**
   * @brief Given a graph prints the maximum matching size.
   * If perfect=true, it prints "YES" iff the matching is perfect.
   * If print_matching=true, it prints the matching that was found.
   */
  void solve(bool perfect = true, bool print_matching = false) {
    vector<pair<int, int>> matching;

    size_t matching_size = 0;
    for (int it = 0; it < MAX_IT; it++) {
      matching = algorithmStrategy->solve(graph, edges);
      if (validate(matching, perfect)) {
        matching_size = max(matching_size, matching.size());
        if (2 * matching_size == graph.size()) break;
      }
    }

    if (not perfect or 2 * matching_size == graph.size()) {
      if (perfect) {
        cout << "YES\n";
      } else {
        cout << matching_size << '\n';
      }

      if (print_matching) {
        for (const auto& [u, v] : matching) {
          cout << u << ' ' << v << '\n';
        }
      }
    } else {
      cout << "NO\n";
    }
  }

  /**
   * @brief Given a adjacency matrix (i.e., 1 if edge[u,v] exists and 0,
   * otherwise), checks if a matching is valid (and if wants perfect matching,
   * checks that as well).
   */
  bool validate(const vector<pair<int, int>>& matching,
                bool want_perfect = false) {
    set<int> matched_vertices;
    for (auto [u, v] : matching) {
      if (graph[u][v]) {
        if (matched_vertices.find(u) != matched_vertices.end()) return false;
        if (matched_vertices.find(v) != matched_vertices.end()) return false;
        matched_vertices.insert(u);
        matched_vertices.insert(v);
      } else {
        return false;
      }
    }
    if (want_perfect and matched_vertices.size() != graph.size()) return false;
    return true;
  }
};

array<std::shared_ptr<IAlgorithmStrategy>, NUM_STRATEGIES>
    MatchingSolver::strategies;