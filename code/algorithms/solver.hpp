#pragma once

#include <algorithm>
#include <iostream>
#include <memory>
#include <set>
#include <vector>

#include "algorithm-strategy-interface.hpp"
#include "rank2/rank2-algorithm.hpp"
#include "simple/simple-algorithm.hpp"
#include "harvey/harvey-algorithm.hpp"

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
    strategies[2] = std::make_shared<HarveyAlgorithmStrategy>();
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
  void PerfectMatching(bool print_matching = false) {
    vector<pair<int, int>> matching;

    size_t matching_size = 0;
    for (int it = 0; it < MAX_IT; it++) {
      matching = algorithmStrategy->solve(graph, edges);
      if (validate(matching, true /* wantsPerfectMatching */)) {
        matching_size = max(matching_size, matching.size());
        if (2 * matching_size == graph.size()) break;
      }
    }

    if (2 * matching_size != graph.size()) {
      cout<<"NO\n";
    }
    cout << "YES\n";
    if (print_matching) {
      for (const auto& [u, v] : matching) {
        cout << u << ' ' << v << '\n';
      }
    }
  }

  void MaximumMatching() {
    const int n = graph.size();

    int sz = get_maximum_matching_size();
    int unmatched_vertices = n - 2 * sz;
    vector<vector<int>> G(n + unmatched_vertices, vector<int>(n + unmatched_vertices));
    vector<pair<int, int>> E = edges;
    for (int i = 0; i < unmatched_vertices; i++) {
      int u = n + i; // Index of the new vertex.
      for (int v = 0; v < n; v++) {
        E.emplace_back(u, v);
      }
    }

    for (auto [u, v] : E) {
      G[u][v] = G[v][u] = 1;
    }

    // Augmented graph G has been created.
    vector<pair<int, int>> matching;
    for (int it = 0; it < MAX_IT; it++) {
      matching = algorithmStrategy->solve(G, E);
      if (validate(matching, false /* wantsPerfectMatching */) and matching.size() == sz) {
        break;
      }
    }

    cout << sz << '\n';
    for (auto [u, v] : matching) {
      if (u >= n or v >= n) continue; // Not an original vertex.
      cout << u << ' ' << v << '\n';
    }
  }

  int get_maximum_matching_size() {
    TutteMatrix<MOD> T = build_tutte_matrix<MOD>(graph.size(), edges);
    return T.rank() / 2;
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
      if (!graph[u][v]) {
        return false;
      } 
      if (matched_vertices.find(u) != matched_vertices.end()) return false;
      if (matched_vertices.find(v) != matched_vertices.end()) return false;
      matched_vertices.insert(u);
      matched_vertices.insert(v);
    }
    if (want_perfect and matched_vertices.size() != graph.size()) return false;
    return true;
  }
};

array<std::shared_ptr<IAlgorithmStrategy>, NUM_STRATEGIES>
    MatchingSolver::strategies;