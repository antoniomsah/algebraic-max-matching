#pragma once

#include "algorithm-strategy-interface.hpp"

using namespace std;

/** 
 * Implements the Blossom algorithm.
 * Solves maximum matching in O(n^2m).
 *
 * This implementation is based on the Capstone project from
 * Giovana Gomes Delfino, see https://bccdev.ime.usp.br/tccs/2017/gigd/.
 **/
class BlossomAlgorithm {
  int n, m;

  struct Edge {
      int v, next;
      Edge(int _v = -1, int _next = -1) : v(_v), next(_next) {}
  };

  vector<int> adj;
  vector<Edge> edges;

  int qh, qt;
  vector<int> match, q, p, base, inQ, inB;

  void add_edge(int u, int v) {
    edges.emplace_back(v, adj[u]);
    adj[u] = edges.size() - 1;

    edges.emplace_back(u, adj[v]);
    adj[v] = edges.size() - 1;
  }

  int LCA(int root, int u, int v) {
    vector<bool> inP(n, false);
    while (42) {
      u = base[u];
      inP[u] = true;
      if (u == root) break;
      u = p[match[u]];
    }

    while (42) {
      v = base[v];
      if (inP[v]) return v;
      else v = p[match[v]];
    }
  }

  void mark(int lca, int u) {
    while (base[u] != lca) {
      int v = match[u];
      inB[base[u]] = inB[base[v]] = true;
      u = p[v]; 
      if (base[u] != lca) p[u] = v;
    }
  }

  void contract(int s, int u, int v) {
    int lca = LCA(s, u, v);
    fill(inB.begin(), inB.end(), 0);
    mark(lca, u);
    mark(lca, v);
    
    if (base[u] != lca) p[u] = v;
    if (base[v] != lca) p[v] = u;

    for (int u = 0; u < n; u++) {
      if (!inB[base[u]]) continue;
      base[u] = lca;
      if (!inQ[u]) { 
          q[++qt] = u;
          inQ[u] = true;
      }
    }
  }

  int find_augmenting_path(int s) {
    fill(inQ.begin(), inQ.end(), 0);
    fill(p.begin(), p.end(), -1);
    iota(base.begin(), base.end(), 0);

    qh = qt = 0;
    q[0] = s;
    inQ[s] = true;

    while (qh <= qt) {
        int u = q[qh++];
        for (int e = adj[u]; e != -1; e = edges[e].next) {
            int v = edges[e].v;
            if (base[u] == base[v] || match[u] == v) continue;

            if ((v == s) || (match[v] != -1 && p[match[v]] != -1)) {
              contract(s, u, v);
              continue;
            }

            if (p[v] == -1) {
                p[v] = u;
                if (match[v] == -1) return v;

                if (!inQ[match[v]]) {
                    q[++qt] = match[v];
                    inQ[match[v]] = true;
                }
            }
        }
    }
    return -1;
  }

  int augment(int s, int t) {
    int u = t, v, w;
    while (u != -1) {
        v = p[u];
        w = match[v];
        match[v] = u;
        match[u] = v;
        u = w;
    }
    return t != -1;
  }

public: 
  BlossomAlgorithm(const Graph& _G) : n(V(_G).size()), m(E(_G).size()) {
    match.resize(n);
    q.resize(n);
    p.resize(n);
    base.resize(n);
    inQ.resize(n);
    inB.resize(n);
    adj.resize(n, -1);
    for (auto [u, v]: E(_G)) {
      add_edge(u, v);
    }
  }

  vector<pair<int, int>> solve() {
    int match_new = 0;
    fill(match.begin(), match.end(), -1);
    for (int u = 0; u < n; u++) {
      if (match[u] == -1) {
        match_new += augment(u, find_augmenting_path(u));
      }
    }

    vector<pair<int, int>> M;
    for (int u = 0; u < n; u++) {
      if (u < match[u]) M.emplace_back(match[u], u);
    }
    return M;
  }
};

/**
 * BlossomAlgorithmStrategy uses an algorithm based on Edmonds-Blossom algorithm
 * to solve Maximum Matching in O(n^2m).
 */
class BlossomAlgorithmStrategy : public IAlgorithmStrategy {
 public:
  vector<pair<int, int>> solve(const Graph& G) const override {
    return BlossomAlgorithm(G).solve();
  }

  string name() const override {
    return "Blossom";
  }
};
