#pragma once

#include "algorithm-strategy-interface.hpp"
#include <queue>

using namespace std;

/** 
 * Implements the Blossom algorithm.
 * Solves maximum matching in O(nm).
 *
 * The implementation was based on: https://codeforces.com/blog/entry/92339.
 **/
class BlossomAlgorithm {
  int n, sz;
  vector<int> mate, // mate[u] is the vertex matched to u.
              id,   // id[u] is the blossom containing u.
              p,    // p[u] is the parent of u.
              vis;  // vis[u] is 0 if not visited, 1 if even depth and 2 if odd depth.
  vector<vector<int>> blossom; // blossom[b] is the list of vertices contained in blossom 'b'.
  Matrix<int> g;

  /**
   * Prepares the iteration.
   **/
  queue<int> init() {
    queue<int> q;
    fill(vis.begin(), vis.end(), 0);
    iota(id.begin(), id.end(), 0);
    for (int v = 0; v < n; v++) if (mate[v] == -1) {
      // Exposed vertex.
      q.push(v);
      p[v] = v;
      vis[v] = 1;
    }
    return q;
  }

  // addEdge adds the edge {u, v} to 'g'.
  void addEdge(int u, int v) {
    assert(u < sz && v < sz);
    g(u, v) = v;
    g(v, u) = u;
  }

  // removeEdge removes edge {u, v} from 'g'.
  void removeEdge(int u, int v) {
    assert(u < sz && v < sz);
    g(u, v) = -1;
    g(v, u) = -1;
  }

  // isAdj checks if two vertices are adjacent.
  bool isAdj(int u, int v) {
    assert(u < sz && v < sz);
    return g(u, v) != -1;
  }

  // match matches vertices 'u' and 'v'.
  void match(int u, int v) {
    assert(u < sz && v < sz);
    removeEdge(u, v);
    mate[u] = v;
    mate[v] = u;
  }

  // contract contracts a blossom.
  void contract(int c, int u, int v, vector<int> &vx, vector<int> &vy) {
    blossom[c].clear();
    int r = vx.back();
    while (!vx.empty() && !vy.empty() && vx.back() == vy.back()) {
      r = vx.back();
      vx.pop_back();
      vy.pop_back();
    }
    blossom[c].push_back(r);
    blossom[c].insert(blossom[c].end(), vx.rbegin(), vx.rend());
    blossom[c].insert(blossom[c].end(), vy.begin(), vy.end());
    for (int i = 0; i <= c; i++) {
      removeEdge(i, c);
    }
    for (int b: blossom[c]) {
      id[b] = c;
      for (int i = 0; i < c; i++) {
        if (!isAdj(b, i)) continue;
        g(c, i) = b;
        g(i, c) = g(i, b);
      }
    }
  }

  // trace returns the path from 'u' to the root.
  vector<int> trace(int u) {
    vector<int> vx;
    while (true) {
      while (id[u] != u) u = id[u];
      if (!vx.empty() && vx.back() == u) break;
      vx.push_back(u);
      u = p[u];
    }
    return vx;
  }

  // lift liftes a path from the contracted graph to the original graph.
  vector<int> lift(vector<int> &vx) {
    vector<int> A;
    while (vx.size() >= 2) {
      int z = vx.back(); vx.pop_back();
      if (z < n) {
        A.push_back(z);
        continue;
      }

      int w = vx.back();
      int i = (A.size() % 2 == 0 ? find(blossom[z].begin(), blossom[z].end(), g(z, w)) - blossom[z].begin() : 0);
      int j = (A.size() % 2 == 1 ? find(blossom[z].begin(), blossom[z].end(), g(z, A.back())) - blossom[z].begin() : 0);
      int k = blossom[z].size();
      int dif = (A.size() % 2 == 0 ? i % 2 == 1 : j % 2 == 0) ? 1 : k - 1;
      while (i != j) {
          vx.push_back(blossom[z][i]);
          i = (i + dif) % k;
      }
      vx.push_back(blossom[z][i]);
    }
    return A;
  }

public:
  BlossomAlgorithm(const Graph& _G) : n(V(_G).size()), sz(n + n/2), g(sz, sz, -1) {
    mate.resize(n, -1);
    id.resize(sz);
    p.resize(sz);
    vis.resize(sz);
    blossom.resize(sz);
    for (const auto& [u, v]: E(_G)) {
      addEdge(u, v);
    }
  }

  vector<pair<int, int>> solve() {
    bool aug = true;
    while (aug) {
      queue<int> q = init();
      int c = n;
      aug = false;
      while (!q.empty() && !aug) {
        int u = q.front(); q.pop();
        if (id[u] != u) continue; // Not exposed.
        for (int v = 0; v < c; v++) {
          if (id[v] != v || !isAdj(u, v)) {
            continue;
          }

          if (vis[v] == 0) {
            p[v] = u;
            vis[v] = 2;
            int m = mate[v]; // Matched vertex to v.
            p[m] = v;
            vis[m] = 1;
            q.push(m);
            continue;
          }

          if (vis[v] == 2) {
            continue;
          }

          vector<int> vx = trace(u), vy = trace(v);
          if (vx.back() == vy.back()) {
            contract(c, u, v, vx, vy);
            q.push(c);
            p[c] = p[blossom[c][0]];
            vis[c] = 1;
            c++;
            break;
          }

          aug = true;
          vx.insert(vx.begin(), v);
          vy.insert(vy.begin(), u);
          vector<int> A = lift(vx), B = lift(vy);
          A.insert(A.end(), B.rbegin(), B.rend());
          for (size_t i = 0; i < A.size(); i += 2) {
            match(A[i], A[i+1]);
            if (i+2 < A.size()) addEdge(A[i+1], A[i+2]);
          }
          break;
        }
      }
    }

    vector<pair<int, int>> M;
    for (int u = 0; u < n; u++) {
      if (mate[u] > u) M.emplace_back(u, mate[u]);
    }
    return M;
  }
};

/**
 * BlossomAlgorithmStrategy uses an algorithm based on Edmonds-Blossom algorithm
 * to solve Maximum Matching in O(n^3).
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
