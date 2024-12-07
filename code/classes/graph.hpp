#pragma once

#include <vector>

#include "tutte-matrix.hpp"
#include "../config.hpp"

using std::vector;
using std::pair;

// Graph class. Does NOT support addition of vertices.
class Graph {
    vector<int> V;
    vector<pair<int, int>> E;
    Matrix<int> adj; // Adjacency matrix, adj(u, v) = index of edge connecting u and v.
    TutteMatrix<MOD> T;

    bool hasVertex(int u) { return u < V.size(); }

public:
    Graph(int n) : V(n), adj(n, n), T(TutteMatrix<MOD>(n)) {
        for (size_t i = 0; i < n; i++) {
            V[i] = i;
            for (size_t j = 0; j < n; j++) {
                adj(i, j) = -1;
            }
        }
    }

    int size() { return V.size(); }

    bool isAdj(int u, int v) {
        if (!hasVertex(u) || !hasVertex(v)) return false;
        return adj(u, v) != -1;
    }

    // addEdge adds an edge uv to G.
    // Returns false if it fails to add the edge.
    bool addEdge(int u, int v) {
        if (!hasVertex(u) || !hasVertex(v) || isAdj(u, v)) return false;
        adj(u, v) = adj(v, u) = E.size();
        E.emplace_back(u, v);
        T.addEdge(u, v);
        return true;
    }

    void removeEdge(int u, int v) {
        if (!hasVertex(u) || !hasVertex(v) || !isAdj(u, v)) return;
        adj(u, v) = adj(v, u) = 0;
        T.removeEdge(u, v);
    }

    void removeEdge(int i) {
        const auto& [u, v] = E[i];
        removeEdge(u, v);
    }

    friend vector<int> V(const Graph& G) { return G.V; }
    friend vector<pair<int, int>> E(const Graph& G) { return G.E; }

    // GetTutteMatrix returns a underlying Tutte Matrix of G.
    friend TutteMatrix<MOD> GetTutteMatrix(const Graph& G) { return G.T; }
};