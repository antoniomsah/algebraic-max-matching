#pragma once

#include <bits/stdc++.h>
#include <chrono>
#include <clocale>

using namespace std;

#include "../solver.hpp"

struct Result {
  string algorithmName;
  double averageTime;
};

class Benchmarker {
private:
  const int NUM_IT;
  map<string, vector<Result>> results;
  ofstream outFile;

  Graph readInput(string filename) {
    if (!freopen(filename.c_str(), "r", stdin)) {
      cerr << "Failed to open file: " << filename << endl;
      return Graph(0);
    }

    int n, m;
    cin >> n >> m;

    Graph G(n);
    for (int i = 0; i < m; i++) {
      int u, v;
      cin >> u >> v;
      G.addEdge(u, v);
    }
    return G;
  }
public:
  Benchmarker(string name, int num_it) : NUM_IT(num_it) {
    MatchingSolver::initialize();
    outFile.open(name + ".csv");
    outFile << "Input,Algorithm,Average Time (ms)\n";
  }

  void run(string filename, bool maximumMatching, int alg_id) {
    MatchingSolver solver(readInput(filename));
    solver.setStrategy(alg_id);

    double total_time = 0.0;
    for (int it = 0; it < NUM_IT; it++) {
      auto start = chrono::high_resolution_clock::now();

      // Runs the algorithm.
      if (maximumMatching) {
        solver.MaximumMatching();
      } else {
        solver.PerfectMatching();
      }

      auto end=  std::chrono::high_resolution_clock::now();
      auto duration = chrono::duration_cast<chrono::microseconds>(end - start);

      double time = duration.count() / 1000.0; // Time in ms.
      total_time += time;
    }
    results[filename].emplace_back(solver.algorithmName(), total_time / NUM_IT);
  }

  void write() {
    outFile << fixed << setprecision(4);
    for (const auto& [filename, results]: results) {
      for (const auto& result: results) {
        outFile << filename << "," << result.algorithmName << "," << result.averageTime << "\n";
      }
    }
    outFile.close();
  }

  void print() {
    cout << fixed << setprecision(4);
    cout << "Summary\n";
    for (const auto& [filename, results]: results) {
      for (const auto& result: results) {
        cout << "Input: " << filename << "(" << result.algorithmName << "): " << result.averageTime << "ms average\n";
      }
    }
  }
};
