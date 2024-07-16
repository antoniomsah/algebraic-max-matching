#include <iostream>
#include <map>
#include <vector>

#include "algorithms/algorithm-strategy-interface.hpp"
#include "algorithms/solver.hpp"
#include "matching-verifier.hpp"

using namespace std;

array<std::shared_ptr<IAlgorithmStrategy>, NUM_STRATEGIES>
    MatchingSolver::strategies;

int main(int argc, char* argv[]) {
  MatchingSolver::initialize();

  int n, m;
  cin >> n >> m;

  vector<vector<int>> graph(n, vector<int>(n));
  vector<pair<int, int>> edges(m);
  for (int i = 0; i < m; i++) {
    int u, v;
    cin >> u >> v;
    edges[i] = {u, v};
    graph[u][v] = graph[v][u] = 1;
  }

  assert(argc >= 2);

  MatchingSolver solver;
  solver.set_strategy(stoi(argv[1]));
  solver.solve(graph, edges);

  // if (TEST) {
  //   /* For testing purposes. */
  //   int correct_outputs = 0;
  //   for (auto [matching, ocr] : found) {
  //     if (MatchingValidator::validate(graph, matching, true))
  //       correct_outputs += ocr;
  //     cout << '\n';
  //     for (auto [u, v] : matching) cout << u << ' ' << v << '\n';
  //   }

  //   cout << "found " << found.size() << " different outputs and generated "
  //        << correct_outputs << " correct outputs out of " << MAX_IT
  //        << " iterations (prob = " << double(correct_outputs) / MAX_IT <<
  //        ")\n";
  // } else {
  //   /* Prints a matching. */
  //   vector<pair<int, int>> res;
  //   for (auto [matching, ocr] : found) {
  //     if (MatchingValidator::validate(graph, matching, true)) {
  //       res = matching;
  //       break;
  //     }
  //   }

  //   for (auto [u, v] : res) {
  //     cout << u << ' ' << v << '\n';
  //   }
  // }

  return 0;
}
