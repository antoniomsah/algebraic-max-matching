#include <iostream>
#include <map>
#include <vector>

#include "matching-verifier.hpp"
#include "rabin-vazirani-algorithm.hpp"
#include "simple-algorithm.hpp"

using namespace std;

int main(int argc, char* argv[]) {
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

  map<vector<pair<int, int>>, int> found;
  for (int it = 0; it < MAX_IT; it++) {
    vector<pair<int, int>> matching;
    if (argc == 1) {
      matching = SimpleAlgorithm::solve(graph, edges);
    } else {
      switch (stoi(argv[1])) {
        case 0:
          matching = SimpleAlgorithm::solve(graph, edges);
          break;
        default:
          matching = RabinVaziraniAlgorithm::solve(graph, edges);
      }
    }
    found[matching] += 1;
  }

  if (TEST) {
    /* For testing purposes. */
    int correct_outputs = 0;
    for (auto [matching, ocr] : found) {
      if (MatchingValidator::validate(graph, matching, true))
        correct_outputs += ocr;
    }

    cout << "found " << found.size() << " different outputs and generated "
         << correct_outputs << " correct outputs out of " << MAX_IT
         << " iterations (prob = " << double(correct_outputs) / MAX_IT << ")\n";
  } else {
    /* Prints a matching. */
    vector<pair<int, int>> res;
    for (auto [matching, ocr] : found) {
      if (MatchingValidator::validate(graph, matching, true)) {
        res = matching;
        break;
      }
    }

    for (auto [u, v] : res) {
      cout << u << ' ' << v << '\n';
    }
  }

  return 0;
}
