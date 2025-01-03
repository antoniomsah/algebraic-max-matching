#include <bits/stdc++.h>

using namespace std;

int random_with_probability(double p) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::bernoulli_distribution d(p);
  return d(gen);
}

/**
 * @brief Generates a random test with n vertices and probability argv[2] of
 * adding an edge between two vertices.
 */
int main(int argc, char* argv[]) {
  assert(argc == 3);
  int n = stoi(argv[1]);
  set<pair<int, int>> edges;
  double prob = stod(argv[2]);
  for (int u = 0; u < n; u++) {
    for (int v = u + 1; v < n; v++) {
      if (random_with_probability(prob)) edges.emplace(u, v);
    }
  }

  cout << n << ' ' << edges.size() << '\n';
  for (auto [u, v] : edges) cout << u << ' ' << v << '\n';

  return 0;
}