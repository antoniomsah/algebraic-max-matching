#include <iostream>
#include <map>
#include <vector>

#include "../../algorithms/solver.hpp"

using namespace std;

int main(int argc, char* argv[]) {
  MatchingSolver::initialize();

  int n, m;
  cin >> n >> m;

  Graph G(n);
  for (int i = 0; i < m; i++) {
    int u, v;
    cin >> u >> v;
    G.addEdge(u, v);
  }
  assert(argc >= 2);

  MatchingSolver solver(G);
  solver.setStrategy(stoi(argv[1]));
  auto M = solver.MaximumMatching();
  if (!G.hasMatching(M)) {
    cout << -1 << '\n';
    return 0;
  }
  cout << M.size() << '\n';

  return 0;
}