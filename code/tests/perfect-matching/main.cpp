#include <iostream>
#include <map>
#include <vector>

#include "../../algorithms/solver.hpp"

using namespace std;

const string OUTPUT_DIR = "./output";

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

  auto M = solver.PerfectMatching();
  if (M.empty() || !G.hasMatching(M)) {
    cout << "NO\n";
    return 0;
  }
  cout << "YES\n";

  return 0;
}
