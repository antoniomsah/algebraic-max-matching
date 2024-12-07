#include <bits/stdc++.h>

#include "../../../classes/graph.hpp"

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
  assert(argc == 4);

  int n = stoi(argv[1]);
  double prob = stod(argv[2]);

  int matchingNumber = stoi(argv[3]);
  assert(2*matchingNumber <= n);

  // Every vertex with index smaller than argv[3]
  // is part of S.
  Graph G(n);

  auto createComponent = [&](int firstVertex, size_t size) {
    for (size_t i = 0; i < size; i++) {
      for (size_t j = i+1; j < size; j++) {
        if (random_with_probability(prob)) {
          assert(G.addEdge(firstVertex+i, firstVertex+j)); // Asserts this edge was created.
        }
      }
    }
  };

  int currVertex = 0, 
      numOdd = n - 2*matchingNumber; // I want at least this number of odd components.

  // Will add all with same size, it is simpler :)
  int size = 1;
  while (numOdd > 0 && numOdd * (size+2) <= n) size += 2;

  for (int i = 0; i < numOdd; i++) {
    createComponent(currVertex, size);
    currVertex += size;
  }
  createComponent(currVertex, n-currVertex);

  cout << n << ' ' << E(G).size() << '\n';
  for (auto [u, v] : E(G)) cout << u << ' ' << v << '\n';
  return 0;
}