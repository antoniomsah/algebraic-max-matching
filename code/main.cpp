#include <iostream>
#include <vector>

#include "matrix.hpp"
#include "tutte-matrix.hpp"

using namespace std;

// Implementation of the algorithm for perfect matching.

int main() {
  Matrix<int> M(1, 1);
  M.inverse();
  M.determinant();

  vector<pair<int, int>> edges = {{0, 1}, {1, 0}};
  TutteMatrix T(10, edges);
  return 0;
}
