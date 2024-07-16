#include <bits/stdc++.h>

using namespace std;

mt19937 rng((int)chrono::steady_clock::now()
                .time_since_epoch()
                .count());  // random number generator

int uniform(int l, int r) {
  uniform_int_distribution<int> uid(l, r);
  return uid(rng);
}

/**
 * @brief Generates a matrix of size n x n with random entries in range [0,
 * argv[2]].
 */
int main(int argc, char* argv[]) {
  assert(argc == 3);
  int n = stoi(argv[1]);
  cout << n << '\n';
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      cout << uniform(0, stoi(argv[2])) << ' ';
    }
    cout << '\n';
  }
  return 0;
}