#include <bits/stdc++.h>
#include <filesystem>

#include "../benchmarker.hpp"

using namespace std;

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
  const int NUM_IT = 50;

  Benchmarker benchmarker("maximum_matching", NUM_IT);
  try {
    fs::path inputDir("input");
    if (!fs::exists(inputDir)) {
      cerr << "Input directory not found" << endl;
      return 1;
    }

    bool filesFound = false;
    for (const auto& entry: fs::directory_iterator(inputDir)) {
      if (!entry.is_regular_file()) continue;
      filesFound = true;
      cout << "Running benchmark for file: " << entry.path().filename() << endl;
      for (int alg_id = 1; alg_id < NUM_STRATEGIES; alg_id += 1) {
        benchmarker.run(entry.path(), true, alg_id);
      }
    }

    if (!filesFound) {
      cerr << "Input directory is empty." << endl;
      return 1;
    }
    benchmarker.write();
  } catch (const exception& e) {
    cerr << "Error: " << e.what() << endl;
    return 1;
  }
  return 0;
}
