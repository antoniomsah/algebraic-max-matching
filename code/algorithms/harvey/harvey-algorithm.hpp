#pragma once

#include <vector>
#include <array>

#include "../../config.hpp"
#include "../algorithm-strategy-interface.hpp"

using namespace std;

class HarveyAlgorithmStrategy : public IAlgorithmStrategy {
  array<vector<int>, 2> divide_in_half(const vector<int>& V) const {
      const size_t n = V.size();

      array<vector<int>, 2> S;
      S[0].reserve(n / 2);
      S[1].reserve(n - n / 2);
      for (size_t i = 0; i < n; i++) {
        S[i >= n / 2].push_back(V[i]);
      }
      return S;
  }

  template <int MOD>
  void delete_edges_within(const vector<int> &S, TutteMatrix<MOD> &T, TutteMatrix<MOD> &N) const {
      if (S.size() == 1) return;

      array<vector<int>, 2> S_ = divide_in_half(S);
      for (size_t i = 0; i < 2; i++) {
        // save state
        TutteMatrix<MOD> TSi = T(S_[i], S_[i]);
        TutteMatrix<MOD> oldNSS = N;

        delete_edges_within(S_[i], T, N);

        // update N[S,S]
        TutteMatrix<MOD> Delta = T(S_[i], S_[i]) - TSi,
          ISi(S_[i].size(), S_[i].size());

        // Identity matrix, maybe add to matrix constructor?
        for (size_t j = 0; j < S_[i].size(); j++) {
          ISi(j, j) = 1;
        }
    
        // updated N[S,S]
        TutteMatrix<MOD> newNSS = oldNSS(S, S) - oldNSS(S, S_[i]) * (ISi + Delta * oldNSS(S_[i], S_[i])).inverse() * Delta * oldNSS(S_[i], S);

        // update N
        for (size_t j = 0; j < S.size(); j++) {
          for (size_t k = 0; k < S.size(); k++) {
            N(S[j], S[k]) = newNSS(j, k);
          }
        }
      }
      delete_edges_crossing(S_[0], S_[1], T, N);
  }

  template <int MOD>
  void delete_edges_crossing(
    const vector<int> &R, 
    const vector<int> &S,
    TutteMatrix<MOD> &T,
    TutteMatrix<MOD> &N) const {
    // There are no edges in this case
    if (min(R.size(), S.size()) == 0) {
      return;
    }

    // Need to check a SINGLE edge for the update to be valid
    // Not sufficient to be |R| = 1
    if (max(R.size(), S.size()) == 1) {
      for (const int& r : R) {
        for (const int& s : S) {
          if (T(r, s).x != 0 and N(r, s) != T(r, s).inv() * (-1)) {
            N(r, s) = N(r, s) * (T(r, s) * N(r, s) * (-1) + 1) / (T(r, s) * N(r, s) + 1);
            N(s, r) = N(r, s) * (-1);

            // remove edge
            T(r, s) = T(s, r) = 0;
          }
        }
      }
      return;
    } 

    vector<int> RS(R.size() + S.size());
    for (size_t i = 0; i < RS.size(); i++) {
      if (i < R.size()) {
        RS[i] = R[i];
      } else {
        RS[i] = S[i - R.size()];
      }
    }

    array<vector<int>, 2> RM = divide_in_half(R),
                          SM = divide_in_half(S);

    for (size_t i = 0; i < 2; i++) {
      for (size_t j = 0; j < 2; j++) {
        // save state
        vector<int> RMSM(RM[i].size() + SM[j].size());
        for (size_t k = 0; k < RMSM.size(); k++) {
          if (k < RM[i].size()) {
            RMSM[k] = RM[i][k];
          } else {
            RMSM[k] = SM[j][k - RM[i].size()];
          }
        }

        TutteMatrix<MOD> TRMSM = T(RMSM, RMSM);
        TutteMatrix<MOD> oldNRS = N;

        delete_edges_crossing(RM[i], SM[j], T, N);

        TutteMatrix<MOD> Delta = T(RMSM, RMSM) - TRMSM,
                          I(RMSM.size(), RMSM.size());

        for (size_t k = 0; k < RMSM.size(); k++) {
          I(k, k) = 1;
        }
        
        // update N[R \cup S, R \cup S]
        TutteMatrix<MOD> newNRS = oldNRS(RS, RS) - oldNRS(RS, RMSM) * (I + Delta * oldNRS(RMSM, RMSM)).inverse() * Delta * oldNRS(RMSM, RS);

        // update N
        // This update is FAILING
        for (size_t k = 0; k < RS.size(); k++) {
          for (size_t l = 0; l < RS.size(); l++) {
            N(RS[k], RS[l]) = newNRS(k, l);
          }
        }
      }
    }
  }

 public:
  /**
   * @brief Finds a perfect matching.
   * Complexity: O(n^{\omega}).
   */
  vector<pair<int, int>> solve(
      const vector<vector<int>>& graph,
      const vector<pair<int, int>>& edges) const override {
    const size_t n = graph.size(), m = edges.size();

    TutteMatrix T = build_tutte_matrix<MOD>(n, edges);

    std::cerr << "Running Harvey's algorithm" << std::endl;

    if (T.is_singular()) {
      return {};
    }

    TutteMatrix<MOD> N = T.inverse();

    vector<int> V(n);
    std::iota(V.begin(), V.end(), 0);
    delete_edges_within(V, T, N);

    vector<pair<int, int>> matching;
    for (size_t u = 0; u < n; u++) {
      for (size_t v = u+1; v < n; v++) {
        if (T(u, v).x != 0) matching.emplace_back(u, v);
      }
    }
    return matching;
  }
};