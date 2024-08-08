#pragma once

#include <vector>

#include "../../config.hpp"
#include "../algorithm-strategy-interface.hpp"

using namespace std;

class HarveyAlgorithmStrategy : public IAlgorithmStrategy {
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

    std::function<array<vector<int>, 2>(const vector<int>& V)> 
    divide_in_half = [](const vector<int>& V) {
      const size_t n = V.size();

      array<vector<int>, 2> S;
      S[0].reserve(n / 2);
      S[1].reserve(n - n / 2);
      for (size_t i = 0; i < n; i++) {
        S[i < n / 2].push_back(V[i]);
      }
      return S;
    };

    std::function<void(const vector<int> &R, const vector<int>&S)> 
    delete_edges_crossing = [&](const vector<int> &R, const vector<int> &S) {
      if (R.size() == 1) {
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
      } else {
        array<vector<int>, 2> RM = divide_in_half(R),
                              SM = divide_in_half(S);
        
        vector<int> RS(R.size() + S.size());
        for (size_t i = 0; i < RS.size(); i++) {
          if (i < R.size()) {
            RS[i] = R[i];
          } else {
            RS[i] = S[i - R.size()];
          }
        }

        for (size_t i = 0; i < 2; i++) {
          for (size_t j = 0; j < 2; j++) {
            // save state
            vector<int> RMSM(RM[i].size() + SM[j].size());
            for (size_t k = 0; k < RMSM.size(); k++) {
              if (k < RM.size()) {
                RMSM[k] = RM[i][k];
              } else {
                RMSM[k] = SM[j][k - RM[i].size()];
              }
            }

            TutteMatrix<MOD> TRMSM = T(RMSM, RMSM);

            delete_edges_crossing(RM[i], SM[j]);

            TutteMatrix<MOD> Delta = T(RS, RS) - TRMSM,
                             I(RMSM.size(), RMSM.size());

            for (size_t k = 0; k < RMSM.size(); k++) {
              I(k, k) = 1;
            }
            
            // update N[R \cup S, R \cup S]
            TutteMatrix<MOD> NRS = NRS - N(RS, RMSM) * (I + Delta * N(RMSM, RMSM)).inverse() * Delta * N(RMSM, RS);

            // update N
            for (size_t k = 0; k < RS.size(); k++) {
              for (size_t l = 0; l < RS.size(); l++) {
                N(RS[k], RS[l]) = NRS(k, l);
              }
            }
          }
        }
      }
    };

    std::function<void(const vector<int>&V)> 
    delete_edges_within = [&](const vector<int> &V) {
      if (V.size() == 1) return;

      array<vector<int>, 2> S = divide_in_half(V);
      for (size_t i = 0; i < 2; i++) {
        // save state
        TutteMatrix<MOD> TS = T(S[i], S[i]);

        delete_edges_within(S[i]);

        // update N[S, S]
        TutteMatrix<MOD> Delta = T(S[i], S[i]) - TS,
                        I(S[i].size(), S[i].size());

        // Identity matrix, maybe add to matrix constructor?
        for (size_t j = 0; j < S[i].size(); j++) {
          I(j, j) = 1;
        }
    
        // updated N[S, S]
        TutteMatrix<MOD> NVV = N(V, V) - N(V, S[i]) * (I + Delta * (TS)).inverse() * Delta * N(S[i], V);

        // update N
        for (size_t j = 0; j < V.size(); j++) {
          for (size_t k = 0; k < V.size(); k++) {
            N(V[j], V[k]) = NVV(j, k);
          }
        }
      }

      delete_edges_crossing(S[0], S[1]);
    };

    vector<int> V(n);
    std::iota(V.begin(), V.end(), 0);
    delete_edges_within(V);

    vector<pair<int, int>> matching;
    for (size_t u = 0; u < n; u++) {
      for (size_t v = u+1; v < n; v++) {
        if (T(u, v).x != 0) matching.emplace_back(u, v);
      }
    }
    return matching;
  }
};