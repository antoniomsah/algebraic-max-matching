#pragma once

#include <vector>
#include <array>

#include "../../config.hpp"
#include "../algorithm-strategy-interface.hpp"

using namespace std;

/**
 * @brief Strategy that uses the algorithm proposed by Harvey. 
 * Solves perfect matching in O(n^\omega). 
 */
class HarveyAlgorithmStrategy : public IAlgorithmStrategy {
  array<vector<int>, 2> DivideInTwo(const vector<int>& V) const {
      const int n = V.size();

      array<vector<int>, 2> S;
      S[0].reserve(n / 2);
      S[1].reserve(n - n / 2);
      for (int i = 0; i < n; i++) {
        S[i >= n / 2].push_back(V[i]);
      }
      return S;
  }

  template <int MOD>
  void DeleteEdgesWithin(const vector<int> &S, TutteMatrix<MOD> &T, TutteMatrix<MOD> &N) const {
      if (S.size() == 1) return;

      array<vector<int>, 2> S_ = DivideInTwo(S);
      for (int i = 0; i < 2; i++) {
        // Save the states.
        TutteMatrix<MOD> TSi = T(S_[i], S_[i]);
        TutteMatrix<MOD> oldNSS = N(S, S);

        DeleteEdgesWithin(S_[i], T, N);

        // Update N[S,S]
        TutteMatrix<MOD> Delta = T(S_[i], S_[i]) - TSi,
          ISi(S_[i].size());

        // Identity matrix
        for (int j = 0; j < S_[i].size(); j++) {
          ISi(j, j) = 1;
        }

        // Array mapping elements to their indices in S,
        // e.g. S = {1, 2, 3, 4} and S_[i] = {2, 4}, then Si = {1, 3}.
        vector<int> Si(S_[i].size());
        for (int j = 0, k = 0; j < S.size() and k < S_[i].size(); j++) {
          if (S[j] == S_[i][k]) Si[k++] = j;
        }
    
        // Updated N[S,S]
        TutteMatrix<MOD> newNSS = oldNSS - oldNSS('*', Si) * (ISi + Delta * oldNSS(Si, Si)).inverse() * Delta * oldNSS(Si, '*');

        // Update N
        for (int j = 0; j < S.size(); j++) {
          for (int k = 0; k < S.size(); k++) {
            N(S[j], S[k]) = newNSS(j, k);
          }
        }
      }
      DeleteEdgesCrossing(S_[0], S_[1], T, N);
  }

  template <int MOD>
  void DeleteEdgesCrossing(
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

    // RS := R \cup S.
    vector<int> RS(R.size() + S.size());
    for (int i = 0; i < RS.size(); i++) {
      if (i < R.size()) {
        RS[i] = R[i];
      } else {
        RS[i] = S[i - R.size()];
      }
    }

    array<vector<int>, 2> RM = DivideInTwo(R),
                          SM = DivideInTwo(S);

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        vector<int> RMSM(RM[i].size() + SM[j].size());
        for (int k = 0; k < RMSM.size(); k++) {
          if (k < RM[i].size()) {
            RMSM[k] = RM[i][k];
          } else {
            RMSM[k] = SM[j][k - RM[i].size()];
          }
        }

        TutteMatrix<MOD> TRMSM = T(RMSM, RMSM);
        TutteMatrix<MOD> oldNRS = N(RS, RS);

        DeleteEdgesCrossing(RM[i], SM[j], T, N);

        TutteMatrix<MOD> Delta = T(RMSM, RMSM) - TRMSM,
                          I(RMSM.size());

        // Identity matrix
        for (int k = 0; k < RMSM.size(); k++) {
          I(k, k) = 1;
        }

        // Array mapping elements to their indices in RS,
        // e.g. RS = {2, 3, 4} and RMSM = {3, 4} then RMSMi = {1, 2}.
        vector<int> RMSMi(RMSM.size());
        for (int j = 0, k = 0; j < RS.size() and k < RMSM.size(); j++) {
          if (RS[j] == RMSM[k]) RMSMi[k++] = j;
        }

        // Update N[R \cup S, R \cup S]
        TutteMatrix<MOD> newNRS = oldNRS - oldNRS('*', RMSMi) * (I + Delta * oldNRS(RMSMi, RMSMi)).inverse() * Delta * oldNRS(RMSMi, '*');

        // Update N
        for (int k = 0; k < RS.size(); k++) {
          for (int l = 0; l < RS.size(); l++) {
            N(RS[k], RS[l]) = newNRS(k, l);
          }
        }
      }
    }
  }

 public:
  /**
   * Finds a perfect matching.
   * Complexity: O(n^{\omega}).
   * @return A perfect matching, if one exists. Else, returns an empty set.
   */
  vector<pair<int, int>> solve(const Graph& G) const override {
    const int n = V(G).size(), m = E(G).size();

    TutteMatrix T = GetTutteMatrix(G);

    if (T.isSingular()) {
      return {};
    }

    TutteMatrix<MOD> N = T.getInverse();
    DeleteEdgesWithin(V(G), T, N);

    vector<pair<int, int>> matching;
    for (int u = 0; u < n; u++) {
      for (int v = u+1; v < n; v++) {
        if (T(u, v).x != 0) matching.emplace_back(u, v);
      }
    }
    return matching;
  }
};