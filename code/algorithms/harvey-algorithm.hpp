#pragma once

#include <vector>
#include <array>

#include "algorithm-strategy-interface.hpp"
#include "../config.hpp"
#include "../utils.hpp"

using namespace std;

/**
 * Implements HarveyAlgorithm. 
 * Time complexity: O(O(multiply)), where O(multiply) is the time complexity of the matrix multiplication algorithm used.
 */
class HarveyAlgorithmStrategy : public IAlgorithmStrategy {
  template <int MOD>
  void DeleteEdgesWithin(const vector<int> &S, TutteMatrix<MOD> &T, TutteMatrix<MOD> &N) const {
    if (S.size() == 1) return;

    auto S_ = DivideInTwo(S);
    for (int i = 0; i < 2; i++) {
      // Save the states.
      TutteMatrix<MOD> oldTSi = T(S_[i], S_[i]);
      TutteMatrix<MOD> oldNSS = N(S, S);

      DeleteEdgesWithin(S_[i], T, N);

      // Array mapping elements to their indices in S,
      // e.g. S = {1, 2, 3, 4} and S_[i] = {2, 4}, then Si = {1, 3}.
      vector<int> Si(S_[i].size());
      for (int j = 0, k = 0; j < S.size() and k < S_[i].size(); j++) {
        if (S[j] == S_[i][k]) Si[k++] = j;
      }

      auto Delta = T(S_[i], S_[i]) - oldTSi;
      auto ISi = Identity<MOD>(S_[i].size());
    
      // Update N[S,S].
      TutteMatrix<MOD> newNSS = oldNSS - oldNSS('*', Si) * (ISi + Delta * oldNSS(Si, Si)).inverse() * Delta * oldNSS(Si, '*');

      // Update N.
      for (int j = 0; j < S.size(); j++)
        for (int k = 0; k < S.size(); k++)
          N(S[j], S[k]) = newNSS(j, k);
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
    if (min(R.size(), S.size()) == 0) return;

    // Need to check a SINGLE edge for the update to be valid
    // Not sufficient to be |R| = 1
    if (max(R.size(), S.size()) == 1) {
      int r = R[0], s = S[0];
      if (T(r, s) != 0 and N(r, s) != -1 / T(r, s)) {
        N(r, s) = N(r, s) * (1 - T(r, s) * N(r, s)) / (T(r, s) * N(r, s) + 1);
        N(s, r) = -N(r, s);
        T.removeEdge(r, s);
      }
      return;
    } 

    vector<int> RS = Unite(R, S);
    auto RM = DivideInTwo(R);
    auto SM = DivideInTwo(S);
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        vector<int> RMSM = Unite(RM[i], SM[j]);
        TutteMatrix oldT = T(RMSM, RMSM);
        TutteMatrix oldNRS = N(RS, RS);

        DeleteEdgesCrossing(RM[i], SM[j], T, N);

        auto Delta = T(RMSM, RMSM) - oldT;
        auto I = Identity<MOD>(RMSM.size());

        // Array mapping elements to their indices in RS,
        // e.g. RS = {2, 3, 4} and RMSM = {3, 4} then RMSMi = {1, 2}.
        vector<int> RMSMi(RMSM.size());
        for (int j = 0, k = 0; j < RS.size() and k < RMSM.size(); j++) {
          if (RS[j] == RMSM[k]) RMSMi[k++] = j;
        }

        // Update N[R \cup S, R \cup S]
        TutteMatrix newNRS = oldNRS - oldNRS('*', RMSMi) * (I + Delta * oldNRS(RMSMi, RMSMi)).inverse() * Delta * oldNRS(RMSMi, '*');

        // Update N
        for (int k = 0; k < RS.size(); k++)
          for (int l = 0; l < RS.size(); l++)
            N(RS[k], RS[l]) = newNRS(k, l);
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
    TutteMatrix T = GetTutteMatrix(G);
    if (T.isSingular()) {
      return {};
    }
    TutteMatrix N = T.inverse();
    DeleteEdgesWithin(V(G), T, N);
    return E(T);
  }
};
