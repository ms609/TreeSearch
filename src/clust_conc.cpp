#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <numeric>

using namespace Rcpp;

// ---------- small utilities ----------

static inline int pop8(unsigned x) {
  // popcount for a single byte (0..255)
  // Using builtin would be fine, but a tiny LUT avoids any compiler/arch quirks.
  static bool init = false;
  static unsigned char lut[256];
  if (!init) {
    for (int i = 0; i < 256; ++i) {
      unsigned v = (unsigned)i, c = 0;
      while (v) { v &= v - 1; ++c; }
      lut[i] = (unsigned char)c;
    }
    init = true;
  }
  return lut[x & 0xFFu];
}

static inline int popcount_rawrow(const RawMatrix& M, int row, int nBytes) {
  int s = 0;
  for (int b = 0; b < nBytes; ++b) s += pop8((unsigned)M(row, b));
  return s;
}

static inline int popcount_and(const RawMatrix& A, int arow,
                               const std::vector<unsigned char>& B, int nBytes) {
  int s = 0;
  for (int b = 0; b < nBytes; ++b) s += pop8(((unsigned)A(arow, b)) & ((unsigned)B[b]));
  return s;
}

static inline int popcount_and2(const std::vector<unsigned char>& A,
                                const std::vector<unsigned char>& B, int nBytes) {
  int s = 0;
  for (int b = 0; b < nBytes; ++b) s += pop8(((unsigned)A[b]) & ((unsigned)B[b]));
  return s;
}

static inline void bitset_or(std::vector<unsigned char>& dst,
                             const std::vector<unsigned char>& src, int nBytes) {
  for (int b = 0; b < nBytes; ++b) dst[b] |= src[b];
}

static inline void bitset_set(std::vector<unsigned char>& bs, int tip) {
  bs[tip >> 3] |= (1u << (tip & 7));
}

static inline double entropy_from_counts(const std::vector<int>& cnt) {
  // Shannon entropy in bits
  int n = 0;
  for (int c : cnt) n += c;
  if (n <= 0) return 0.0;
  double invn = 1.0 / (double)n;
  double H = 0.0;
  for (int c : cnt) {
    if (c > 0) {
      double p = c * invn;
      H -= p * std::log2(p);
    }
  }
  return H;
}

static inline double entropy_2way(int a, int b) {
  if (a < 0 || b < 0) return 0.0;
  int n = a + b;
  if (n <= 0) return 0.0;
  if (a == 0 || b == 0) return 0.0;
  double pa = (double)a / (double)n;
  double pb = 1.0 - pa;
  return - (pa * std::log2(pa) + pb * std::log2(pb));
}

// ---------- main exported function ----------

// [[Rcpp::export]]
double cc_tree_normalized_cpp(const IntegerMatrix& mat,   // nTips x nChar
                              const RawMatrix& splitsRaw, // nSplits x nBytes (bitmasks)
                              int nSamples) {
  const int nTips   = mat.nrow();
  const int nChar   = mat.ncol();
  const int nSplits = splitsRaw.nrow();
  const int nBytes  = splitsRaw.ncol();
  
  if (nTips <= 0 || nChar <= 0 || nSplits <= 0) return NA_REAL;
  
  // ---- Precompute per-character structures ----
  struct CharInfo {
    int nActive = 0;
    int chMax = 0;
    double Hx = 0.0;
    std::vector<int> stateCount;                    // size = chMax
    std::vector< std::vector<unsigned char> > S;    // S[k] bitmask where state == (k+1)
    std::vector<unsigned char> activeMask;          // bitmask of non-NA tips
  };
  std::vector<CharInfo> C(nChar);
  
  for (int j = 0; j < nChar; ++j) {
    // find chMax and active
    int chMax = 0, nAct = 0;
    for (int i = 0; i < nTips; ++i) {
      int v = mat(i, j);
      if (v != NA_INTEGER) {
        if (v > chMax) chMax = v;
        ++nAct;
      }
    }
    if (nAct == 0) {
      C[j].nActive = 0;
      C[j].chMax = 0;
      C[j].Hx = 0.0;
      C[j].activeMask.assign(nBytes, 0);
      continue;
    }
    
    C[j].nActive = nAct;
    C[j].chMax = chMax;
    C[j].stateCount.assign(chMax, 0);
    C[j].S.assign(chMax, std::vector<unsigned char>(nBytes, 0));
    C[j].activeMask.assign(nBytes, 0);
    
    for (int i = 0; i < nTips; ++i) {
      int v = mat(i, j);
      if (v != NA_INTEGER) {
        bitset_set(C[j].activeMask, i);
        if (v >= 1 && v <= chMax) {
          ++C[j].stateCount[v - 1];
          bitset_set(C[j].S[v - 1], i);
        }
      }
    }
    C[j].Hx = entropy_from_counts(C[j].stateCount);
  }
  
  // ---- Observed weighted sums ----
  double sumW = 0.0;        // sum of hBest_ij
  double sumW2 = 0.0;       // sum of hBest_ij^2  (for Hmax)
  double sumW_MI = 0.0;     // sum of hBest_ij * MI_ij
  
  for (int s = 0; s < nSplits; ++s) {
    for (int j = 0; j < nChar; ++j) {
      if (C[j].nActive == 0) continue;
      
      // counts on split among active
      int n1 = popcount_and(splitsRaw, s, C[j].activeMask, nBytes); // side "1"
      int n0 = C[j].nActive - n1;
      
      if (n0 < 2 || n1 < 2) {
        // trivial by your rule -> MI = 0, hBest = 0 => contributes nothing
        continue;
      }
      
      double HB = entropy_2way(n0, n1);
      double HX = C[j].Hx;
      double hBest = std::min(HX, HB);
      
      // joint counts: [B=0, states...] then [B=1, states...]
      std::vector<int> joint(2 * C[j].chMax, 0);
      for (int k = 0; k < C[j].chMax; ++k) {
        int c1 = popcount_and2(C[j].S[k], C[j].activeMask, nBytes); // sanity: equals stateCount[k]
        // but we need side=1 among active on the split:
        c1 = popcount_and(splitsRaw, s, C[j].S[k], nBytes); // (B=1, X=k+1)
        int tot = C[j].stateCount[k];
        int c0 = tot - c1;
        joint[k] = c0;
        joint[C[j].chMax + k] = c1;
      }
      double HBX = entropy_from_counts(joint);
      double MI  = HX + HB - HBX;
      
      sumW     += hBest;
      sumW2    += hBest * hBest;
      sumW_MI  += hBest * MI;
    }
  }
  
  if (sumW == 0.0) return NA_REAL; // no informative pairs
  
  const double M_obs = sumW_MI / sumW;
  const double Hmax  = sumW2    / sumW;
  
  // ---- If no numeric normalization requested, just return M_obs / Hmax (on [0,1]) ----
  if (nSamples <= 0) {
    return (Hmax > 0.0) ? (M_obs / Hmax) : NA_REAL;
  }
  
  // ---- Randomization: tip relabelling via multivariate hypergeometric sampling ----
  // Precompute total "1"s per split (across all tips)
  std::vector<int> splitOne(nSplits, 0);
  for (int s = 0; s < nSplits; ++s) {
    splitOne[s] = popcount_rawrow(splitsRaw, s, nBytes);
  }
  
  // Accumulate expected weighted mean MI under randomization
  double sumWeightedRandMeans = 0.0;
  
  RNGScope scope; // use R's RNG (rhyper) safely
  
  for (int rep = 0; rep < nSamples; ++rep) {
    double acc = 0.0; // sum w_ij * MI_rand_ij for this replicate
    
    for (int s = 0; s < nSplits; ++s) {
      const int L = splitOne[s]; // number of tips on side 1 (fixed)
      const int N = nTips;
      
      for (int j = 0; j < nChar; ++j) {
        if (C[j].nActive == 0) continue;
        
        // Draw how many active tips fall on side 1:
        // rhyper(nn1=L, nn2=N-L, k = nActive) -> successes among draws
        int n1 = (int) ::Rf_rhyper((double)L, (double)(N - L), (double)C[j].nActive);
        int n0 = C[j].nActive - n1;
        
        if (n0 < 2 || n1 < 2) {
          // trivial -> MI_rand = 0 (matches your rule)
          continue;
        }
        
        // Now distribute state counts to side 1 by multivariate hypergeometric
        // Sequential HG draws:
        int drawsLeft = n1;
        int totalLeft = C[j].nActive;
        std::vector<int> c1_state(C[j].chMax, 0);
        
        for (int k = 0; k < C[j].chMax - 1; ++k) {
          int m = C[j].stateCount[k];
          if (m <= 0) { c1_state[k] = 0; continue; }
          // rhyper(nn1 = m, nn2 = totalLeft - m, k = drawsLeft)
          int take = (int) ::Rf_rhyper((double)m, (double)(totalLeft - m), (double)drawsLeft);
          if (take > m) take = m;
          if (take > drawsLeft) take = drawsLeft;
          c1_state[k] = take;
          drawsLeft  -= take;
          totalLeft  -= m;
        }
        c1_state[C[j].chMax - 1] = std::max(0, drawsLeft); // leftover to last state
        
        // Entropies for MI
        const double HB = entropy_2way(n0, n1);
        const double HX = C[j].Hx;
        const double hBest = std::min(HX, HB);
        
        // joint entropy
        std::vector<int> joint(2 * C[j].chMax, 0);
        for (int k = 0; k < C[j].chMax; ++k) {
          int c1 = c1_state[k];
          int tot = C[j].stateCount[k];
          int c0 = tot - c1;
          joint[k] = c0;
          joint[C[j].chMax + k] = c1;
        }
        double HBX = entropy_from_counts(joint);
        double MIr = HX + HB - HBX;
        
        acc += hBest * MIr;
      }
    }
    
    sumWeightedRandMeans += (acc / sumW);
  }
  
  const double M_rand = sumWeightedRandMeans / (double)nSamples;
  
  // Final normalized score on [0,1]: 0 at random expectation; 1 at "perfect"
  const double denom = (Hmax - M_rand);
  if (denom <= 0.0) return NA_REAL;
  return (M_obs - M_rand) / denom;
}
