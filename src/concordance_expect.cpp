#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// Expected concordance under the fixed-marginal (hypergeometric) null used by
// QuartetConcordance(normalize = TRUE).  For one character each token is
// reassigned at random across the scored leaves, holding the state counts and
// the split sizes fixed.  For a state-pair (i, j) with counts nI, nJ over t
// scored leaves, of which M are on side A of a split, the 2x2 cell vector
// (p, q, r, s) is trivariate-hypergeometric:
//   P(p, r) = C(nI,p) C(nJ,r) C(t-nI-nJ, M-p-r) / C(t, M)
// with q = nI - p, s = nJ - r, and mA = p + r on side A.  These kernels sum the
// per-pair contributions exactly (double sum over p, r), mirroring the R
// reference `.ExpectedQuartet`/`.QuartetExpect` and `.ExpectedTrit`/
// `.CharTritExpect` in R/Concordance.R.  The p-outer / r-inner accumulation
// order matches the R code so the results agree to machine precision.

// C(z, 2) with the same z < 2 -> 0 guard as the R `choose2`.
static inline double choose2(double z) {
  return z < 2.0 ? 0.0 : z * (z - 1.0) / 2.0;
}

// (z - 1) floored at 0, matching the R `pos(z - 1)`.
static inline double posm1(double z) {
  return z < 1.0 ? 0.0 : z - 1.0;
}

// Per-character scored state counts + per-split side-A size, shared by both
// kernels.  `stateOf[t]` is the compact 0-based state index of a scored leaf, or
// -1 if the leaf is NA for this character.
struct CharCounts {
  int tc;                       // number of scored leaves
  std::vector<int> cnt;         // count of each compact state
  std::vector<int> mSideA;      // per split: scored leaves on side A
};

static CharCounts char_counts(const LogicalMatrix& splits,
                              const IntegerMatrix& characters, int c) {
  const int n_taxa = splits.nrow();
  const int n_splits = splits.ncol();
  CharCounts out;
  out.tc = 0;
  out.mSideA.assign(n_splits, 0);

  // Map raw state values to compact 0-based indices (sorted ascending, as R's
  // sort(unique(...)) does; the pair loop is symmetric so order only affects
  // which is "i" vs "j", not the sums).
  std::vector<int> raw;
  raw.reserve(n_taxa);
  std::vector<int> stateOf(n_taxa, -1);
  for (int t = 0; t < n_taxa; ++t) {
    int st = characters(t, c);
    if (IntegerVector::is_na(st)) continue;
    ++out.tc;
    raw.push_back(st);
  }
  std::sort(raw.begin(), raw.end());
  raw.erase(std::unique(raw.begin(), raw.end()), raw.end());
  out.cnt.assign(raw.size(), 0);

  for (int t = 0; t < n_taxa; ++t) {
    int st = characters(t, c);
    if (IntegerVector::is_na(st)) continue;
    int idx = int(std::lower_bound(raw.begin(), raw.end(), st) - raw.begin());
    stateOf[t] = idx;
    ++out.cnt[idx];
    for (int s = 0; s < n_splits; ++s) {
      if (splits(t, s)) ++out.mSideA[s];
    }
  }
  return out;
}

// Exact E[conc], E[dec] for one state-pair over all needed side-A sizes.
// Returns, for each M in 0..tc, the pair's (Econc, Edec); callers pick M[split].
static void expected_quartet_by_M(int nI, int nJ, int tc,
                                  std::vector<double>& eConc,
                                  std::vector<double>& eDec) {
  const int nOther = tc - nI - nJ;
  eConc.assign(tc + 1, 0.0);
  eDec.assign(tc + 1, 0.0);
  // Hoist the log-choose terms out of the M / p / r loops (the dominant cost).
  std::vector<double> lcI(nI + 1), lcJ(nJ + 1), lcO(nOther + 1), lcT(tc + 1);
  for (int p = 0; p <= nI; ++p) lcI[p] = R::lchoose(nI, p);
  for (int r = 0; r <= nJ; ++r) lcJ[r] = R::lchoose(nJ, r);
  for (int o = 0; o <= nOther; ++o) lcO[o] = R::lchoose(nOther, o);
  for (int M = 0; M <= tc; ++M) lcT[M] = R::lchoose(tc, M);
  for (int M = 0; M <= tc; ++M) {
    const double lt = lcT[M];
    double aConc = 0.0, aDec = 0.0;
    const int pMax = std::min(nI, M);
    for (int p = 0; p <= pMax; ++p) {
      const int q = nI - p;
      const double c_p = choose2(p), c_q = choose2(q), lp = lcI[p];
      const int rMax = std::min(nJ, M - p);
      for (int r = 0; r <= rMax; ++r) {
        const int oA = M - p - r;
        if (oA < 0 || oA > nOther) continue;
        const double prob = std::exp(lp + lcJ[r] + lcO[oA] - lt);
        if (prob <= 0.0) continue;
        const int s = nJ - r;
        const double conc = c_p * choose2(s) + c_q * choose2(r);
        aConc += prob * conc;
        aDec += prob * (conc + double(p) * q * r * s);
      }
    }
    eConc[M] = aConc;
    eDec[M] = aDec;
  }
}

// [[Rcpp::export]]
List quartet_expect(const LogicalMatrix splits, const IntegerMatrix characters) {
  const int n_splits = splits.ncol();
  const int n_chars = characters.ncol();
  NumericMatrix eConc(n_splits, n_chars), eDec(n_splits, n_chars);

  for (int c = 0; c < n_chars; ++c) {
    const CharCounts cc = char_counts(splits, characters, c);
    const int nStates = int(cc.cnt.size());
    if (nStates < 2) continue;
    std::vector<double> ec, ed;
    for (int a = 0; a < nStates - 1; ++a) {
      for (int b = a + 1; b < nStates; ++b) {
        expected_quartet_by_M(cc.cnt[a], cc.cnt[b], cc.tc, ec, ed);
        for (int s = 0; s < n_splits; ++s) {
          const int M = cc.mSideA[s];
          eConc(s, c) += ec[M];
          eDec(s, c) += ed[M];
        }
      }
    }
  }
  return List::create(_["concordant"] = eConc, _["decisive"] = eDec);
}

// Exact E[m], E[m*A/wk], E[m*A/wc] for one trit state-pair over all side-A
// sizes.  wc = (nI-1)+ (nJ-1)+ is fixed; wk = (mA-1)+ (tP-mA-1)+ and
// m = min(wc, wk) vary with mA = p + r.  Mirrors R `.ExpectedTrit`.
static void expected_trit_by_M(int nI, int nJ, int tc,
                               std::vector<double>& eM,
                               std::vector<double>& eMAwk,
                               std::vector<double>& eMAwc) {
  const int nOther = tc - nI - nJ;
  const int tP = nI + nJ;
  const double wc = posm1(nI) * posm1(nJ);
  eM.assign(tc + 1, 0.0);
  eMAwk.assign(tc + 1, 0.0);
  eMAwc.assign(tc + 1, 0.0);
  // Hoist the log-choose terms out of the M / p / r loops (the dominant cost).
  std::vector<double> lcI(nI + 1), lcJ(nJ + 1), lcO(nOther + 1), lcT(tc + 1);
  for (int p = 0; p <= nI; ++p) lcI[p] = R::lchoose(nI, p);
  for (int r = 0; r <= nJ; ++r) lcJ[r] = R::lchoose(nJ, r);
  for (int o = 0; o <= nOther; ++o) lcO[o] = R::lchoose(nOther, o);
  for (int M = 0; M <= tc; ++M) lcT[M] = R::lchoose(tc, M);
  for (int M = 0; M <= tc; ++M) {
    const double lt = lcT[M];
    double aM = 0.0, aWk = 0.0, aWc = 0.0;
    const int pMax = std::min(nI, M);
    for (int p = 0; p <= pMax; ++p) {
      const int q = nI - p;
      const double lp = lcI[p];
      const int rMax = std::min(nJ, M - p);
      for (int r = 0; r <= rMax; ++r) {
        const int oA = M - p - r;
        if (oA < 0 || oA > nOther) continue;
        const double prob = std::exp(lp + lcJ[r] + lcO[oA] - lt);
        if (prob <= 0.0) continue;
        const int s = nJ - r;
        const int mA = p + r;
        const double A = posm1(p) * posm1(s) + posm1(q) * posm1(r);
        const double wk = posm1(mA) * posm1(tP - mA);
        const double m = std::min(wc, wk);
        aM += prob * m;
        if (wk > 0.0) aWk += prob * m * A / wk;
        if (wc > 0.0) aWc += prob * m * A / wc;
      }
    }
    eM[M] = aM;
    eMAwk[M] = aWk;
    eMAwc[M] = aWc;
  }
}

// [[Rcpp::export]]
List trit_expect(const LogicalMatrix splits, const IntegerMatrix characters) {
  const int n_splits = splits.ncol();
  const int n_chars = characters.ncol();
  NumericMatrix numEdge(n_splits, n_chars);
  NumericMatrix numChar(n_splits, n_chars);
  NumericMatrix denM(n_splits, n_chars);

  for (int c = 0; c < n_chars; ++c) {
    const CharCounts cc = char_counts(splits, characters, c);
    const int nStates = int(cc.cnt.size());
    if (nStates < 2) continue;
    std::vector<double> eM, eWk, eWc;
    for (int a = 0; a < nStates - 1; ++a) {
      for (int b = a + 1; b < nStates; ++b) {
        expected_trit_by_M(cc.cnt[a], cc.cnt[b], cc.tc, eM, eWk, eWc);
        for (int s = 0; s < n_splits; ++s) {
          const int M = cc.mSideA[s];
          denM(s, c) += eM[M];
          numEdge(s, c) += eWk[M];
          numChar(s, c) += eWc[M];
        }
      }
    }
  }
  return List::create(_["numEdge"] = numEdge, _["numChar"] = numChar,
                      _["denM"] = denM);
}
