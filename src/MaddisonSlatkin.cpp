// MaddisonSlatkin.cpp
#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

static const double NEG_INF = -std::numeric_limits<double>::infinity();

// ----- Custom hash for std::vector<int> with cached hash
struct vec_int_hash {
  std::size_t operator()(const std::vector<int>& v) const noexcept {
    std::size_t h = 0;
    for(int x : v) {
      h ^= std::hash<int>{}(x) + 0x9e3779b9 + (h << 6) + (h >> 2);
    }
    return h;
  }
};

// ----- Utilities
// Streaming log-sum-exp accumulator: identical semantics to log_sum_exp on a container,
// but avoids allocating vectors when accumulating many terms incrementally.
struct LSEAccumulator {
  long double maxv;
  long double acc; // accumulates sum of exp(v - maxv)
  bool empty;
  
  LSEAccumulator() : maxv(-std::numeric_limits<long double>::infinity()), acc(0.0L), empty(true) {}
  
  inline void add(long double v) {
    // ignore non-finite terms (same semantics as log_sum_exp on a container
    // which effectively ignores -Inf unless all terms are -Inf).
    if (!std::isfinite((double)v)) return;
    
    if (empty) {
      maxv = v;
      acc = 1.0L;        // exp(v - maxv) = exp(0) = 1
      empty = false;
    } else {
      if (v > maxv) {
        // bring old accumulator to the scale of the new max
        acc = acc * expl(maxv - v);
        maxv = v;
      }
      acc += expl(v - maxv);
    }
  }
  
  inline double result() const {
    if (empty) return NEG_INF;
    long double res = maxv + logl(acc);
    return (double)res;
  }
  
  inline bool is_empty() const { return empty; }
};

static inline double log_sum_exp(const std::vector<double>& xs) {
  if (xs.empty()) return NEG_INF;
  double m = -std::numeric_limits<double>::infinity();
  for (double v : xs) if (v > m) m = v;
  if (!R_finite(m)) return m;
  long double acc = 0.0L;
  for (double v : xs) acc += std::exp((long double)(v - m));
  return m + std::log((double)acc);
}

static inline double log_prod_sum(const std::vector<double>& xs) {
  for (double v : xs) if (!R_finite(v)) return NEG_INF;
  long double s = 0.0L;
  for (double v : xs) s += v;
  return (double)s;
}

static inline int sum_int(const std::vector<int>& v) {
  int s = 0;
  for (int x : v) s += x;
  return s;
}

static inline int token_mask(int tokenIdx) { return tokenIdx + 1; }

static inline int lex_compare(const std::vector<int>& a, const std::vector<int>& b) {
  for (size_t i = 0; i < a.size(); ++i) {
    if (a[i] < b[i]) return -1;
    if (a[i] > b[i]) return 1;
  }
  return 0;
}

// ----- LnRooted(n)
struct LnRootedCache {
  std::vector<double> lnR;
  
  explicit LnRootedCache(int n_max) {
    lnR.assign(std::max(2, n_max + 1), 0.0);
    for (int n = 2; n <= n_max; ++n) {
      lnR[n] = lnR[n-1] + std::log((double)(2*n - 3));
    }
  }
  
  inline double operator()(int n) const {
    if (n < 0) return NEG_INF;
    if (n < (int)lnR.size()) return lnR[n];
    double acc = lnR.back();
    for (int i = (int)lnR.size(); i <= n; ++i) {
      acc += std::log((double)(2*i - 3));
    }
    return acc;
  }
};

static inline double lchoose_log(int n, int k) {
  if (k < 0 || k > n) return NEG_INF;
  return R::lgammafn(n + 1.0) - R::lgammafn(k + 1.0) - R::lgammafn(n - k + 1.0);
}

// ----- Downpass
struct Downpass {
  int nStates;
  std::vector<int> dp;
  std::vector<uint8_t> step;
  
  explicit Downpass(int nLevels) {
    nStates = (1 << nLevels) - 1;
    dp.resize(nStates * nStates);
    step.resize(nStates * nStates);
    
    for (int i = 0; i < nStates; ++i) {
      int a = i + 1;
      for (int j = 0; j < nStates; ++j) {
        int b = j + 1;
        int inter = (a & b);
        bool newStep = (inter == 0);
        int outMask = newStep ? (a | b) : inter;
        dp[i * nStates + j] = outMask;
        step[i * nStates + j] = newStep ? 1u : 0u;
      }
    }
  }
  
  inline int dp_at(int i, int j) const { return dp[i * nStates + j]; }
  inline bool step_at(int i, int j) const { return step[i * nStates + j] != 0u; }
};

// ----- Pair structures
struct Pair { int a; int b; };

struct TokenPairs {
  std::vector< std::vector<Pair> > noStep;
  std::vector< std::vector<Pair> > yesStep;
  
  TokenPairs(const Downpass& D, int presentBits) {
    int nStates = D.nStates;
    noStep.assign(nStates, {});
    yesStep.assign(nStates, {});
    for (int i = 0; i < nStates; ++i) {
      for (int j = 0; j < nStates; ++j) {
        int tok = D.dp_at(i, j);
        int outIdx = tok - 1;
        if ((tok & ~presentBits) != 0) continue;
        if (D.step_at(i, j)) yesStep[outIdx].push_back({i, j});
        else                 noStep[outIdx].push_back({i, j});
      }
    }
  }
};

// ----- ValidDraws
struct DrawPair {
  std::vector<int> drawn;
  std::vector<int> undrawn;
};

class ValidDrawsCache {
  std::unordered_map<std::vector<int>, std::vector<DrawPair>, vec_int_hash> cache;
  
  void rec(const std::vector<int>& leaves, std::vector<int>& drawn, int idx,
           int n, int half, int curSum, std::vector<DrawPair>& out) {
    if (idx == (int)leaves.size()) {
      if (curSum == 0 || curSum > half) return;
      if (curSum * 2 == n) {
        std::vector<int> undrawn_local(leaves.size());
        for (size_t i = 0; i < leaves.size(); ++i) undrawn_local[i] = leaves[i] - drawn[i];
        int cmp = lex_compare(drawn, undrawn_local);
        if (cmp > 0) return;
      }
      // compute undrawn once and push both
      DrawPair dp;
      dp.drawn = drawn;
      dp.undrawn.resize(leaves.size());
      for (size_t i = 0; i < leaves.size(); ++i) dp.undrawn[i] = leaves[i] - drawn[i];
      out.push_back(std::move(dp));
      return;
    }
    int maxTake = std::min(leaves[idx], half - curSum);
    for (int k = 0; k <= maxTake; ++k) {
      drawn[idx] = k;
      rec(leaves, drawn, idx + 1, n, half, curSum + k, out);
    }
    drawn[idx] = 0;
  }
  
public:
  // Now returns a vector of DrawPair (drawn + undrawn), precomputed
  const std::vector<DrawPair>& get(const std::vector<int>& leaves) {
    auto it = cache.find(leaves);
    if (it != cache.end()) return it->second;
    int n = sum_int(leaves);
    int half = n / 2;
    std::vector<DrawPair> out;
    std::vector<int> drawn(leaves.size(), 0);
    rec(leaves, drawn, 0, n, half, 0, out);
    auto ins = cache.emplace(leaves, std::move(out));
    return ins.first->second;
  }
};

class GlobalValidDraws {
public:
  static ValidDrawsCache& get() {
    static ValidDrawsCache instance;   // constructed once per R session
    return instance;
  }
};

// ----- LogRD using pair of vectors as key (with cached hash)
struct LogRDKey {
  std::vector<int> drawn;
  std::vector<int> leaves;
  mutable size_t cached_hash;
  mutable bool hash_computed;
  
  LogRDKey() : cached_hash(0), hash_computed(false) {}
  LogRDKey(const std::vector<int>& d, const std::vector<int>& l) 
    : drawn(d), leaves(l), cached_hash(0), hash_computed(false) {}
  
  bool operator==(const LogRDKey& other) const {
    return drawn == other.drawn && leaves == other.leaves;
  }
  
  size_t hash() const {
    if (!hash_computed) {
      cached_hash = 0;
      for (int x : drawn) cached_hash ^= std::hash<int>{}(x) + 0x9e3779b9 + (cached_hash << 6) + (cached_hash >> 2);
      for (int x : leaves) cached_hash ^= std::hash<int>{}(x) + 0x9e3779b9 + (cached_hash << 6) + (cached_hash >> 2);
      hash_computed = true;
    }
    return cached_hash;
  }
};

struct LogRDKeyHash {
  size_t operator()(const LogRDKey& k) const noexcept {
    return k.hash();
  }
};

class LogRDCache {
  LnRootedCache& lnRooted;
  std::unordered_map<LogRDKey, double, LogRDKeyHash> cache;
public:
  explicit LogRDCache(LnRootedCache& lnr) : lnRooted(lnr) {
    cache.reserve(1024);
  }
  
  double compute(const std::vector<int>& drawn, const std::vector<int>& leaves) {
    LogRDKey key(drawn, leaves);
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;
    
    const int m = sum_int(drawn);
    const int n = sum_int(leaves);
    double bal = (n == 2*m) ? std::log(0.5) : 0.0;
    long double lc = 0.0L;
    for (size_t i = 0; i < leaves.size(); ++i) lc += lchoose_log(leaves[i], drawn[i]);
    double val = bal + lnRooted(m) + lnRooted(n - m) - lnRooted(n) + (double)lc;
    cache.emplace(std::move(key), val);
    return val;
  }
};

// ----- LogPKey (with cached hash)
struct LogPKey {
  int s;
  int token;
  std::vector<int> leaves;
  mutable size_t cached_hash;
  mutable bool hash_computed;
  
  LogPKey() : s(0), token(0), cached_hash(0), hash_computed(false) {}
  LogPKey(int s_, int t_, const std::vector<int>& l_) 
    : s(s_), token(t_), leaves(l_), cached_hash(0), hash_computed(false) {}
  
  bool operator==(const LogPKey& other) const {
    return s == other.s && token == other.token && leaves == other.leaves;
  }
  
  size_t hash() const {
    if (!hash_computed) {
      cached_hash = std::hash<int>{}(s);
      cached_hash ^= std::hash<int>{}(token) + 0x9e3779b9 + (cached_hash << 6) + (cached_hash >> 2);
      for (int x : leaves) {
        cached_hash ^= std::hash<int>{}(x) + 0x9e3779b9 + (cached_hash << 6) + (cached_hash >> 2);
      }
      hash_computed = true;
    }
    return cached_hash;
  }
};

struct LogPKeyHash {
  size_t operator()(const LogPKey& k) const noexcept {
    return k.hash();
  }
};

// ----- Solver
class Solver {
  const Downpass& D;
  const TokenPairs& pairs;
  LogRDCache& logRD;
  int presentBits;
  
  std::unordered_map<std::vector<int>, std::vector<double>, vec_int_hash> logB_cache;
  std::unordered_map<LogPKey, double, LogPKeyHash> logP_cache;
  
  inline bool all_equal(const std::vector<int>& a, const std::vector<int>& b) const {
    for (size_t i = 0; i < a.size(); ++i) if (a[i] != b[i]) return false;
    return true;
  }
  
  double LogB(int token0, const std::vector<int>& leaves) {
    if ((token_mask(token0) & ~presentBits) != 0) return NEG_INF;
    
    auto it = logB_cache.find(leaves);
    if (it == logB_cache.end()) {
      std::vector<double> v(D.nStates, std::numeric_limits<double>::quiet_NaN());
      it = logB_cache.emplace(leaves, std::move(v)).first;
    }
    std::vector<double>& row = it->second;
    double& slot = row[token0];
    if (R_finite(slot)) return slot;
    if (ISNAN(slot) == false && !R_finite(slot)) return slot;
    
    const int n = sum_int(leaves);
    if (n == 1) {
      slot = (leaves[token0] == 1) ? 0.0 : NEG_INF;
      return slot;
    }
    if (n == 2) {
      int twice = -1, a = -1, b = -1;
      for (int i = 0; i < (int)leaves.size(); ++i) {
        if (leaves[i] == 2) twice = i;
        else if (leaves[i] == 1) { if (a < 0) a = i; else b = i; }
      }
      int resultTok = (twice >= 0) ? twice : (D.dp_at(a, b) - 1);
      slot = (token0 == resultTok) ? 0.0 : NEG_INF;
      return slot;
    }
    
    auto &validDraws = GlobalValidDraws::get();
    const auto& drawpairs = validDraws.get(leaves);
    std::vector<double> terms;
    terms.reserve(drawpairs.size());
    
    for (const auto& dp : drawpairs) {
      const std::vector<int>& drawn = dp.drawn;
      const std::vector<int>& undrawn = dp.undrawn;
      int m = sum_int(drawn);
      
      double balancedCorrection = (2*m == n) ? std::log(2.0) : 0.0;
      if (all_equal(drawn, undrawn)) balancedCorrection -= std::log(2.0);
      
      // accumulate inner log-sum-exp without allocating an inner vector
      LSEAccumulator innerAcc;
      for (const auto& pr : pairs.noStep[token0]) {
        double val = LogB(pr.a, drawn) + LogB(pr.b, undrawn);
        innerAcc.add((long double)val);
      }
      for (const auto& pr : pairs.yesStep[token0]) {
        double val = LogB(pr.a, drawn) + LogB(pr.b, undrawn);
        innerAcc.add((long double)val);
      }
      double innerSum = innerAcc.result();
      
      double acc = balancedCorrection + logRD.compute(drawn, leaves) + innerSum;
      terms.push_back(acc);
    }
    slot = log_sum_exp(terms);
    return slot;
    
  }
  
  double LogP(int s, const std::vector<int>& leaves, int token0) {
    LogPKey key(s, token0, leaves);
    auto it = logP_cache.find(key);
    if (it != logP_cache.end()) return it->second;
    
    const int n = sum_int(leaves);
    
    if (n == 1) {
      double out = (leaves[token0] == 1) ? ((s == 0) ? 0.0 : NEG_INF) : NEG_INF;
      logP_cache.emplace(std::move(key), out);
      return out;
    }
    if (n == 2) {
      int twice = -1;
      int a = -1, b = -1;
      for (int i = 0; i < (int)leaves.size(); ++i) {
        if (leaves[i] == 2) twice = i;
        else if (leaves[i] == 1) { if (a < 0) a = i; else b = i; }
      }
      double out;
      if (twice >= 0) {
        out = (s == 0) ? 0.0 : NEG_INF;
      } else {
        bool stepAdds = D.step_at(a, b);
        int needed = stepAdds ? 1 : 0;
        out = (s == needed) ? 0.0 : NEG_INF;
      }
      logP_cache.emplace(std::move(key), out);
      return out;
    }
    
    double denom = LogB(token0, leaves);
    if (!R_finite(denom)) {
      logP_cache.emplace(std::move(key), denom);
      return denom;
    }
    
    auto &validDraws = GlobalValidDraws::get();
    const auto& drawpairs = validDraws.get(leaves);
    std::vector<double> outerTerms;
    outerTerms.reserve(drawpairs.size());
    
    for (const auto& dp : drawpairs) {
      const std::vector<int>& drawn = dp.drawn;
      const std::vector<int>& undrawn = dp.undrawn;
      const int m = sum_int(drawn);
      
      double sizeCorrection = ((m + m == n) && !all_equal(drawn, undrawn)) ? std::log(2.0) : 0.0;
      
      // noStep: accumulate over r (0..s)
      LSEAccumulator noStepAcc; // accumulates log-sum-exp over r of (log_sum_exp over pairs)
      for (int r = 0; r <= s; ++r) {
        // for this r, compute log-sum-exp across pairs (L)
        LSEAccumulator pairAcc;
        const auto& L = pairs.noStep[token0];
        for (const auto& pr : L) {
          double t = log_prod_sum({
            LogP(r, drawn, pr.a),
            LogB(pr.a, drawn),
            LogP(s - r, undrawn, pr.b),
            LogB(pr.b, undrawn)
          });
          pairAcc.add((long double)t);
        }
        double val_r = pairAcc.result();
        noStepAcc.add((long double)val_r);
      }
      double noStepSum = noStepAcc.result();
      
      // yesStep: accumulate over r (0..s-1) if applicable
      double yesStepSum = NEG_INF;
      if (!pairs.yesStep[token0].empty() && s >= 1) {
        LSEAccumulator yesStepAcc;
        for (int r = 0; r <= s - 1; ++r) {
          LSEAccumulator pairAcc;
          const auto& L2 = pairs.yesStep[token0];
          for (const auto& pr : L2) {
            double t = log_prod_sum({
              LogP(r, drawn, pr.a),
              LogB(pr.a, drawn),
              LogP(s - r - 1, undrawn, pr.b),
              LogB(pr.b, undrawn)
            });
            pairAcc.add((long double)t);
          }
          double val_r = pairAcc.result();
          yesStepAcc.add((long double)val_r);
        }
        yesStepSum = yesStepAcc.result();
      }
      
      // combine noStepSum and yesStepSum (log-sum-exp of the two)
      LSEAccumulator bothAcc;
      bothAcc.add((long double)noStepSum);
      bothAcc.add((long double)yesStepSum);
      double combined = bothAcc.result();
      
      double inner = logRD.compute(drawn, leaves) + sizeCorrection + combined;
      outerTerms.push_back(inner);
    }
    
    double result = log_sum_exp(outerTerms) - denom;
    logP_cache.emplace(std::move(key), result);
    return result;
    
  }
  
public:
  Solver(const Downpass& D_, const TokenPairs& p, LogRDCache& rd, int presentBits_)
    : D(D_), pairs(p), logRD(rd), presentBits(presentBits_) {
    logB_cache.reserve(256);
    logP_cache.reserve(1024);
  }
  
  double run(int steps, const std::vector<int>& states) {
    std::vector<double> terms;
    terms.reserve(D.nStates);
    for (int token0 = 0; token0 < D.nStates; ++token0) {
      double b = LogB(token0, states);
      double p = LogP(steps, states, token0);
      terms.push_back(log_prod_sum({b, p}));
    }
    return log_sum_exp(terms);
  }
};

//' @rdname Carter1
//' @examples
//' # Number of trees with 2 steps for character 0011122
//' MaddisonSlatkin(2, c("0" = 2, "1" = 3, "01" = 0, "2" = 2)) * NUnrooted(7)
//' 
//' @export
// [[Rcpp::export]]
NumericVector MaddisonSlatkin(IntegerVector steps, IntegerVector states) {
  int len = states.size();
  if (len <= 0) stop("`states` must have positive length.");
  int nLevels = (int)std::floor(std::log2((double)len)) + 1;
  int nStates = (1 << nLevels) - 1;
  
  std::vector<int> leaves(nStates, 0);
  for (int i = 0; i < std::min(len, nStates); ++i) {
   int v = states[i];
   if (IntegerVector::is_na(v)) v = 0;
   if (v < 0) stop("`states` must be non-negative counts.");
   leaves[i] = v;
  }
  
  int nTaxa = sum_int(leaves);
  LnRootedCache lnRooted(nTaxa);
  LogRDCache logRD(lnRooted);
  Downpass D(nLevels);
  
  int presentBits = 0;
  for (int t = 0; t < nStates; ++t) if (leaves[t] > 0) presentBits |= token_mask(t);
  TokenPairs pairs(D, presentBits);
  
  Solver solver(D, pairs, logRD, presentBits);
  
  int k = steps.size();
  NumericVector out(k);
  for (int i = 0; i < k; ++i) {
    int s = steps[i];
    if (IntegerVector::is_na(s))
      out[i] = NA_REAL;
    else
      out[i] = solver.run(s, leaves);
  }
  
  return out;
}

//' @export
 //' @keywords internal
 // [[Rcpp::export]]
 NumericVector MaddisonSlatkin_steps(int s_min, int s_max, IntegerVector states) {
   if (s_min < 0 || s_max < s_min) stop("Invalid steps range.");
   
   int len = states.size();
   int nLevels = (int)std::floor(std::log2((double)len)) + 1;
   int nStates = (1 << nLevels) - 1;
   
   std::vector<int> leaves(nStates, 0);
   for (int i = 0; i < std::min(len, nStates); ++i) {
     int v = states[i];
     if (IntegerVector::is_na(v)) v = 0;
     if (v < 0) stop("`states` must be non-negative counts.");
     leaves[i] = v;
   }
   
   int nTaxa = sum_int(leaves);
   LnRootedCache lnRooted(nTaxa);
   LogRDCache logRD(lnRooted);
   Downpass D(nLevels);
   int presentBits = 0;
   for (int t = 0; t < nStates; ++t) if (leaves[t] > 0) presentBits |= token_mask(t);
   TokenPairs pairs(D, presentBits);
   
   Solver solver(D, pairs, logRD, presentBits);
   
   int K = s_max - s_min + 1;
   NumericVector out(K);
   for (int k = 0; k < K; ++k) {
     out[k] = solver.run(s_min + k, leaves);
   }
   return out;
 }
