// MaddisonSlatkin.cpp
#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include <memory>
#include <cstring>
#include <array>

using namespace Rcpp;

static const double NEG_INF = -std::numeric_limits<double>::infinity();

// ============================================================================
// Fixed-size Key to avoid malloc in Map lookups
// Supports up to 16 states and 65535 taxa per state.
// Fits in 32 bytes (friendly to CPU cache lines).
// ============================================================================
constexpr int MAX_STATES_OPT = 16;
struct StateKey {
  // Total size required for 16 states (32 bytes) + len (1 byte) = 33 bytes.
  // We aim for 32 or 40 bytes (multiple of 8/32).
  
  // Data (32 bytes)
  uint16_t data[MAX_STATES_OPT];
  
  // Length (1 byte)
  uint8_t len;
  
  // Padding to ensure consistent size/hash: 
  // (We need at least 7 bytes of padding if we want 40 bytes total, 
  // or rely on implicit padding for 32 bytes if len is moved to the start).
  // Let's rely on memsetting to zero for hashing consistency and define explicit padding for 40 bytes.
  // Size = 16*2 + 1 = 33 bytes. Padding 7 bytes to get 40 bytes.
  uint8_t padding[7]; 
  
  StateKey() {
    std::memset(this, 0, sizeof(StateKey));
    len = 0; // Explicitly initialize len
  }
  
  // Fast construction from vector
  explicit StateKey(const std::vector<int>& v) {
    // We memset 0 to ensure padding bytes are 0 for hashing
    std::memset(this, 0, sizeof(StateKey));
    len = (uint8_t)v.size();
    for(size_t i=0; i<v.size(); ++i) {
      data[i] = (uint16_t)v[i];
    }
  }
  
  // Construct from two StateKeys (subtraction)
  StateKey(const StateKey& total, const StateKey& drawn) {
    std::memset(this, 0, sizeof(StateKey));
    len = total.len;
    for(int i=0; i<len; ++i) {
      data[i] = total.data[i] - drawn.data[i];
    }
  }
  
  bool operator==(const StateKey& other) const {
    // Total size is now 40 bytes
    return std::memcmp(this, &other, sizeof(StateKey)) == 0;
  }
  
  inline int sum() const {
    int s = 0;
    for(int i=0; i<len; ++i) s += data[i];
    return s;
  }
  
  inline int get(int idx) const { return data[idx]; }
  
  bool matches_vec(const std::vector<int>& v) const {
    if (v.size() != len) return false;
    for(size_t i=0; i<v.size(); ++i) if(v[i] != data[i]) return false;
    return true;
  }
};

// StateKeyHash remains the same, relying on sizeof(StateKey)
struct StateKeyHash {
  std::size_t operator()(const StateKey& k) const noexcept {
    uint64_t hash = 14695981039346656037ULL;
    const unsigned char* p = reinterpret_cast<const unsigned char*>(&k);
    // Hash the bytes of the entire struct (now 40 bytes)
    for (size_t i = 0; i < sizeof(StateKey); ++i) {
      hash ^= p[i];
      hash *= 1099511628211ULL;
    }
    return (std::size_t)hash;
  }
};

// Key for LogP cache
struct LogPKeyOpt {
  StateKey leaves;
  int s;
  int token;
  
  LogPKeyOpt(int s_, int t_, const StateKey& l_) : leaves(l_), s(s_), token(t_) {}
  
  bool operator==(const LogPKeyOpt& other) const {
    return s == other.s && token == other.token && leaves == other.leaves;
  }
};

struct LogPKeyOptHash {
  std::size_t operator()(const LogPKeyOpt& k) const noexcept {
    // Combine hash of StateKey with s and token
    size_t h = StateKeyHash{}(k.leaves);
    h ^= std::hash<int>{}(k.s) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= std::hash<int>{}(k.token) + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
  }
};


// ============================================================================
// Utilities
// ============================================================================

struct LSEAccumulator {
  long double maxv;
  long double acc; 
  bool empty;
  
  LSEAccumulator() : maxv(-std::numeric_limits<long double>::infinity()), acc(0.0L), empty(true) {}
  
  inline void add(long double v) {
    if (!std::isfinite((double)v)) return;
    
    if (empty) {
      maxv = v;
      acc = 1.0L;
      empty = false;
    } else {
      if (v > maxv) {
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
};

static inline double log_prod_sum(const std::vector<double>& xs) {
  long double s = 0.0L;
  for (double v : xs) {
    if (!R_finite(v)) return NEG_INF;
    s += v;
  }
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
    // Fallback for unexpected size
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
// Note: ValidDraws still returns vectors to be compatible with recursion logic,
// but the cache key is now StateKey.
struct DrawPair {
  StateKey drawn;
  StateKey undrawn;
  // We keep basic int counts for logic checks if needed, but primarily use keys
  // Storing pre-calculated info helps avoid lookups
  int m; 
};

class ValidDrawsCache {
  std::unordered_map<StateKey, std::vector<DrawPair>, StateKeyHash> cache;
  
  // Helper to convert StateKey back to vector for generation (only done once per key)
  std::vector<int> keyToVec(const StateKey& k) {
    std::vector<int> v(k.len);
    for(int i=0; i<k.len; ++i) v[i] = k.data[i];
    return v;
  }
  
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
      
      DrawPair dp;
      dp.drawn = StateKey(drawn);
      
      // Create undrawn directly from keys
      StateKey total(leaves);
      dp.undrawn = StateKey(total, dp.drawn);
      dp.m = curSum;
      
      out.push_back(dp);
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
  void clear() {
    cache.clear();
  }
  
  const std::vector<DrawPair>& get(const StateKey& leavesKey) {
    auto it = cache.find(leavesKey);
    if (it != cache.end()) return it->second;
    
    std::vector<int> leaves = keyToVec(leavesKey);
    int n = sum_int(leaves);
    int half = n / 2;
    
    std::vector<DrawPair> out;
    std::vector<int> drawn(leaves.size(), 0);
    rec(leaves, drawn, 0, n, half, 0, out);
    
    auto ins = cache.emplace(leavesKey, std::move(out));
    return ins.first->second;
  }
};

static ValidDrawsCache VALID_DRAWS_GLOBAL;

// ----- LogRD using StateKey
class LogRDCache {
  LnRootedCache& lnRooted;
  // LogRDKey is essentially (Drawn, Total), but since Total is constant in a context...
  // We combine two StateKeys.
  struct LogRDKeyOpt {
    StateKey drawn;
    StateKey leaves;
    bool operator==(const LogRDKeyOpt& o) const { return drawn == o.drawn && leaves == o.leaves; }
  };
  struct LogRDKeyOptHash {
    size_t operator()(const LogRDKeyOpt& k) const {
      return StateKeyHash{}(k.drawn) ^ (StateKeyHash{}(k.leaves) << 1);
    }
  };
  
  std::unordered_map<LogRDKeyOpt, double, LogRDKeyOptHash> cache;
  
public:
  explicit LogRDCache(LnRootedCache& lnr) : lnRooted(lnr) {
    cache.reserve(1024);
  }
  
  double compute(const StateKey& drawn, const StateKey& leaves) {
    LogRDKeyOpt key{drawn, leaves};
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;
    
    const int m = drawn.sum();
    const int n = leaves.sum();
    
    double bal = (n == 2*m) ? std::log(0.5) : 0.0;
    long double lc = 0.0L;
    for (int i = 0; i < leaves.len; ++i) {
      lc += lchoose_log(leaves.data[i], drawn.data[i]);
    }
    double val = bal + lnRooted(m) + lnRooted(n - m) - lnRooted(n) + (double)lc;
    cache.emplace(std::move(key), val);
    return val;
  }
};

// Global caches indexed by nTaxa
static std::unordered_map<int, std::shared_ptr<LogRDCache>> LOGRD_CACHE;
static std::unordered_map<int, std::shared_ptr<LnRootedCache>> LNROOT_CACHE;
static std::unordered_map<int, std::shared_ptr<Downpass>> DP_CACHE;
static std::unordered_map<long long, std::shared_ptr<TokenPairs>> TP_CACHE;

// Solver caches using Optimized Keys
struct SolverCaches {
  std::unordered_map<StateKey, std::vector<double>, StateKeyHash> logB_cache;
  std::unordered_map<LogPKeyOpt, double, LogPKeyOptHash> logP_cache;
};

// Key to solver cache is still the root config (packed leaves)
static std::unordered_map<uint64_t, std::shared_ptr<SolverCaches>> SOLVER_CACHE;

inline uint64_t pack_leaves(const std::vector<int> &v) {
  uint64_t h = 146527;
  for (int x : v) {
    uint64_t z = (uint64_t)x;
    z ^= z >> 33; z *= 0xff51afd7ed558ccdULL;
    z ^= z >> 33; z *= 0xc4ceb9fe1a85ec53ULL;
    z ^= z >> 33;
    h ^= z + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
  }
  return h;
}

// ----- Solver
class Solver {
  const Downpass& D;
  const TokenPairs& pairs;
  ValidDrawsCache& validDraws;
  LogRDCache& logRD;
  int presentBits;
  
  std::unordered_map<StateKey, std::vector<double>, StateKeyHash>& logB_cache;
  std::unordered_map<LogPKeyOpt, double, LogPKeyOptHash>& logP_cache;
  
  double LogB(int token0, const StateKey& leaves) {
    // Check Trivial Cases BEFORE HashMap Lookup.
    
    // n is computed fast via StateKey
    const int n = leaves.sum(); 
    
    if (n == 1) {
      return (leaves.get(token0) == 1) ? 0.0 : NEG_INF;
    }
    
    if (n == 2) {
      // Unroll logic for n=2 without vector iterators
      int twice = -1, a = -1, b = -1;
      for (int i = 0; i < leaves.len; ++i) {
        int count = leaves.get(i);
        if (count == 2) { twice = i; break; } // Optimization: break early
        else if (count == 1) { if (a < 0) a = i; else b = i; }
      }
      int resultTok = (twice >= 0) ? twice : (D.dp_at(a, b) - 1);
      return (token0 == resultTok) ? 0.0 : NEG_INF;
    }
    
    // Check if token allowed
    if ((token_mask(token0) & ~presentBits) != 0) return NEG_INF;
    
    // Only NOW do we check the cache
    auto it = logB_cache.find(leaves);
    if (it == logB_cache.end()) {
      std::vector<double> v(D.nStates, std::numeric_limits<double>::quiet_NaN());
      it = logB_cache.emplace(leaves, std::move(v)).first;
    }
    
    double& slot = it->second[token0];
    if (R_finite(slot)) return slot;
    if (ISNAN(slot) == false && !R_finite(slot)) return slot;
    
    const auto& drawpairs = validDraws.get(leaves);
    
    LSEAccumulator outerAcc;
    
    for (const auto& dp : drawpairs) {
      const StateKey& drawn = dp.drawn;
      const StateKey& undrawn = dp.undrawn;
      int m = dp.m;
      
      double balancedCorrection = (2*m == n) ? std::log(2.0) : 0.0;
      if (drawn == undrawn) balancedCorrection -= std::log(2.0);
      
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
      outerAcc.add((long double)acc);
    }
    slot = outerAcc.result();
    return slot;
  }
  
  double LogP(int s, const StateKey& leaves, int token0) {
    // OPTIMIZATION 2: Logic Inversion (Base cases first)
    const int n = leaves.sum();
    
    if (n == 1) {
      return (leaves.get(token0) == 1) ? ((s == 0) ? 0.0 : NEG_INF) : NEG_INF;
    }
    
    if (n == 2) {
      int twice = -1, a = -1, b = -1;
      for (int i = 0; i < leaves.len; ++i) {
        int count = leaves.get(i);
        if (count == 2) { twice = i; break; }
        else if (count == 1) { if (a < 0) a = i; else b = i; }
      }
      if (twice >= 0) {
        return (s == 0) ? 0.0 : NEG_INF;
      } else {
        bool stepAdds = D.step_at(a, b);
        int needed = stepAdds ? 1 : 0;
        return (s == needed) ? 0.0 : NEG_INF;
      }
    }
    
    // Cache Check
    LogPKeyOpt key(s, token0, leaves);
    auto it = logP_cache.find(key);
    if (it != logP_cache.end()) return it->second;
    
    double denom = LogB(token0, leaves);
    if (!R_finite(denom)) {
      logP_cache.emplace(std::move(key), denom);
      return denom;
    }
    
    const auto& drawpairs = validDraws.get(leaves);
    
    LSEAccumulator outerAcc;
    
    for (const auto& dp : drawpairs) {
      const StateKey& drawn = dp.drawn;
      const StateKey& undrawn = dp.undrawn;
      const int m = dp.m;
      
      double sizeCorrection = ((m + m == n) && !(drawn == undrawn)) ? std::log(2.0) : 0.0;
      
      LSEAccumulator noStepAcc;
      for (int r = 0; r <= s; ++r) {
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
        noStepAcc.add((long double)pairAcc.result());
      }
      double noStepSum = noStepAcc.result();
      
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
          yesStepAcc.add((long double)pairAcc.result());
        }
        yesStepSum = yesStepAcc.result();
      }
      
      LSEAccumulator bothAcc;
      bothAcc.add((long double)noStepSum);
      bothAcc.add((long double)yesStepSum);
      double combined = bothAcc.result();
      
      double inner = logRD.compute(drawn, leaves) + sizeCorrection + combined;
      outerAcc.add((long double)inner);
    }
    
    double result = outerAcc.result() - denom;
    logP_cache.emplace(std::move(key), result);
    return result;
  }
  
public:
  Solver(const Downpass& D_, const TokenPairs& p,
         ValidDrawsCache& vd, LogRDCache& rd, int presentBits_,
         SolverCaches& SC)
    : D(D_), pairs(p), validDraws(vd), logRD(rd), presentBits(presentBits_),
      logB_cache(SC.logB_cache), logP_cache(SC.logP_cache) {
    logB_cache.reserve(256);
    logP_cache.reserve(1024);
  }
  
  double run(int steps, const StateKey& states) {
    LSEAccumulator acc;
    for (int token0 = 0; token0 < D.nStates; ++token0) {
      double b = LogB(token0, states);
      double p = LogP(steps, states, token0);
      acc.add((long double)log_prod_sum({b, p}));
    }
    return acc.result();
  }
};

//' @rdname Carter1
//' @examples
//' MaddisonSlatkin(2, c("0" = 2, "1" = 3, "01" = 0, "2" = 2)) * NUnrooted(7)
//' 
//' @export
// [[Rcpp::export]]
NumericVector MaddisonSlatkin(IntegerVector steps, IntegerVector states) {
  int len = states.size();
  if (len <= 0) stop("`states` must have positive length.");
  
  // Performance Guard:
  if (len > MAX_STATES_OPT) {
    stop("This implementation supports a maximum of 16 distinct state tokens.");
  }
  
  int nLevels = (int)std::floor(std::log2((double)len)) + 1;
  int nStates = (1 << nLevels) - 1;
  
  std::vector<int> leavesVec(nStates, 0);
  for (int i = 0; i < std::min(len, nStates); ++i) {
    int v = states[i];
    if (IntegerVector::is_na(v)) v = 0;
    if (v < 0) stop("`states` must be non-negative counts.");
    leavesVec[i] = v;
  }
  
  // Create optimized key structure
  StateKey rootKey(leavesVec);
  
  int nTaxa = rootKey.sum();
  
  // --- Setup Shared Pointers (Same as before) ---
  std::shared_ptr<Downpass> Dptr;
  {
    auto it = DP_CACHE.find(nLevels);
    if (it != DP_CACHE.end()) Dptr = it->second;
    else {
      Dptr = std::make_shared<Downpass>(nLevels);
      DP_CACHE[nLevels] = Dptr;
    }
  }
  Downpass& D = *Dptr;
  
  int presentBits = 0;
  for (int t = 0; t < nStates; ++t) if (leavesVec[t] > 0) presentBits |= token_mask(t);
  
  long long tp_key = ((long long)nLevels << 20) | presentBits;
  std::shared_ptr<TokenPairs> Tpptr;
  {
    auto it = TP_CACHE.find(tp_key);
    if (it != TP_CACHE.end()) Tpptr = it->second;
    else {
      Tpptr = std::make_shared<TokenPairs>(D, presentBits);
      TP_CACHE[tp_key] = Tpptr;
    }
  }
  TokenPairs& pairs = *Tpptr;
  
  std::shared_ptr<LnRootedCache> LNRptr;
  {
    auto it = LNROOT_CACHE.find(nTaxa);
    if (it != LNROOT_CACHE.end()) LNRptr = it->second;
    else {
      LNRptr = std::make_shared<LnRootedCache>(nTaxa);
      LNROOT_CACHE[nTaxa] = LNRptr;
    }
  }
  LnRootedCache& lnRooted = *LNRptr;
  
  std::shared_ptr<LogRDCache> LRptr;
  {
    auto it = LOGRD_CACHE.find(nTaxa);
    if (it != LOGRD_CACHE.end()) LRptr = it->second;
    else {
      LRptr = std::make_shared<LogRDCache>(lnRooted);
      LOGRD_CACHE[nTaxa] = LRptr;
    }
  }
  LogRDCache& logRD = *LRptr;
  
  ValidDrawsCache& validDraws = VALID_DRAWS_GLOBAL;
  
  std::shared_ptr<SolverCaches> sc;
  uint64_t key = pack_leaves(leavesVec);
  
  {
    auto it = SOLVER_CACHE.find(key);
    if (it != SOLVER_CACHE.end()) sc = it->second;
    else {
      sc = std::make_shared<SolverCaches>();
      SOLVER_CACHE[key] = sc;
    }
  }
  
  Solver solver(D, pairs, validDraws, logRD, presentBits, *sc);
  
  int k = steps.size();
  NumericVector out(k);
  for (int i = 0; i < k; ++i) {
    int s = steps[i];
    if (IntegerVector::is_na(s))
      out[i] = NA_REAL;
    else
      out[i] = solver.run(s, rootKey);
  }
  
  return out;
}
//' @export
//' @keywords internal
// [[Rcpp::export]]
NumericVector MaddisonSlatkin_steps(int s_min, int s_max, IntegerVector states) {
  if (s_min < 0 || s_max < s_min) stop("Invalid steps range.");
  
  int len = states.size();
  
  // Performance Guard
  if (len > MAX_STATES_OPT) {
    stop("This implementation supports a maximum of 16 distinct state tokens.");
  }
  
  int nLevels = (int)std::floor(std::log2((double)len)) + 1;
  int nStates = (1 << nLevels) - 1;
  
  // 1. Prepare leaves vector for key creation
  std::vector<int> leavesVec(nStates, 0);
  for (int i = 0; i < std::min(len, nStates); ++i) {
    int v = states[i];
    if (IntegerVector::is_na(v)) v = 0;
    if (v < 0) stop("`states` must be non-negative counts.");
    leavesVec[i] = v;
  }
  
  // --- Key Creation ---
  // Use the optimized key structure for the solver run
  StateKey rootKey(leavesVec);
  int nTaxa = rootKey.sum();
  // --------------------
  
  // --- Downpass (persistent)
  std::shared_ptr<Downpass> Dptr;
  {
    auto it = DP_CACHE.find(nLevels);
    if (it != DP_CACHE.end()) Dptr = it->second;
    else {
      Dptr = std::make_shared<Downpass>(nLevels);
      DP_CACHE[nLevels] = Dptr;
    }
  }
  Downpass& D = *Dptr;
  
  int presentBits = 0;
  for (int t = 0; t < nStates; ++t) if (leavesVec[t] > 0) presentBits |= token_mask(t);
  
  // --- TokenPairs (persistent)
  long long tp_key = ((long long)nLevels << 20) | presentBits;
  std::shared_ptr<TokenPairs> Tpptr;
  {
    auto it = TP_CACHE.find(tp_key);
    if (it != TP_CACHE.end()) Tpptr = it->second;
    else {
      Tpptr = std::make_shared<TokenPairs>(D, presentBits);
      TP_CACHE[tp_key] = Tpptr;
    }
  }
  TokenPairs& pairs = *Tpptr;
  
  // --- LnRootedCache (persistent)
  std::shared_ptr<LnRootedCache> LNRptr;
  {
    auto it = LNROOT_CACHE.find(nTaxa);
    if (it != LNROOT_CACHE.end()) LNRptr = it->second;
    else {
      LNRptr = std::make_shared<LnRootedCache>(nTaxa);
      LNROOT_CACHE[nTaxa] = LNRptr;
    }
  }
  LnRootedCache& lnRooted = *LNRptr;
  
  // --- LogRDCache (persistent per nTaxa)
  std::shared_ptr<LogRDCache> LRptr;
  {
    auto it = LOGRD_CACHE.find(nTaxa);
    if (it != LOGRD_CACHE.end()) LRptr = it->second;
    else {
      LRptr = std::make_shared<LogRDCache>(lnRooted);
      LOGRD_CACHE[nTaxa] = LRptr;
    }
  }
  LogRDCache& logRD = *LRptr;
  
  // --- ValidDrawsCache is global
  ValidDrawsCache& validDraws = VALID_DRAWS_GLOBAL;
  
  // --- Solver-level caches per leaves
  std::shared_ptr<SolverCaches> sc;
  // NOTE: pack_leaves still requires the original vector for its key hashing algorithm
  uint64_t key = pack_leaves(leavesVec); 
  
  {
    auto it = SOLVER_CACHE.find(key);
    if (it != SOLVER_CACHE.end()) sc = it->second;
    else {
      sc = std::make_shared<SolverCaches>();
      SOLVER_CACHE[key] = sc;
    }
  }
  
  // Initialize solver using optimized caches
  Solver solver(D, pairs, validDraws, logRD, presentBits, *sc);
  
  int K = s_max - s_min + 1;
  NumericVector out(K);
  
  // Call the solver's run method using the optimized StateKey
  for (int k = 0; k < K; ++k) {
    out[k] = solver.run(s_min + k, rootKey);
  }
  
  return out;
}
//' @export
//' @keywords internal
// [[Rcpp::export]]
void MaddisonSlatkin_clear_cache() {
  SOLVER_CACHE.clear();
  LNROOT_CACHE.clear();
  LOGRD_CACHE.clear();
  DP_CACHE.clear();
  TP_CACHE.clear();
  VALID_DRAWS_GLOBAL.clear();
}