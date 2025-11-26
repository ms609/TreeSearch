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
constexpr int MAX_STATES_OPT = 32;
// ---------- replace existing StateKey definition ----------
struct StateKey {
  
  // cached 64-bit fingerprint (FNV-like). Mutable so maps can read it
  // even when keys are stored as const in std::unordered_map.
  mutable uint64_t cached_hash;
  uint16_t data[MAX_STATES_OPT];
  int cached_sum;
  uint8_t len;
  uint8_t padding[3];
  
  StateKey() : cached_sum(0), len(0), cached_hash(0) {
    std::memset(data, 0, sizeof(data));
    std::memset(padding, 0, sizeof(padding));
  }
  
  explicit StateKey(const std::vector<int>& v) : cached_sum(0), cached_hash(0) {
    std::memset(data, 0, sizeof(data));
    std::memset(padding, 0, sizeof(padding));
    len = (uint8_t)v.size();
    for (size_t i = 0; i < v.size(); ++i) {
      data[i] = (uint16_t)v[i];
      cached_sum += v[i];
    }
    // compute cached_hash once
    cached_hash = compute_hash_prefix();
  }
  
  StateKey(const StateKey& total, const StateKey& drawn) : cached_sum(0), cached_hash(0) {
    std::memset(data, 0, sizeof(data));
    std::memset(padding, 0, sizeof(padding));
    len = total.len;
    for (int i = 0; i < len; ++i) {
      data[i] = total.data[i] - drawn.data[i];
      cached_sum += data[i];
    }
    cached_hash = compute_hash_prefix();
  }
  
  inline int sum() const { return cached_sum; }
  inline int get(int idx) const { return data[idx]; }
  
  bool operator==(const StateKey& other) const {
    if (cached_sum != other.cached_sum) return false;
    if (len != other.len) return false;
    return std::memcmp(data, other.data, sizeof(uint16_t) * len) == 0;
  }
  
private:
  // FNV-1a-like over exactly the meaningful bytes (2*len bytes)
  uint64_t compute_hash_prefix() const noexcept {
    const uint8_t* bytes = reinterpret_cast<const uint8_t*>(data);
    const int used_bytes = len * sizeof(uint16_t);
    uint64_t h = 14695981039346656037ULL;
    constexpr uint64_t P = 1099511628211ULL;
    // simple byte loop (len <= 16 so this is tiny)
    for (int i = 0; i < used_bytes; ++i) {
      h ^= (uint64_t)bytes[i];
      h *= P;
    }
    // mix in metadata
    h = (h ^ (uint64_t)cached_sum) * P;
    h = (h ^ (uint64_t)len) * P;
    return h;
  }
};

struct StateKeyHash {
  std::size_t operator()(const StateKey& k) const noexcept {
    // Return the cached fingerprint (already mixed with sum/len).
    return (std::size_t)k.cached_hash;
  }
};


// ============================================================================
// Template-specialized StateKey for different token counts
// ============================================================================


// helper: core FNV compute over bytes (len*2 bytes) - repeated but tiny (len small)
static inline uint64_t compute_key_hash_from_bytes(const uint8_t* bytes, int used_bytes, int cached_sum, int len) noexcept {
  uint64_t h = 14695981039346656037ULL;
  constexpr uint64_t P = 1099511628211ULL;
  for (int i = 0; i < used_bytes; ++i) {
    h ^= (uint64_t)bytes[i];
    h *= P;
  }
  h = (h ^ (uint64_t)cached_sum) * P;
  h = (h ^ (uint64_t)len) * P;
  return h;
}
template<int nTokens>
struct StateKeyT;

// 2 tokens -> 3 states
template<>
struct StateKeyT<2> {
  mutable uint64_t cached_hash;
  uint16_t data[3];
  int cached_sum;
  uint8_t len;
  uint8_t padding;
  
  StateKeyT() : cached_sum(0), len(0), padding(0), cached_hash(0) {
    std::memset(data, 0, sizeof(data));
    cached_hash = compute_key_hash_from_bytes(reinterpret_cast<const uint8_t*>(data), 0, cached_sum, len);
  }
  
  explicit StateKeyT(const std::vector<int>& v) : cached_sum(0), padding(0), cached_hash(0) {
    std::memset(data, 0, sizeof(data));
    len = (uint8_t)v.size();
    for (size_t i = 0; i < v.size(); ++i) { data[i] = (uint16_t)v[i]; cached_sum += v[i]; }
    cached_hash = compute_key_hash_from_bytes(reinterpret_cast<const uint8_t*>(data), len * sizeof(uint16_t), cached_sum, len);
  }
  
  StateKeyT(const StateKeyT& total, const StateKeyT& drawn) : cached_sum(0), padding(0), cached_hash(0) {
    std::memset(data, 0, sizeof(data));
    len = total.len;
    for (int i = 0; i < len; ++i) { data[i] = total.data[i] - drawn.data[i]; cached_sum += data[i]; }
    cached_hash = compute_key_hash_from_bytes(reinterpret_cast<const uint8_t*>(data), len * sizeof(uint16_t), cached_sum, len);
  }
  
  inline int sum() const { return cached_sum; }
  inline int get(int idx) const { return data[idx]; }
  
  bool operator==(const StateKeyT& other) const {
    if (cached_hash != other.cached_hash) return false;             // fast reject
    if (cached_sum != other.cached_sum) return false;
    if (len != other.len) return false;
    return std::memcmp(data, other.data, sizeof(uint16_t) * len) == 0; // compare meaningful bytes only
  }
  
  bool matches_vec(const std::vector<int>& v) const {
    if ((int)v.size() != (int)len) return false;
    for (size_t i = 0; i < v.size(); ++i) if (v[i] != data[i]) return false;
    return true;
  }
};

// 3 tokens -> 7 states
template<>
struct StateKeyT<3> {
  mutable uint64_t cached_hash;
  uint16_t data[7];
  int cached_sum;
  uint8_t len;
  uint8_t padding;
  
  StateKeyT() : cached_sum(0), len(0), padding(0), cached_hash(0) { std::memset(data,0,sizeof(data)); cached_hash = compute_key_hash_from_bytes(reinterpret_cast<const uint8_t*>(data), 0, cached_sum, len); }
  explicit StateKeyT(const std::vector<int>& v) : cached_sum(0), padding(0), cached_hash(0) {
    std::memset(data,0,sizeof(data));
    len = (uint8_t)v.size();
    for (size_t i=0;i<v.size();++i){ data[i]=(uint16_t)v[i]; cached_sum += v[i]; }
    cached_hash = compute_key_hash_from_bytes(reinterpret_cast<const uint8_t*>(data), len * sizeof(uint16_t), cached_sum, len);
  }
  StateKeyT(const StateKeyT& total, const StateKeyT& drawn) : cached_sum(0), padding(0), cached_hash(0) {
    std::memset(data,0,sizeof(data));
    len = total.len;
    for (int i=0;i<len;++i){ data[i] = total.data[i] - drawn.data[i]; cached_sum += data[i]; }
    cached_hash = compute_key_hash_from_bytes(reinterpret_cast<const uint8_t*>(data), len * sizeof(uint16_t), cached_sum, len);
  }
  inline int sum() const { return cached_sum; }
  inline int get(int idx) const { return data[idx]; }
  bool operator==(const StateKeyT& other) const {
    if (cached_hash != other.cached_hash) return false;
    if (cached_sum != other.cached_sum) return false;
    if (len != other.len) return false;
    return std::memcmp(data, other.data, sizeof(uint16_t) * len) == 0;
  }
  bool matches_vec(const std::vector<int>& v) const {
    if ((int)v.size() != (int)len) return false;
    for (size_t i=0;i<v.size();++i) if (v[i] != data[i]) return false;
    return true;
  }
};

// 4 tokens -> 15 states
template<>
struct StateKeyT<4> {
  mutable uint64_t cached_hash;
  uint16_t data[15];
  int cached_sum;
  uint8_t len;
  uint8_t padding;
  
  StateKeyT() : cached_sum(0), len(0), padding(0), cached_hash(0) { std::memset(data,0,sizeof(data)); cached_hash = compute_key_hash_from_bytes(reinterpret_cast<const uint8_t*>(data), 0, cached_sum, len); }
  explicit StateKeyT(const std::vector<int>& v) : cached_sum(0), padding(0), cached_hash(0) {
    std::memset(data,0,sizeof(data));
    len = (uint8_t)v.size();
    for (size_t i=0;i<v.size();++i){ data[i] = (uint16_t)v[i]; cached_sum += v[i]; }
    cached_hash = compute_key_hash_from_bytes(reinterpret_cast<const uint8_t*>(data), len * sizeof(uint16_t), cached_sum, len);
  }
  StateKeyT(const StateKeyT& total, const StateKeyT& drawn) : cached_sum(0), padding(0), cached_hash(0) {
    std::memset(data,0,sizeof(data));
    len = total.len;
    for (int i=0;i<len;++i){ data[i] = total.data[i] - drawn.data[i]; cached_sum += data[i]; }
    cached_hash = compute_key_hash_from_bytes(reinterpret_cast<const uint8_t*>(data), len * sizeof(uint16_t), cached_sum, len);
  }
  inline int sum() const { return cached_sum; }
  inline int get(int idx) const { return data[idx]; }
  bool operator==(const StateKeyT& other) const {
    if (cached_hash != other.cached_hash) return false;
    if (cached_sum != other.cached_sum) return false;
    if (len != other.len) return false;
    return std::memcmp(data, other.data, sizeof(uint16_t) * len) == 0;
  }
  bool matches_vec(const std::vector<int>& v) const {
    if ((int)v.size() != (int)len) return false;
    for (size_t i=0;i<v.size();++i) if (v[i] != data[i]) return false;
    return true;
  }
};

// 5 tokens -> 31 states
template<>
struct StateKeyT<5> {
  mutable uint64_t cached_hash;
  uint16_t data[31];
  int cached_sum;
  uint8_t len;
  uint8_t padding;
  
  StateKeyT() : cached_sum(0), len(0), padding(0), cached_hash(0) { std::memset(data,0,sizeof(data)); cached_hash = compute_key_hash_from_bytes(reinterpret_cast<const uint8_t*>(data), 0, cached_sum, len); }
  explicit StateKeyT(const std::vector<int>& v) : cached_sum(0), padding(0), cached_hash(0) {
    std::memset(data,0,sizeof(data));
    len = (uint8_t)v.size();
    for (size_t i=0;i<v.size();++i){ data[i] = (uint16_t)v[i]; cached_sum += v[i]; }
    cached_hash = compute_key_hash_from_bytes(reinterpret_cast<const uint8_t*>(data), len * sizeof(uint16_t), cached_sum, len);
  }
  StateKeyT(const StateKeyT& total, const StateKeyT& drawn) : cached_sum(0), padding(0), cached_hash(0) {
    std::memset(data,0,sizeof(data));
    len = total.len;
    for (int i=0;i<len;++i){ data[i] = total.data[i] - drawn.data[i]; cached_sum += data[i]; }
    cached_hash = compute_key_hash_from_bytes(reinterpret_cast<const uint8_t*>(data), len * sizeof(uint16_t), cached_sum, len);
  }
  inline int sum() const { return cached_sum; }
  inline int get(int idx) const { return data[idx]; }
  bool operator==(const StateKeyT& other) const {
    if (cached_hash != other.cached_hash) return false;
    if (cached_sum != other.cached_sum) return false;
    if (len != other.len) return false;
    return std::memcmp(data, other.data, sizeof(uint16_t) * len) == 0;
  }
  bool matches_vec(const std::vector<int>& v) const {
    if ((int)v.size() != (int)len) return false;
    for (size_t i=0;i<v.size();++i) if (v[i] != data[i]) return false;
    return true;
  }
};

template<int nTokens>
struct StateKeyHashT;

template<>
struct StateKeyHashT<2> {
  std::size_t operator()(const StateKeyT<2>& k) const noexcept { return (std::size_t)k.cached_hash; }
};
template<>
struct StateKeyHashT<3> {
  std::size_t operator()(const StateKeyT<3>& k) const noexcept { return (std::size_t)k.cached_hash; }
};
template<>
struct StateKeyHashT<4> {
  std::size_t operator()(const StateKeyT<4>& k) const noexcept { return (std::size_t)k.cached_hash; }
};
template<>
struct StateKeyHashT<5> {
  std::size_t operator()(const StateKeyT<5>& k) const noexcept { return (std::size_t)k.cached_hash; }
};

// ============================================================================
// Fixed-size probability array instead of vector allocations
// Max states for 16 leaves = (1<<5)-1 = 31. We round to 32.
// ============================================================================
struct FixedProbList {
  double val[32]; // Fixed size, no malloc needed
  
  FixedProbList() {
    // Initialize with quiet_NaN to match original logic
    double nan = std::numeric_limits<double>::quiet_NaN();
    for(int i=0; i<32; ++i) val[i] = nan;
  }
  
  // Array access operator for convenience
  inline double& operator[](int idx) { return val[idx]; }
  inline const double& operator[](int idx) const { return val[idx]; }
};

// ============================================================================
// Custom Pool Allocator to kill malloc/free overhead
// This allocator is for the largest map (logP_cache)
// ============================================================================
template <typename T>
class MallocPoolAllocator {
private:
  static constexpr size_t BLOCK_SIZE = 1024 * 1024; // 1MB block size
  char* current_block = nullptr;
  char* current_pos = nullptr;
  char* end_pos = nullptr;
  std::vector<char*> blocks; // Stores pointers to all allocated blocks
  std::allocator<T> fallback_allocator;
  
  void allocate_new_block() {
    // Use standard C malloc/free for the large blocks
    current_block = (char*)std::malloc(BLOCK_SIZE);
    if (!current_block) throw std::bad_alloc();
    current_pos = current_block;
    end_pos = current_block + BLOCK_SIZE;
    blocks.push_back(current_block);
  }
  
public:
  using value_type = T;
  // Standard required typedefs and constructors
  MallocPoolAllocator() noexcept { allocate_new_block(); }
  template <class U> MallocPoolAllocator(const MallocPoolAllocator<U>&) noexcept : MallocPoolAllocator() {}
  
  // Destructor frees all memory blocks in one go
  ~MallocPoolAllocator() noexcept {
    for (char* block : blocks) {
      std::free(block);
    }
  }
  
  // Allocate method: Check n. Use pool for n=1, fallback for n > 1.
  T* allocate(size_t n) {
    if (n == 1) {
      // POOL ALLOCATION (for single node)
      size_t size = sizeof(T);
      
      // Ensure alignment:
      size_t alignment = alignof(T);
      size_t aligned_pos = (size_t)current_pos;
      size_t padding = (alignment - (aligned_pos % alignment)) % alignment;
      char* next_pos = current_pos + padding;
      
      if (next_pos + size > end_pos) {
        allocate_new_block();
        // Re-calculate alignment for the new block
        aligned_pos = (size_t)current_pos;
        padding = (alignment - (aligned_pos % alignment)) % alignment;
        next_pos = current_pos + padding;
      }
      
      current_pos = next_pos + size;
      return (T*)next_pos;
    } else {
      // FALLBACK ALLOCATION (for bucket arrays, etc.)
      // Delegates the large allocation to the standard C++ allocator
      return fallback_allocator.allocate(n);
    }
  }
  
  // Deallocate method: must delegate to the correct method
  void deallocate(T* p, size_t n) noexcept {
    if (n > 1) {
      // Deallocate memory allocated by the fallback allocator
      fallback_allocator.deallocate(p, n);
    } 
    // else n == 1: NO-OP, memory is pool-allocated and freed in the destructor.
  }
  
  // Required for C++ standard containers
  template <class U> struct rebind {
    using other = MallocPoolAllocator<U>;
  };
  bool operator!=(const MallocPoolAllocator& other) const noexcept { return this != &other; }
  bool operator==(const MallocPoolAllocator& other) const noexcept { return this == &other; }
};

// Key for LogP cache
struct LogPKeyOpt {
  StateKey leaves;
  int s;
  int token;
  
  LogPKeyOpt(int s_, int t_, const StateKey& l_) : leaves(l_), s(s_), token(t_) {}
  
  bool operator==(const LogPKeyOpt& other) const {
    // Fast check on integers first
    if (s != other.s) return false;
    if (token != other.token) return false;
    // Then delegate to StateKey's optimized comparison
    return leaves == other.leaves;
  }
};

struct LogPKeyOptHash {
  std::size_t operator()(const LogPKeyOpt& k) const noexcept {
    // 1. Get the fast StateKey hash
    uint64_t hash = StateKeyHash{}(k.leaves);
    
    // 2. Simple, fast incorporation of s and token using known mixing primes.
    // This is significantly faster than re-running the full FNV algorithm.
    
    // Hash for s: XOR and multiply with a good prime
    hash = (hash ^ (uint64_t)k.s) * 3935559000370003845ULL;
    
    // Hash for token: XOR and multiply with a different prime
    hash = (hash ^ (uint64_t)k.token) * 4488902095908611103ULL;
    
    // Final size_t cast is all that's needed.
    return (std::size_t)hash;
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

static inline double log_prod_sum_4(double v1, double v2, double v3, double v4) {
  if (!std::isfinite(v1) || !std::isfinite(v2) || 
      !std::isfinite(v3) || !std::isfinite(v4)) return NEG_INF;
      return (double)((long double)v1 + (long double)v2 + (long double)v3 + (long double)v4);
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
};// ============================================================================
// Optimization 3: Fixed-Size DrawPairs to eliminate vector allocation
// ============================================================================

struct DrawPair {
  StateKey drawn;
  StateKey undrawn;
  int m; 
};

struct FixedDraws {
  DrawPair draws[256];
  uint8_t count = 0; // Tracks the actual number of elements used
};

struct ValidDrawsCache {
  std::unordered_map<StateKey, FixedDraws, StateKeyHash> cache;
  
  // Helper to convert StateKey back to vector for generation (only done once per key)
  std::vector<int> keyToVec(const StateKey& k) {
    std::vector<int> v(k.len);
    for(int i=0; i<k.len; ++i) v[i] = k.data[i];
    return v;
  }
  
  void rec(const std::vector<int>& leaves, std::vector<int>& drawn, int idx,
           int n, int half, int curSum, FixedDraws& out) {
    
    if (idx == (int)leaves.size()) {
      if (curSum == 0 || curSum > half) return;
      if (curSum * 2 == n) {
        std::vector<int> undrawn_local(leaves.size());
        for (size_t i = 0; i < leaves.size(); ++i) undrawn_local[i] = leaves[i] - drawn[i];
        
        // Assuming lex_compare is defined elsewhere and handles symmetry breaking
        int cmp = lex_compare(drawn, undrawn_local); 
        if (cmp > 0) return;
      }
      
      // Safety check for the fixed array size (unlikely to be hit)
      if (out.count >= 256) throw std::runtime_error("FixedDraws array capacity exceeded.");
      
      DrawPair dp;
      dp.drawn = StateKey(drawn);
      
      // Create undrawn directly from keys
      StateKey total(leaves);
      dp.undrawn = StateKey(total, dp.drawn);
      dp.m = curSum;
      
      out.draws[out.count++] = dp; 
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
  
  const FixedDraws& get(const StateKey& leavesKey) {
    auto it = cache.find(leavesKey);
    if (it != cache.end()) return it->second;
    
    std::vector<int> leaves = keyToVec(leavesKey);
    int n = sum_int(leaves);
    int half = n / 2;
    
    FixedDraws out;
    
    std::vector<int> drawn(leaves.size(), 0);
    rec(leaves, drawn, 0, n, half, 0, out);
    
    // Use std::move to efficiently place 'out' into the cache
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
  // We kept FixedProbList for logB_cache from the last step
  std::unordered_map<StateKey, FixedProbList, StateKeyHash> logB_cache;
  
  // CHANGED: logP_cache now uses the custom pool allocator
  using LogPMapType = std::unordered_map<
    LogPKeyOpt, 
    double, 
    LogPKeyOptHash, 
    std::equal_to<LogPKeyOpt>, 
    MallocPoolAllocator<std::pair<const LogPKeyOpt, double>> 
  >;
  
  LogPMapType logP_cache;
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
// ============================================================================
// Template wrapper to dispatch to specialized implementations
// ============================================================================

template<int nTokens>
class SolverT {
  using KeyType = StateKeyT<nTokens>;
  using HashType = StateKeyHashT<nTokens>;
  
  const Downpass& D;
  const TokenPairs& pairs;
  int presentBits;
  
  // Specialized caches
  std::unordered_map<KeyType, FixedProbList, HashType> logB_cache;
  
  struct LogPKeyOpt {
    KeyType leaves;
    int s;
    int token;
    
    LogPKeyOpt(int s_, int t_, const KeyType& l_) : leaves(l_), s(s_), token(t_) {}
    
    bool operator==(const LogPKeyOpt& other) const {
      if (s != other.s) return false;
      if (token != other.token) return false;
      return leaves == other.leaves;
    }
  };
  
  struct LogPKeyOptHash {
    std::size_t operator()(const LogPKeyOpt& k) const noexcept {
      uint64_t hash = HashType{}(k.leaves);
      hash = (hash ^ (uint64_t)k.s) * 3935559000370003845ULL;
      hash = (hash ^ (uint64_t)k.token) * 4488902095908611103ULL;
      return (std::size_t)hash;
    }
  };
  
  using LogPMapType = std::unordered_map<
  LogPKeyOpt, 
  double, 
  LogPKeyOptHash, 
  std::equal_to<LogPKeyOpt>, 
  MallocPoolAllocator<std::pair<const LogPKeyOpt, double>>
    >;
    
    LogPMapType logP_cache;
    
    // Specialized ValidDraws cache
    struct DrawPairT {
      KeyType drawn;
      KeyType undrawn;
      int m;
    };
    
    struct FixedDrawsT {
      DrawPairT draws[256];
      uint8_t count = 0;
    };
    
    std::unordered_map<KeyType, FixedDrawsT, HashType> validDraws_cache;
    
    // LogRD cache
    struct LogRDKeyOpt {
      KeyType drawn;
      KeyType leaves;
      bool operator==(const LogRDKeyOpt& o) const { 
        return drawn == o.drawn && leaves == o.leaves; 
      }
    };
    
    struct LogRDKeyOptHash {
      size_t operator()(const LogRDKeyOpt& k) const {
        return HashType{}(k.drawn) ^ (HashType{}(k.leaves) << 1);
      }
    };
    
    std::unordered_map<LogRDKeyOpt, double, LogRDKeyOptHash> logRD_cache;
    LnRootedCache& lnRooted;
    
    // Helper methods (same logic as before, but with KeyType)
    void generateValidDraws(const KeyType& leavesKey, FixedDrawsT& out) {
      std::vector<int> leaves(leavesKey.len);
      for(int i=0; i<leavesKey.len; ++i) leaves[i] = leavesKey.data[i];
      
      int n = leavesKey.sum();
      int half = n / 2;
      std::vector<int> drawn(leaves.size(), 0);
      
      recDraws(leaves, drawn, 0, n, half, 0, out);
    }
    
    void recDraws(const std::vector<int>& leaves, std::vector<int>& drawn, 
                  int idx, int n, int half, int curSum, FixedDrawsT& out) {
      if (idx == (int)leaves.size()) {
        if (curSum == 0 || curSum > half) return;
        if (curSum * 2 == n) {
          std::vector<int> undrawn_local(leaves.size());
          for (size_t i = 0; i < leaves.size(); ++i) 
            undrawn_local[i] = leaves[i] - drawn[i];
          
          int cmp = lex_compare(drawn, undrawn_local);
          if (cmp > 0) return;
        }
        
        if (out.count >= 256) 
          throw std::runtime_error("FixedDraws array capacity exceeded.");
        
        DrawPairT dp;
        dp.drawn = KeyType(drawn);
        KeyType total(leaves);
        dp.undrawn = KeyType(total, dp.drawn);
        dp.m = curSum;
        
        out.draws[out.count++] = dp;
        return;
      }
      
      int maxTake = std::min(leaves[idx], half - curSum);
      for (int k = 0; k <= maxTake; ++k) {
        drawn[idx] = k;
        recDraws(leaves, drawn, idx + 1, n, half, curSum + k, out);
      }
      drawn[idx] = 0;
    }
    
    const FixedDrawsT& getValidDraws(const KeyType& leavesKey) {
      auto it = validDraws_cache.find(leavesKey);
      if (it != validDraws_cache.end()) return it->second;
      
      FixedDrawsT out;
      generateValidDraws(leavesKey, out);
      auto ins = validDraws_cache.emplace(leavesKey, std::move(out));
      return ins.first->second;
    }
    
    double computeLogRD(const KeyType& drawn, const KeyType& leaves) {
      LogRDKeyOpt key{drawn, leaves};
      auto it = logRD_cache.find(key);
      if (it != logRD_cache.end()) return it->second;
      
      const int m = drawn.sum();
      const int n = leaves.sum();
      
      double bal = (n == 2*m) ? std::log(0.5) : 0.0;
      long double lc = 0.0L;
      for (int i = 0; i < leaves.len; ++i) {
        lc += lchoose_log(leaves.data[i], drawn.data[i]);
      }
      double val = bal + lnRooted(m) + lnRooted(n - m) - lnRooted(n) + (double)lc;
      logRD_cache.emplace(std::move(key), val);
      return val;
    }
    
    double LogB(int token0, const KeyType& leaves) {
      const int n = leaves.sum();
      
      if (n == 1) {
        return (leaves.get(token0) == 1) ? 0.0 : NEG_INF;
      }
      
      if (n == 2) {
        int twice = -1, a = -1, b = -1;
        for (int i = 0; i < leaves.len; ++i) {
          int count = leaves.get(i);
          if (count == 2) { twice = i; break; }
          else if (count == 1) { if (a < 0) a = i; else b = i; }
        }
        int resultTok = (twice >= 0) ? twice : (D.dp_at(a, b) - 1);
        return (token0 == resultTok) ? 0.0 : NEG_INF;
      }
      
      if ((token_mask(token0) & ~presentBits) != 0) return NEG_INF;
      
      auto it = logB_cache.find(leaves);
      if (it == logB_cache.end()) {
        it = logB_cache.emplace(leaves, FixedProbList{}).first;
      }
      
      double& slot = it->second[token0];
      if (std::isfinite(slot)) return slot;
      if (std::isnan(slot) == false && !std::isfinite(slot)) return slot;
      
      const auto& drawpairs = getValidDraws(leaves);
      LSEAccumulator outerAcc;
      
      for (uint8_t i = 0; i < drawpairs.count; ++i) {
        const auto& dp = drawpairs.draws[i];
        const KeyType& drawn = dp.drawn;
        const KeyType& undrawn = dp.undrawn;
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
        
        double acc = balancedCorrection + computeLogRD(drawn, leaves) + innerSum;
        outerAcc.add((long double)acc);
      }
      slot = outerAcc.result();
      return slot;
    }
    
    double LogP(int s, const KeyType& leaves, int token0) {
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
      
      LogPKeyOpt key(s, token0, leaves);
      auto it = logP_cache.find(key);
      if (it != logP_cache.end()) return it->second;
      
      double denom = LogB(token0, leaves);
      if (!std::isfinite(denom)) {
        logP_cache.emplace(std::move(key), denom);
        return denom;
      }
      
      const auto& drawpairs = getValidDraws(leaves);
      LSEAccumulator outerAcc;
      
      for (uint8_t i = 0; i < drawpairs.count; ++i) {
        const auto& dp = drawpairs.draws[i];
        const KeyType& drawn = dp.drawn;
        const KeyType& undrawn = dp.undrawn;
        const int m = dp.m;
        
        double sizeCorrection = ((m + m == n) && !(drawn == undrawn)) ? std::log(2.0) : 0.0;
        
        LSEAccumulator noStepAcc;
        for (int r = 0; r <= s; ++r) {
          LSEAccumulator pairAcc;
          const auto& L = pairs.noStep[token0];
          for (const auto& pr : L) {
            double t = log_prod_sum_4(
              LogP(r, drawn, pr.a),
              LogB(pr.a, drawn),
              LogP(s - r, undrawn, pr.b),
              LogB(pr.b, undrawn)
            );
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
              double t = log_prod_sum_4(
                LogP(r, drawn, pr.a),
                LogB(pr.a, drawn),
                LogP(s - r - 1, undrawn, pr.b),
                LogB(pr.b, undrawn)
              );
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
        
        double inner = computeLogRD(drawn, leaves) + sizeCorrection + combined;
        outerAcc.add((long double)inner);
      }
      
      double result = outerAcc.result() - denom;
      logP_cache.emplace(std::move(key), result);
      return result;
    }
    
public:
  SolverT(const Downpass& D_, const TokenPairs& p, int presentBits_,
          LnRootedCache& lnr)
    : D(D_), pairs(p), presentBits(presentBits_), lnRooted(lnr) {
    
    logB_cache.reserve(16384);
    logB_cache.max_load_factor(0.5f);
    
    logP_cache.reserve(1 << 17);
    logP_cache.max_load_factor(0.5f);
    
    logRD_cache.reserve(1024);
    validDraws_cache.reserve(256);
  }
  
  double run(int steps, const KeyType& states) {
    LSEAccumulator acc;
    for (int token0 = 0; token0 < D.nStates; ++token0) {
      double b = LogB(token0, states);
      double p = LogP(steps, states, token0);
      double val;
      if (!std::isfinite(b) || !std::isfinite(p)) {
        val = NEG_INF;
      } else {
        val = b + p;
      }
      acc.add((long double)val);
    }
    return acc.result();
  }
};
// ----- Solver
class Solver {
  const Downpass& D;
  const TokenPairs& pairs;
  ValidDrawsCache& validDraws;
  LogRDCache& logRD;
  int presentBits;
  
  std::unordered_map<StateKey, FixedProbList, StateKeyHash>& logB_cache;
  SolverCaches::LogPMapType& logP_cache;
  
  double LogB(int token0, const StateKey& leaves) {
    // 1. Trivial Case Checks (Logic Inversion)
    const int n = leaves.sum(); 
    
    if (n == 1) {
      return (leaves.get(token0) == 1) ? 0.0 : NEG_INF;
    }
    
    if (n == 2) {
      int twice = -1, a = -1, b = -1;
      for (int i = 0; i < leaves.len; ++i) {
        int count = leaves.get(i);
        if (count == 2) { twice = i; break; }
        else if (count == 1) { if (a < 0) a = i; else b = i; }
      }
      int resultTok = (twice >= 0) ? twice : (D.dp_at(a, b) - 1);
      return (token0 == resultTok) ? 0.0 : NEG_INF;
    }
    
    // Check if token allowed
    if ((token_mask(token0) & ~presentBits) != 0) return NEG_INF;
    
    // 2. Cache Lookup
    // Note: emplace creates the FixedProbList (and its NaNs) automatically
    auto it = logB_cache.find(leaves);
    if (it == logB_cache.end()) {
      it = logB_cache.emplace(leaves, FixedProbList{}).first;
    }
    
    // 3. Check for computed value
    // Use operator[] on our new struct
    double& slot = it->second[token0]; 
    if (std::isfinite(slot)) return slot;
    if (std::isnan(slot) == false && !std::isfinite(slot)) return slot;
    
    // 4. Compute
    const auto& drawpairs = validDraws.get(leaves);
    LSEAccumulator outerAcc;
    
    for (uint8_t i = 0; i < drawpairs.count; ++i) { 
      const auto& dp = drawpairs.draws[i];
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
    if (!std::isfinite(denom)) {
      logP_cache.emplace(std::move(key), denom);
      return denom;
    }
    
    const auto& drawpairs = validDraws.get(leaves);
    
    LSEAccumulator outerAcc;
    
    for (uint8_t i = 0; i < drawpairs.count; ++i) { 
      const auto& dp = drawpairs.draws[i];
      const StateKey& drawn = dp.drawn;
      const StateKey& undrawn = dp.undrawn;
      const int m = dp.m;
      
      double sizeCorrection = ((m + m == n) && !(drawn == undrawn)) ? std::log(2.0) : 0.0;
      
      LSEAccumulator noStepAcc;
      for (int r = 0; r <= s; ++r) {
        LSEAccumulator pairAcc;
        const auto& L = pairs.noStep[token0];
        for (const auto& pr : L) {
          double t = log_prod_sum_4(
            LogP(r, drawn, pr.a),
            LogB(pr.a, drawn),
            LogP(s - r, undrawn, pr.b),
            LogB(pr.b, undrawn)
          );
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
            double t = log_prod_sum_4(
              LogP(r, drawn, pr.a),
              LogB(pr.a, drawn),
              LogP(s - r - 1, undrawn, pr.b),
              LogB(pr.b, undrawn)
            );
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
    
    logB_cache.reserve(16384);
    logB_cache.max_load_factor(0.7f);
    
    logP_cache.reserve(65536);
    logP_cache.max_load_factor(0.7f);
  }
  
  double run(int steps, const StateKey& states) {
    LSEAccumulator acc;
    for (int token0 = 0; token0 < D.nStates; ++token0) {
      double b = LogB(token0, states);
      double p = LogP(steps, states, token0);
      double val;
      if (!std::isfinite(p) || !std::isfinite(p)) {
        val = NEG_INF;
      } else {
        val = b + p;
      }
      acc.add((long double)val);
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
  
  // Determine nTokens from length
  int nTokens = (int)std::floor(std::log2((double)len)) + 1;
  
  if (nTokens < 2 || nTokens > 5) {
    stop("MaddisonSlatkin() supports 2-5 tokens (3-31 states).");
  }
  
  int nLevels = nTokens;
  int nStates = (1 << nLevels) - 1;
  
  std::vector<int> leavesVec(nStates, 0);
  for (int i = 0; i < std::min(len, nStates); ++i) {
    int v = states[i];
    if (IntegerVector::is_na(v)) v = 0;
    if (v < 0) stop("`states` must be non-negative counts.");
    leavesVec[i] = v;
  }
  
  int nTaxa = sum_int(leavesVec);
  
  // Setup shared structures
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
  for (int t = 0; t < nStates; ++t) {
    if (leavesVec[t] > 0) presentBits |= token_mask(t);
  }
    
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
  
  int k = steps.size();
  NumericVector out(k);
  
  // DISPATCH based on nTokens
  if (nTokens == 2) {
    StateKeyT<2> rootKey(leavesVec);
    SolverT<2> solver(D, pairs, presentBits, lnRooted);
    for (int i = 0; i < k; ++i) {
      int s = steps[i];
      if (IntegerVector::is_na(s))
        out[i] = NA_REAL;
      else
        out[i] = solver.run(s, rootKey);
    }
  } else if (nTokens == 3) {
    StateKeyT<3> rootKey(leavesVec);
    SolverT<3> solver(D, pairs, presentBits, lnRooted);
    for (int i = 0; i < k; ++i) {
      int s = steps[i];
      if (IntegerVector::is_na(s))
        out[i] = NA_REAL;
      else
        out[i] = solver.run(s, rootKey);
    }
  } else if (nTokens == 4) {
    StateKeyT<4> rootKey(leavesVec);
    SolverT<4> solver(D, pairs, presentBits, lnRooted);
    for (int i = 0; i < k; ++i) {
      int s = steps[i];
      if (IntegerVector::is_na(s))
        out[i] = NA_REAL;
      else
        out[i] = solver.run(s, rootKey);
    }
  } else { // nTokens == 5
    StateKeyT<5> rootKey(leavesVec);
    SolverT<5> solver(D, pairs, presentBits, lnRooted);
    for (int i = 0; i < k; ++i) {
      int s = steps[i];
      if (IntegerVector::is_na(s))
        out[i] = NA_REAL;
      else
        out[i] = solver.run(s, rootKey);
    }
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
