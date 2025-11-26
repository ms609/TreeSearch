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
  
  StateKey() : cached_hash(0), cached_sum(0), len(0) {
    std::memset(data, 0, sizeof(data));
    std::memset(padding, 0, sizeof(padding));
  }
  
  explicit StateKey(const std::vector<int>& v) : cached_hash(0), cached_sum(0) {
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
  
  StateKey(const StateKey& total, const StateKey& drawn) : cached_hash(0), cached_sum(0) {
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
  
  StateKeyT() : cached_hash(0), cached_sum(0), len(0), padding(0) {
    std::memset(data, 0, sizeof(data));
    cached_hash = compute_key_hash_from_bytes(reinterpret_cast<const uint8_t*>(data), 0, cached_sum, len);
  }
  
  explicit StateKeyT(const std::vector<int>& v) : cached_hash(0), cached_sum(0), padding(0) {
    std::memset(data, 0, sizeof(data));
    len = (uint8_t)v.size();
    for (size_t i = 0; i < v.size(); ++i) { data[i] = (uint16_t)v[i]; cached_sum += v[i]; }
    cached_hash = compute_key_hash_from_bytes(reinterpret_cast<const uint8_t*>(data), len * sizeof(uint16_t), cached_sum, len);
  }
  
  StateKeyT(const StateKeyT& total, const StateKeyT& drawn) : cached_hash(0), cached_sum(0), padding(0) {
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
  
  StateKeyT() : cached_hash(0), cached_sum(0), len(0), padding(0) {
    std::memset(data,0,sizeof(data)); cached_hash = compute_key_hash_from_bytes(reinterpret_cast<const uint8_t*>(data), 0, cached_sum, len); }
  explicit StateKeyT(const std::vector<int>& v) : cached_hash(0), cached_sum(0), padding(0) {
    std::memset(data,0,sizeof(data));
    len = (uint8_t)v.size();
    for (size_t i=0;i<v.size();++i){ data[i]=(uint16_t)v[i]; cached_sum += v[i]; }
    cached_hash = compute_key_hash_from_bytes(reinterpret_cast<const uint8_t*>(data), len * sizeof(uint16_t), cached_sum, len);
  }
  StateKeyT(const StateKeyT& total, const StateKeyT& drawn) : cached_hash(0), cached_sum(0), padding(0) {
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
  
  StateKeyT() : cached_hash(0), cached_sum(0), len(0), padding(0) {
    std::memset(data,0,sizeof(data));
    cached_hash = compute_key_hash_from_bytes(reinterpret_cast<const uint8_t*>(data), 0, cached_sum, len); }
  explicit StateKeyT(const std::vector<int>& v) : cached_hash(0), cached_sum(0), padding(0) {
    std::memset(data,0,sizeof(data));
    len = (uint8_t)v.size();
    for (size_t i=0;i<v.size();++i){ data[i] = (uint16_t)v[i]; cached_sum += v[i]; }
    cached_hash = compute_key_hash_from_bytes(reinterpret_cast<const uint8_t*>(data), len * sizeof(uint16_t), cached_sum, len);
  }
  StateKeyT(const StateKeyT& total, const StateKeyT& drawn) : cached_hash(0), cached_sum(0), padding(0) {
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
  
  StateKeyT() : cached_hash(0), cached_sum(0), len(0), padding(0) {
    std::memset(data,0,sizeof(data)); cached_hash = compute_key_hash_from_bytes(reinterpret_cast<const uint8_t*>(data), 0, cached_sum, len); }
  explicit StateKeyT(const std::vector<int>& v) : cached_hash(0), cached_sum(0), padding(0) {
    std::memset(data,0,sizeof(data));
    len = (uint8_t)v.size();
    for (size_t i=0;i<v.size();++i){ data[i] = (uint16_t)v[i]; cached_sum += v[i]; }
    cached_hash = compute_key_hash_from_bytes(reinterpret_cast<const uint8_t*>(data), len * sizeof(uint16_t), cached_sum, len);
  }
  StateKeyT(const StateKeyT& total, const StateKeyT& drawn) : cached_hash(0), cached_sum(0), padding(0) {
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

// ============================================================================
// PART 2: FIXED TREE TOPOLOGY SOLVER
// ============================================================================

// 1. Structure to hold the pre-processed tree topology
//    This allows us to traverse the tree efficiently without re-parsing R matrices.
struct FixedTreeNode {
  int id;           // 0-based index
  int left;         // Child index (or -1 if tip)
  int right;        // Child index (or -1 if tip)
  int nTips;        // Total tips in this subtree
  bool isTip;
};

class FixedTreeStructure {
public:
  std::vector<FixedTreeNode> nodes;
  int rootIdx;
  int nTaxa;
  
  FixedTreeStructure(IntegerMatrix edge, int nTips_in) : nTaxa(nTips_in) {
    int nEdges = edge.nrow();
    // Max node index in R is usually nTips + Nnode. 
    // We map R indices (1-based) to C++ (0-based).
    // R tips: 1..nTips. R internals: nTips+1..max.
    
    // Determine max index
    int maxId = 0;
    for(int i=0; i < nEdges*2; ++i) {
      if(edge[i] > maxId) maxId = edge[i];
    }
    
    nodes.resize(maxId + 1); // direct mapping for simplicity
    
    // Initialize
    for(int i=0; i<=maxId; ++i) {
      nodes[i].id = i;
      nodes[i].left = -1;
      nodes[i].right = -1;
      nodes[i].nTips = (i < nTips_in) ? 1 : 0; // Check R 0-based logic carefully below
      nodes[i].isTip = (i < nTips_in);         // 1..nTips in R -> 0..nTips-1 in C++
    }
    
    // Identify Root: The node that never appears in column 2 (child)
    std::vector<bool> isChild(maxId + 1, false);
    
    // Build Adjacency (assuming binary)
    for(int i=0; i < nEdges; ++i) {
      int p = edge(i, 0) - 1; // Parent (0-based)
      int c = edge(i, 1) - 1; // Child (0-based)
      
      isChild[c] = true;
      
      if (nodes[p].left == -1) nodes[p].left = c;
      else if (nodes[p].right == -1) nodes[p].right = c;
      else stop("Tree is not binary: Node %d has >2 children.", p+1);
    }
    
    // Find root
    rootIdx = -1;
    // Standard ape trees: internals start at nTips.
    // However, we check explicitly.
    for(int i=0; i < nTaxa; ++i) {
      // Tips are strictly 0..nTaxa-1? usually yes in ape.
      // If edge matrix allows reordering, we must be careful.
      // We assume 'ape' convention: 1..nTaxa are tips.
      nodes[i].isTip = true;
      nodes[i].nTips = 1;
    }
    
    for(int i=0; i <= maxId; ++i) {
      if (!isChild[i] && (nodes[i].left != -1 || nodes[i].right != -1)) {
        rootIdx = i;
        break;
      }
    }
    
    if (rootIdx == -1) stop("Could not identify root (or tree is a single node).");
    
    // Compute nTips recursively
    calcNTips(rootIdx);
  }
  
private:
  void calcNTips(int u) {
    if (nodes[u].isTip) {
      nodes[u].nTips = 1;
      return;
    }
    // Recurse
    if (nodes[u].left != -1) {
      calcNTips(nodes[u].left);
      nodes[u].nTips += nodes[nodes[u].left].nTips;
    }
    if (nodes[u].right != -1) {
      calcNTips(nodes[u].right);
      nodes[u].nTips += nodes[nodes[u].right].nTips;
    }
  }
};

// Global cache for the Pre-processed tree
// Key: A pointer or ID provided by R? 
// For safety/simplicity in this context, we will instantiate it per call or manage via an external pointer if speed is critical.
// Given "hundreds of times", passing an ExternalPtr is best.

// ============================================================================
// Fixed Tree Specific Generators
// ============================================================================

// We need a generator that returns draws specifically summing to 'target_m'.
// Reusing ValidDrawsCache logic but filtering for m.

struct FixedSplitGenerator {
  // We use the existing StateKey logic.
  // We will populate a FixedDraws struct but ONLY with draws where sum == target_m.
  
  template<typename KeyType, typename DrawsType>
  static void generate(const KeyType& leavesKey, int target_m, DrawsType& out) {
    std::vector<int> leaves(leavesKey.len);
    for(int i=0; i<leavesKey.len; ++i) leaves[i] = leavesKey.data[i];
    
    std::vector<int> drawn(leaves.size(), 0);
    out.count = 0;
    
    // We only want sum == target_m.
    // The "half" optimization from M&S (lex_compare) applies if left and right are
    // isomorphic, but here Left and Right children are distinct topological entities.
    // Therefore, (Drawn=A, Undrawn=B) is DISTINCT from (Drawn=B, Undrawn=A) 
    // unless the children themselves are structurally identical and we want to account for that symmetry.
    //
    // However, the user asked for "labelling the constituent taxa".
    // In a fixed tree, Left Child and Right Child are distinct slots.
    // So we must generate ALL splits of size target_m, without symmetry reduction 
    // based on drawn vs undrawn.
    
    rec(leaves, drawn, 0, target_m, 0, out);
  }
  
private:
  template<typename DrawsType>
  static void rec(const std::vector<int>& leaves, std::vector<int>& drawn, 
                  int idx, int target_m, int curSum, DrawsType& out) {
    
    if (curSum > target_m) return; // Pruning
    
    // Optimization: check if remaining items are enough to reach target
    int potential = 0;
    for(size_t i=idx; i<leaves.size(); ++i) potential += leaves[i];
    if (curSum + potential < target_m) return;
    
    if (idx == (int)leaves.size()) {
      if (curSum == target_m) {
        if (out.count >= 256) throw std::runtime_error("FixedDraws capacity exceeded in FixedSplit.");
        
        auto& dp = out.draws[out.count++];
        dp.drawn = decltype(dp.drawn)(drawn);
        
        // Calculate undrawn
        // We can't use the simple constructor that subtracts if we want to be generic,
        // but here we know the types.
        using KeyType = decltype(dp.drawn);
        KeyType total(leaves);
        dp.undrawn = KeyType(total, dp.drawn);
        dp.m = target_m;
      }
      return;
    }
    
    // Iterate
    int maxTake = std::min(leaves[idx], target_m - curSum);
    for (int k = 0; k <= maxTake; ++k) {
      drawn[idx] = k;
      rec(leaves, drawn, idx + 1, target_m, curSum + k, out);
    }
    drawn[idx] = 0;
  }
};


// ============================================================================
// Fixed Tree Solver Template
// ============================================================================

template<int nTokens>
class FixedTreeSolverT {
  using KeyType = StateKeyT<nTokens>;
  using HashType = StateKeyHashT<nTokens>;
  
  const FixedTreeStructure& tree;
  const Downpass& D;
  const TokenPairs& pairs;
  LnRootedCache& lnRooted;
  int presentBits;
  
  // Memoization: (NodeIndex, StateKey) -> Result
  // We use a combined hash.
  struct NodeKey {
    int u;
    KeyType key;
    bool operator==(const NodeKey& o) const { return u == o.u && key == o.key; }
  };
  
  struct NodeKeyHash {
    std::size_t operator()(const NodeKey& k) const {
      // Rotate the key hash and XOR with u
      size_t h = HashType{}(k.key);
      return h ^ (std::hash<int>{}(k.u) + 0x9e3779b9 + (h<<6) + (h>>2));
    }
  };
  
  // Caches
  std::unordered_map<NodeKey, FixedProbList, NodeKeyHash> b_cache;
  std::unordered_map<NodeKey, FixedProbList, NodeKeyHash> p_cache; // Storing all steps [0..max] roughly? No, steps are specific.
  
  // For P, we need (Node, Key, s, token). That's too big.
  // M&S Solver computed P one 's' at a time? 
  // M&S Solver stored P(s, key, token).
  // Here we have Node too.
  // Given "hundreds of times", we might just want to clear cache between trees, 
  // but keep it for the recursion within one tree.
  
  struct LogPKeyFixed {
    int u;
    KeyType leaves;
    int s;
    int token;
    bool operator==(const LogPKeyFixed& o) const {
      return u == o.u && s == o.s && token == o.token && leaves == o.leaves;
    }
  };
  
  struct LogPKeyFixedHash {
    std::size_t operator()(const LogPKeyFixed& k) const {
      size_t h = HashType{}(k.leaves);
      h = (h ^ k.u) * 3935559000370003845ULL;
      h = (h ^ k.s) * 4488902095908611103ULL;
      h = (h ^ k.token);
      return h;
    }
  };
  
  using LogPMapType = std::unordered_map<
    LogPKeyFixed, double, LogPKeyFixedHash, std::equal_to<LogPKeyFixed>,
    MallocPoolAllocator<std::pair<const LogPKeyFixed, double>>
  >;
  
  LogPMapType p_map;
  
  // Split Cache: (Key, TargetM) -> Splits
  // Since TargetM is implied by Node, we could just cache by (Key, TargetM).
  struct SplitKey {
    KeyType key;
    int m;
    bool operator==(const SplitKey& o) const { return m==o.m && key==o.key; }
  };
  struct SplitKeyHash {
    size_t operator()(const SplitKey& k) const { return HashType{}(k.key) ^ k.m; }
  };
  
  struct FixedDrawsT {
    struct { KeyType drawn; KeyType undrawn; int m; } draws[256];
    uint8_t count = 0;
  };
  
  std::unordered_map<SplitKey, FixedDrawsT, SplitKeyHash> split_cache;
  
  // Log Multinomial helper for logic
  // Returns log( n! / (prod n_i!) )
  double logMultinomial(const KeyType& k) {
    double res = R::lgammafn(k.sum() + 1.0);
    for(int i=0; i<k.len; ++i) {
      if(k.get(i) > 1) res -= R::lgammafn(k.get(i) + 1.0);
    }
    return res;
  }
  
  // We also need the number of ways to assign specific taxa to the splits.
  // In M&S (random tree), this was handled by LogRD.
  // In Fixed Tree, the logic is:
  // Node U has children Left (nL tips) and Right (nR tips).
  // We have a bag of states 'S' (sum = nL+nR).
  // We split 'S' into 'SL' (sum=nL) and 'SR' (sum=nR).
  // The number of ways to distribute the *taxa* to satisfy this split is:
  //   Ways = (Combinations of choosing SL states from S) * (Ways for Left) * (Ways for Right) ??
  //
  // Actually, standard Felsenstein pruning sums probabilities. 
  // If we want "Number of Labelings", it is:
  // Total Labelings * P(Data | Tree).
  //
  // P(Data | Tree) calculation:
  // At node u, for state x:
  //   L_u(x) = sum_{split S -> SL, SR} [ P(split) * (sum_{y,z} Trans(x->y)*L_v(y) * Trans(x->z)*L_w(z)) ]
  //
  // In our case, "P(split)" is determined by the permutations.
  // The probability that the left child (size nL) gets state-set SL given parent set S is:
  // Hypergeometric probability?
  //   Prob = choose(nL, SL_counts) * choose(nR, SR_counts) / choose(nL+nR, S_counts) ?? 
  //   No, that's for selecting balls.
  //
  // Let's stick to the definition:
  // We are calculating unnormalized conditional likelihoods (LogB / LogP).
  //
  // The M&S LogB function essentially calculates:
  // sum_{splits} ( Weight(split) * RecurseLeft * RecurseRight )
  //
  // For fixed tree, the Weight(split) is the number of ways to partition the set of tip-labels
  // such that the set of states mapping to the left child matches 'SL'.
  //
  // Given set S of states at U.
  // We choose nL specific taxa to go left.
  // The probability that a random choice of nL taxa has exactly state counts SL is:
  //   Prob = [ (Product_i choose(S_i, SL_i)) ] / choose(Total, nL)
  // This is the Multivariate Hypergeometric Distribution.
  //
  // So:
  // LogProbSplit = sum_i lchoose(S.get(i), SL.get(i)) - lchoose(S.sum(), nL)
  
  double computeSplitProb(const KeyType& parent, const KeyType& left) {
    double logNum = 0.0;
    for(int i=0; i<parent.len; ++i) {
      logNum += lchoose_log(parent.get(i), left.get(i));
    }
    double logDenom = lchoose_log(parent.sum(), left.sum());
    return logNum - logDenom;
  }
  
public:
  FixedTreeSolverT(const FixedTreeStructure& t, const Downpass& D_, const TokenPairs& p, LnRootedCache& lnr, int pb)
    : tree(t), D(D_), pairs(p), lnRooted(lnr), presentBits(pb) 
  {
    b_cache.reserve(4096);
    p_map.reserve(16384);
  }
  
  const FixedDrawsT& getSplits(const KeyType& key, int m) {
    SplitKey sk{key, m};
    auto it = split_cache.find(sk);
    if(it != split_cache.end()) return it->second;
    
    FixedDrawsT out;
    FixedSplitGenerator::generate(key, m, out);
    return split_cache.emplace(sk, out).first->second;
  }
  
  // LogB: Log Probability of observing subtree shapes given root state 'token0' and leaf-counts 'leaves'
  double LogB(int u, int token0, const KeyType& leaves) {
    const FixedTreeNode& node = tree.nodes[u];
    
    // 1. Base Case: Tip
    if (node.isTip) {
      // If we are at a tip, the 'leaves' key must have exactly count 1.
      assert(leaves.sum() == 1);
      return (leaves.get(token0) == 1) ? 0.0 : NEG_INF;
    }
    
    // 2. Check mask
    if ((token_mask(token0) & ~presentBits) != 0) return NEG_INF;
    
    // 3. Cache
    NodeKey nk{u, leaves};
    auto it = b_cache.find(nk);
    if (it == b_cache.end()) {
      it = b_cache.emplace(nk, FixedProbList{}).first;
    }
    double& slot = it->second[token0];
    if (std::isfinite(slot)) return slot;
    if (std::isnan(slot) == false && !std::isfinite(slot)) return slot;
    
    // 4. Recursion
    // Node u has children Left and Right.
    int leftId = node.left;
    int rightId = node.right;
    int nL = tree.nodes[leftId].nTips;
    // int nR = tree.nodes[rightId].nTips; // Implied
    
    const auto& splits = getSplits(leaves, nL);
    LSEAccumulator acc;
    
    for(int i=0; i<splits.count; ++i) {
      const auto& dp = splits.draws[i];
      // dp.drawn goes Left, dp.undrawn goes Right
      
      double splitProb = computeSplitProb(leaves, dp.drawn);
      
      LSEAccumulator inner;
      
      // Pairs Logic (same as M&S)
      // Iterate over valid transitions from token0 -> (a, b)
      for(const auto& pr : pairs.noStep[token0]) {
        double v = LogB(leftId, pr.a, dp.drawn) + LogB(rightId, pr.b, dp.undrawn);
        inner.add(v);
      }
      for(const auto& pr : pairs.yesStep[token0]) {
        double v = LogB(leftId, pr.a, dp.drawn) + LogB(rightId, pr.b, dp.undrawn);
        inner.add(v);
      }
      
      acc.add(splitProb + inner.result());
    }
    
    slot = acc.result();
    return slot;
  }
  
  double LogP(int s, int u, int token0, const KeyType& leaves) {
    const FixedTreeNode& node = tree.nodes[u];
    
    // 1. Base Case
    if (node.isTip) {
      if (leaves.sum() != 1) return NEG_INF;
      return (leaves.get(token0) == 1) ? ((s == 0) ? 0.0 : NEG_INF) : NEG_INF;
    }
    
    // 2. Cache
    LogPKeyFixed k{u, leaves, s, token0};
    auto it = p_map.find(k);
    if (it != p_map.end()) return it->second;
    
    // Denom check
    double denom = LogB(u, token0, leaves);
    if (!std::isfinite(denom)) {
      p_map.emplace(k, denom);
      return denom;
    }
    
    int leftId = node.left;
    int rightId = node.right;
    int nL = tree.nodes[leftId].nTips;
    
    const auto& splits = getSplits(leaves, nL);
    LSEAccumulator outerAcc;
    
    for(int i=0; i<splits.count; ++i) {
      const auto& dp = splits.draws[i];
      double splitProb = computeSplitProb(leaves, dp.drawn);
      
      LSEAccumulator noStepAcc;
      for (int r = 0; r <= s; ++r) {
        LSEAccumulator pairAcc;
        for (const auto& pr : pairs.noStep[token0]) {
          double t = log_prod_sum_4(
            LogP(r, leftId, pr.a, dp.drawn),
            LogB(leftId, pr.a, dp.drawn),
            LogP(s - r, rightId, pr.b, dp.undrawn),
            LogB(rightId, pr.b, dp.undrawn)
          );
          pairAcc.add(t);
        }
        noStepAcc.add(pairAcc.result());
      }
      
      LSEAccumulator yesStepAcc;
      if (s >= 1) {
        for (int r = 0; r <= s - 1; ++r) {
          LSEAccumulator pairAcc;
          for (const auto& pr : pairs.yesStep[token0]) {
            double t = log_prod_sum_4(
              LogP(r, leftId, pr.a, dp.drawn),
              LogB(leftId, pr.a, dp.drawn),
              LogP(s - r - 1, rightId, pr.b, dp.undrawn),
              LogB(rightId, pr.b, dp.undrawn)
            );
            pairAcc.add(t);
          }
          yesStepAcc.add(pairAcc.result());
        }
      }
      
      LSEAccumulator both;
      both.add(noStepAcc.result());
      both.add(yesStepAcc.result());
      
      outerAcc.add(splitProb + both.result());
    }
    
    double res = outerAcc.result() - denom;
    p_map.emplace(k, res);
    return res;
  }
  
  double run(int steps, const KeyType& rootLeaves) {
    LSEAccumulator acc;
    int root = tree.rootIdx;
    
    // We sum over all possible root states
    for(int t=0; t<D.nStates; ++t) {
      double b = LogB(root, t, rootLeaves);
      double p = LogP(steps, root, t, rootLeaves);
      if(std::isfinite(b) && std::isfinite(p)) {
        acc.add(b + p);
      }
    }
    
    // The result 'acc.result()' is Log(Probability).
    // To get the Number of Ways, we multiply by the Total Permutations.
    // Total Perms = N! / (n1! n2! ...)
    double logTotalPerms = logMultinomial(rootLeaves);
    return acc.result() + logTotalPerms;
  }
};


// ============================================================================
// EXPORTED FUNCTIONS
// ============================================================================

// Wrap the Tree structure in an XPtr for passing between R and C++
using FixedTreePtr = std::shared_ptr<FixedTreeStructure>;

//' Pre-process a tree for Fixed-Tree Maddison Slatkin
//' 
//' @param edge The edge matrix (1-based, Nedges x 2)
//' @param nTips Number of tips
//' @export
// [[Rcpp::export]]
SEXP FixedTree_Preprocess(IntegerMatrix edge, int nTips) {
  FixedTreePtr ptr = std::make_shared<FixedTreeStructure>(edge, nTips);
  return Rcpp::XPtr<FixedTreeStructure>(new FixedTreeStructure(edge, nTips), true);
}

//' Calculate Log-Number of labelings with score k on a Fixed Tree
//' 
//' @param treePtr External pointer to pre-processed tree
//' @param steps vector of steps
//' @param states vector of states (raw counts)
//' @export
// [[Rcpp::export]]
NumericVector FixedTree_Count(SEXP treePtr, IntegerVector steps, IntegerVector states) {
  
  Rcpp::XPtr<FixedTreeStructure> tree(treePtr);
  int nTaxa = tree->nTaxa;
  int len = states.size();
  
  int nTokens = (int)std::floor(std::log2((double)len)) + 1;
  int nLevels = nTokens;
  int nStates = (1 << nLevels) - 1;
  
  // Parse states
  std::vector<int> leavesVec(nStates, 0);
  for (int i = 0; i < std::min(len, nStates); ++i) {
    int v = states[i];
    if (IntegerVector::is_na(v)) v = 0;
    leavesVec[i] = v;
  }
  
  // Validation
  if (sum_int(leavesVec) != nTaxa) {
    stop("Sum of state counts does not match number of tips in tree.");
  }
  
  // Shared Structs
  std::shared_ptr<Downpass> Dptr;
  {
    if(DP_CACHE.find(nLevels) == DP_CACHE.end()) DP_CACHE[nLevels] = std::make_shared<Downpass>(nLevels);
    Dptr = DP_CACHE[nLevels];
  }
  
  int presentBits = 0;
  for (int t = 0; t < nStates; ++t) if (leavesVec[t] > 0) presentBits |= token_mask(t);
  
  long long tp_key = ((long long)nLevels << 20) | presentBits;
  std::shared_ptr<TokenPairs> Tpptr;
  {
    if(TP_CACHE.find(tp_key) == TP_CACHE.end()) TP_CACHE[tp_key] = std::make_shared<TokenPairs>(*Dptr, presentBits);
    Tpptr = TP_CACHE[tp_key];
  }
  
  std::shared_ptr<LnRootedCache> LNRptr;
  {
    if(LNROOT_CACHE.find(nTaxa) == LNROOT_CACHE.end()) LNROOT_CACHE[nTaxa] = std::make_shared<LnRootedCache>(nTaxa);
    LNRptr = LNROOT_CACHE[nTaxa];
  }
  
  int k = steps.size();
  NumericVector out(k);
  
  // Dispatch
  if (nTokens == 2) {
    FixedTreeSolverT<2> solver(*tree, *Dptr, *Tpptr, *LNRptr, presentBits);
    StateKeyT<2> root(leavesVec);
    for(int i=0; i<k; ++i) out[i] = solver.run(steps[i], root);
  } else if (nTokens == 3) {
    FixedTreeSolverT<3> solver(*tree, *Dptr, *Tpptr, *LNRptr, presentBits);
    StateKeyT<3> root(leavesVec);
    for(int i=0; i<k; ++i) out[i] = solver.run(steps[i], root);
  } else if (nTokens == 4) {
    FixedTreeSolverT<4> solver(*tree, *Dptr, *Tpptr, *LNRptr, presentBits);
    StateKeyT<4> root(leavesVec);
    for(int i=0; i<k; ++i) out[i] = solver.run(steps[i], root);
  } else {
    FixedTreeSolverT<5> solver(*tree, *Dptr, *Tpptr, *LNRptr, presentBits);
    StateKeyT<5> root(leavesVec);
    for(int i=0; i<k; ++i) out[i] = solver.run(steps[i], root);
  }
  
  return out;
}

//' Total Number of Unique Labelings (Log Scale)
//' 
//' @param states Integer vector of counts for each state
//' @export
// [[Rcpp::export]]
double LogMultinomial(IntegerVector states) {
  int n = 0;
  double denom = 0.0;
  for(int i=0; i<states.size(); ++i) {
    int k = states[i];
    if(k > 0) {
      n += k;
      denom += R::lgammafn(k + 1.0);
    }
  }
  return R::lgammafn(n + 1.0) - denom;
}
