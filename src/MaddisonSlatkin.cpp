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
  double maxv;
  double acc;
  bool empty;

  LSEAccumulator() : maxv(NEG_INF), acc(0.0), empty(true) {}

  inline void add(double v) {
    if (!(v > NEG_INF)) return;

    if (empty) {
      maxv = v;
      acc = 1.0;
      empty = false;
    } else {
      if (v > maxv) {
        acc = acc * std::exp(maxv - v);
        maxv = v;
      }
      acc += std::exp(v - maxv);
    }
  }

  inline double result() const {
    if (empty) return NEG_INF;
    return maxv + std::log(acc);
  }
};

static inline double log_prod_sum_4(double v1, double v2, double v3, double v4) {
  if (!(v1 > NEG_INF) || !(v2 > NEG_INF) ||
      !(v3 > NEG_INF) || !(v4 > NEG_INF)) return NEG_INF;
  return v1 + v2 + v3 + v4;
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
    double lc = 0.0;
    for (int i = 0; i < leaves.len; ++i) {
      lc += lchoose_log(leaves.data[i], drawn.data[i]);
    }
    double val = bal + lnRooted(m) + lnRooted(n - m) - lnRooted(n) + lc;
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

  // Vectorized logP cache: key = (token, leaves) → vector over all step counts.
  // Eliminates the s dimension from the cache key, giving ~s_max× fewer entries.
  struct LogPVecKey {
    KeyType leaves;
    int token;
    LogPVecKey(int t_, const KeyType& l_) : leaves(l_), token(t_) {}
    bool operator==(const LogPVecKey& o) const {
      return token == o.token && leaves == o.leaves;
    }
  };
  struct LogPVecKeyHash {
    std::size_t operator()(const LogPVecKey& k) const noexcept {
      uint64_t h = HashType{}(k.leaves);
      h = (h ^ (uint64_t)k.token) * 4488902095908611103ULL;
      return (std::size_t)h;
    }
  };
  std::unordered_map<LogPVecKey, std::vector<double>,
                     LogPVecKeyHash> logPVec_cache;

  // Global s_max for this solver instance (set at first run() call)
  int s_max_global = 0;
    
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
      double lc = 0.0;
      for (int i = 0; i < leaves.len; ++i) {
        lc += lchoose_log(leaves.data[i], drawn.data[i]);
      }
      double val = bal + lnRooted(m) + lnRooted(n - m) - lnRooted(n) + lc;
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
      if ((slot > NEG_INF)) return slot;
      if (std::isnan(slot) == false && !(slot > NEG_INF)) return slot;
      
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
          innerAcc.add(val);
        }
        for (const auto& pr : pairs.yesStep[token0]) {
          double val = LogB(pr.a, drawn) + LogB(pr.b, undrawn);
          innerAcc.add(val);
        }
        double innerSum = innerAcc.result();
        
        double acc = balancedCorrection + computeLogRD(drawn, leaves) + innerSum;
        outerAcc.add(acc);
      }
      slot = outerAcc.result();
      return slot;
    }
    
    // Log-space convolution: C[s] = LogSumExp_r( A[r] + B[s-r] )
    // Only iterates over the active (finite) range of each operand.
    static std::vector<double> logconv(
        const std::vector<double>& A, int a_lo, int a_hi,
        const std::vector<double>& B, int b_lo, int b_hi,
        int c_size) {
      std::vector<double> C(c_size, NEG_INF);
      for (int r = a_lo; r <= a_hi; ++r) {
        const double va = A[r];
        if (!(va > NEG_INF)) continue;
        for (int q = b_lo; q <= b_hi; ++q) {
          const double vb = B[q];
          if (!(vb > NEG_INF)) continue;
          const int s = r + q;
          if (s >= c_size) break;
          const double v = va + vb;
          // Inline log-sum-exp update: C[s] = log(exp(C[s]) + exp(v))
          if (!(C[s] > NEG_INF)) {
            C[s] = v;
          } else {
            double mx = std::max(C[s], v);
            C[s] = mx + std::log1p(std::exp(std::min(C[s], v) - mx));
          }
        }
      }
      return C;
    }

    // Returns LogPVec for (leaves, token0): a vector of length (s_max+1)
    // where entry [s] = log P(subtree on `leaves` has exactly s steps
    //                       | root token = token0).
    // The vector is indexed from 0; entries before the min possible step
    // count are NEG_INF.
    const std::vector<double>& LogPVec(const KeyType& leaves, int token0) {
      static const std::vector<double> empty_neg_inf;
      const int n = leaves.sum();
      const int s_max = s_max_global;
      const int c_size = s_max + 1;

      // --- Base case n == 1 ---
      if (n == 1) {
        LogPVecKey key(token0, leaves);
        auto it = logPVec_cache.find(key);
        if (it != logPVec_cache.end()) return it->second;
        std::vector<double> v(c_size, NEG_INF);
        if (leaves.get(token0) == 1) v[0] = 0.0;
        auto ins = logPVec_cache.emplace(std::move(key), std::move(v));
        return ins.first->second;
      }

      // --- Base case n == 2 ---
      if (n == 2) {
        LogPVecKey key(token0, leaves);
        auto it = logPVec_cache.find(key);
        if (it != logPVec_cache.end()) return it->second;
        std::vector<double> v(c_size, NEG_INF);
        int twice = -1, a = -1, b = -1;
        for (int i = 0; i < leaves.len; ++i) {
          int count = leaves.get(i);
          if (count == 2) { twice = i; break; }
          else if (count == 1) { if (a < 0) a = i; else b = i; }
        }
        int needed;
        if (twice >= 0) {
          needed = 0;
        } else {
          needed = D.step_at(a, b) ? 1 : 0;
        }
        if (needed <= s_max) v[needed] = 0.0;
        auto ins = logPVec_cache.emplace(std::move(key), std::move(v));
        return ins.first->second;
      }

      // --- Cache lookup ---
      LogPVecKey key(token0, leaves);
      {
        auto it = logPVec_cache.find(key);
        if (it != logPVec_cache.end()) return it->second;
      }

      double denom = LogB(token0, leaves);
      if (!(denom > NEG_INF)) {
        std::vector<double> v(c_size, NEG_INF);
        auto ins = logPVec_cache.emplace(std::move(key), std::move(v));
        return ins.first->second;
      }

      // Accumulator: outer[s] = LogSumExp over all draw partitions
      std::vector<double> outerVec(c_size, NEG_INF);

      const auto& drawpairs = getValidDraws(leaves);

      for (uint8_t i = 0; i < drawpairs.count; ++i) {
        const auto& dp = drawpairs.draws[i];
        const KeyType& drawn   = dp.drawn;
        const KeyType& undrawn = dp.undrawn;
        const int m = dp.m;

        double rdCorr = computeLogRD(drawn, leaves) +
                        (((m + m == n) && !(drawn == undrawn)) ? std::log(2.0) : 0.0);

        // Active ranges: min steps = 0 (all tips same state possible),
        // max steps = n_subtree - 1.
        const int d_lo = 0;
        const int d_hi = std::min(drawn.sum() - 1,   s_max);
        const int u_lo = 0;
        const int u_hi = std::min(undrawn.sum() - 1, s_max);

        // Pre-fetch FixedProbList* for drawn/undrawn once per draw pair.
        // Reduces logB_cache::find from 2*(N_no+N_yes) to 2 per draw pair.
        // Prime the cache with the first available pair so logB_cache[drawn]
        // exists before we take the pointer (n<=2 base cases don't use the cache,
        // so fpl_d/fpl_u remain null and LogB() handles those directly).
        const FixedProbList* fpl_d = nullptr;
        const FixedProbList* fpl_u = nullptr;
        {
          const auto& all_no  = pairs.noStep[token0];
          const auto& all_yes = pairs.yesStep[token0];
          if (!all_no.empty()) {
            LogPVec(drawn,   all_no[0].a);
            LogPVec(undrawn, all_no[0].b);
          } else if (!all_yes.empty()) {
            LogPVec(drawn,   all_yes[0].a);
            LogPVec(undrawn, all_yes[0].b);
          }
          { auto it = logB_cache.find(drawn);   if (it != logB_cache.end()) fpl_d = &it->second; }
          { auto it = logB_cache.find(undrawn); if (it != logB_cache.end()) fpl_u = &it->second; }
        }

        // noStep pairs: convolve LogPVec(drawn,a) with LogPVec(undrawn,b)
        std::vector<double> noStepVec(c_size, NEG_INF);
        for (const auto& pr : pairs.noStep[token0]) {
          const std::vector<double>& A = LogPVec(drawn,   pr.a);
          const std::vector<double>& B = LogPVec(undrawn, pr.b);
          // After LogPVec(drawn,pr.a) returns, fpl_d[pr.a] is non-NaN (LogB fills it).
          double logBa = fpl_d ? (*fpl_d)[pr.a] : LogB(pr.a, drawn);
          double logBb = fpl_u ? (*fpl_u)[pr.b] : LogB(pr.b, undrawn);
          if (!(logBa > NEG_INF) || !(logBb > NEG_INF)) continue;
          double pairScale = logBa + logBb;
          // C = conv(A, B), then add pairScale; merge into noStepVec
          auto C = logconv(A, d_lo, d_hi, B, u_lo, u_hi, c_size);
          for (int s = 0; s < c_size; ++s) {
            if (!(C[s] > NEG_INF)) continue;
            double v = C[s] + pairScale;
            if (!(noStepVec[s] > NEG_INF)) {
              noStepVec[s] = v;
            } else {
              double mx = std::max(noStepVec[s], v);
              noStepVec[s] = mx + std::log1p(std::exp(std::min(noStepVec[s], v) - mx));
            }
          }
        }

        // yesStep pairs: same convolution but shifted by 1 (one extra step consumed)
        std::vector<double> yesStepVec(c_size, NEG_INF);
        if (!pairs.yesStep[token0].empty()) {
          for (const auto& pr : pairs.yesStep[token0]) {
            const std::vector<double>& A = LogPVec(drawn,   pr.a);
            const std::vector<double>& B = LogPVec(undrawn, pr.b);
            double logBa = fpl_d ? (*fpl_d)[pr.a] : LogB(pr.a, drawn);
            double logBb = fpl_u ? (*fpl_u)[pr.b] : LogB(pr.b, undrawn);
            if (!(logBa > NEG_INF) || !(logBb > NEG_INF)) continue;
            double pairScale = logBa + logBb;
            auto C = logconv(A, d_lo, d_hi, B, u_lo, u_hi, c_size);
            // Shift by 1: C_shifted[s] = C[s-1]
            for (int s = c_size - 1; s >= 1; --s) {
              if (!(C[s - 1] > NEG_INF)) continue;
              double v = C[s - 1] + pairScale;
              if (!(yesStepVec[s] > NEG_INF)) {
                yesStepVec[s] = v;
              } else {
                double mx = std::max(yesStepVec[s], v);
                yesStepVec[s] = mx + std::log1p(std::exp(std::min(yesStepVec[s], v) - mx));
              }
            }
          }
        }

        // Combine noStep and yesStep, add rdCorr, merge into outerVec
        for (int s = 0; s < c_size; ++s) {
          double combined;
          bool ns = (noStepVec[s] > NEG_INF);
          bool ys = (yesStepVec[s] > NEG_INF);
          if (!ns && !ys) continue;
          else if (!ys) combined = noStepVec[s];
          else if (!ns) combined = yesStepVec[s];
          else {
            double mx = std::max(noStepVec[s], yesStepVec[s]);
            combined = mx + std::log1p(std::exp(std::min(noStepVec[s], yesStepVec[s]) - mx));
          }
          double v = rdCorr + combined;
          if (!(outerVec[s] > NEG_INF)) {
            outerVec[s] = v;
          } else {
            double mx = std::max(outerVec[s], v);
            outerVec[s] = mx + std::log1p(std::exp(std::min(outerVec[s], v) - mx));
          }
        }
      }

      // Subtract denom (LogB normaliser)
      for (int s = 0; s < c_size; ++s) {
        if ((outerVec[s] > NEG_INF)) outerVec[s] -= denom;
      }

      LogPVecKey key2(token0, leaves);
      auto ins = logPVec_cache.emplace(std::move(key2), std::move(outerVec));
      return ins.first->second;
    }
    
public:
  SolverT(const Downpass& D_, const TokenPairs& p, int presentBits_,
          LnRootedCache& lnr)
    : D(D_), pairs(p), presentBits(presentBits_), lnRooted(lnr) {

    logB_cache.reserve(16384);
    logB_cache.max_load_factor(0.5f);

    logPVec_cache.reserve(8192);
    logPVec_cache.max_load_factor(0.5f);

    logRD_cache.reserve(1024);
    validDraws_cache.reserve(256);
  }

  // Run over a vector of step counts.  All counts share the same solver
  // instance so the logPVec_cache is populated once and reused.
  void runAll(const std::vector<int>& steps_vec, const KeyType& states,
              double* out) {
    if (steps_vec.empty()) return;

    // Set global s_max so base-case vectors are correctly sized
    int s_max = *std::max_element(steps_vec.begin(), steps_vec.end());
    s_max_global = s_max;
    const int c_size = s_max + 1;

    // Compute LogB and LogPVec for each root token
    // LogPVec will recursively fill the cache for all sub-configurations
    std::vector<double> rootLogB(D.nStates);
    for (int token0 = 0; token0 < D.nStates; ++token0) {
      rootLogB[token0] = LogB(token0, states);
    }

    // Combine: log P(s steps | character) = LogSumExp_token(LogB(token) + LogPVec(token)[s])
    std::vector<double> combined(c_size, NEG_INF);
    for (int token0 = 0; token0 < D.nStates; ++token0) {
      const double lb = rootLogB[token0];
      if (!(lb > NEG_INF)) continue;
      const std::vector<double>& pv = LogPVec(states, token0);
      for (int s = 0; s < c_size; ++s) {
        if (!(pv[s] > NEG_INF)) continue;
        double v = lb + pv[s];
        if (!(combined[s] > NEG_INF)) {
          combined[s] = v;
        } else {
          double mx = std::max(combined[s], v);
          combined[s] = mx + std::log1p(std::exp(std::min(combined[s], v) - mx));
        }
      }
    }

    for (int i = 0; i < (int)steps_vec.size(); ++i) {
      int s = steps_vec[i];
      out[i] = (s >= 0 && s <= s_max) ? combined[s] : NEG_INF;
    }
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
    if ((slot > NEG_INF)) return slot;
    if (std::isnan(slot) == false && !(slot > NEG_INF)) return slot;
    
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
        innerAcc.add(val);
      }
      for (const auto& pr : pairs.yesStep[token0]) {
        double val = LogB(pr.a, drawn) + LogB(pr.b, undrawn);
        innerAcc.add(val);
      }
      double innerSum = innerAcc.result();
      
      double acc = balancedCorrection + logRD.compute(drawn, leaves) + innerSum;
      outerAcc.add(acc);
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
    if (!(denom > NEG_INF)) {
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
          pairAcc.add(t);
        }
        noStepAcc.add(pairAcc.result());
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
            pairAcc.add(t);
          }
          yesStepAcc.add(pairAcc.result());
        }
        yesStepSum = yesStepAcc.result();
      }
      
      LSEAccumulator bothAcc;
      bothAcc.add(noStepSum);
      bothAcc.add(yesStepSum);
      double combined = bothAcc.result();
      
      double inner = logRD.compute(drawn, leaves) + sizeCorrection + combined;
      outerAcc.add(inner);
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
      if (!(p > NEG_INF) || !(p > NEG_INF)) {
        val = NEG_INF;
      } else {
        val = b + p;
      }
      acc.add(val);
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
  NumericVector out(k, NA_REAL);

  // Collect non-NA step counts; map positions
  std::vector<int> valid_steps;
  std::vector<int> valid_idx;
  valid_steps.reserve(k);
  valid_idx.reserve(k);
  for (int i = 0; i < k; ++i) {
    if (!IntegerVector::is_na(steps[i])) {
      valid_steps.push_back(steps[i]);
      valid_idx.push_back(i);
    }
  }

  if (valid_steps.empty()) return out;

  std::vector<double> results(valid_steps.size());

  // DISPATCH based on nTokens — call runAll() once for all step counts
  if (nTokens == 2) {
    StateKeyT<2> rootKey(leavesVec);
    SolverT<2> solver(D, pairs, presentBits, lnRooted);
    solver.runAll(valid_steps, rootKey, results.data());
  } else if (nTokens == 3) {
    StateKeyT<3> rootKey(leavesVec);
    SolverT<3> solver(D, pairs, presentBits, lnRooted);
    solver.runAll(valid_steps, rootKey, results.data());
  } else if (nTokens == 4) {
    StateKeyT<4> rootKey(leavesVec);
    SolverT<4> solver(D, pairs, presentBits, lnRooted);
    solver.runAll(valid_steps, rootKey, results.data());
  } else { // nTokens == 5
    StateKeyT<5> rootKey(leavesVec);
    SolverT<5> solver(D, pairs, presentBits, lnRooted);
    solver.runAll(valid_steps, rootKey, results.data());
  }

  for (int i = 0; i < (int)valid_idx.size(); ++i) {
    out[valid_idx[i]] = results[i];
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
