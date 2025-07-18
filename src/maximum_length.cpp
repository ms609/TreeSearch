#include <Rcpp.h>
#include <bitset>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>

constexpr int max_token_bits = 64;

struct ActiveTokenLog {
  const int totally_ambiguous;
  std::unordered_set<int> active_set;
  
  // Initialize with all tokens active except ambiguity token
  ActiveTokenLog(int nToken) : totally_ambiguous(nToken - 1) {}
  
  bool isActive(int i) const {
    if (i == totally_ambiguous) return false;
    return active_set.find(i) != active_set.end();
  }
  
  void deactivate(int i) {
    active_set.erase(i);
  }
  
  void activate(int i) {
    if (i != totally_ambiguous) active_set.insert(i);
  }
  
  bool empty() const {
    return active_set.empty();
  }
  
  template <typename Func>
  void forEachActive(Func f) const {
    // Iterate over each element in the set *now* without updating.
    std::vector<int> snapshot(active_set.begin(), active_set.end());
    for (int i : snapshot) {
      f(i);
    }
  }
};

struct Counts {
  ActiveTokenLog& token_log;
  std::unordered_map<int, int> data;
  
  Counts(ActiveTokenLog& log) : token_log(log) {}
  
  // Get count at index
  int get(int index) const {
    auto it = data.find(index);
    return (it != data.end()) ? it->second : 0;
  }
  
  // Increment by one
  void increment(int index) {
    token_log.activate(index);
    ++data[index];
  }
  
  // Decrement by one, with safety check
  void decrement(int index) {
    auto it = data.find(index);
    if (it != data.end() && it->second > 0) {
      --(it->second);
      if (it->second == 0) {
        data.erase(it);  // optional: reduce memory
      }
    }
  }
  
  // Set a count explicitly
  void set(int index, int value) {
    if (value <= 0) {
      data.erase(index);
    } else {
      data[index] = value;
    }
  }
  
  // Total remaining tokens
  int total() const {
    int sum = 0;
    for (const auto& kv : data) sum += kv.second;
    return sum;
  }
  
  int num_positive() const {
    return static_cast<int>(data.size());
  }
  
  // Return all active indices
  std::vector<int> indices() const {
    std::vector<int> result;
    result.reserve(data.size());
    for (const auto& kv : data) result.push_back(kv.first);
    return result;
  }
  
  // Check if index exists with positive count
  bool has(int index) const {
    return get(index) > 0;
  }
};

struct Token {
  static inline uint64_t token(int i) {
    return static_cast<uint64_t>(i + 1);
  }
  
  static inline int sum(int i) {
    return __builtin_popcountll(token(i));
  }
  
  static inline bool intersect(int i, int j) {
    return token(i) & token(j);
  }
  
  static inline int union_size(int i, int j) {
    return __builtin_popcountll(token(i) | token(j));
  }
  
  static inline int merge(int i, int j) {
    return static_cast<int>(token(i) | token(j)) - 1;
  }
};

// Draft by ChatGPT
// [[Rcpp::export]]
int maximum_length(const Rcpp::IntegerVector& x) {
  
  // Step 1: Preprocessing
  if (x.length() < 1) {
    Rcpp::stop("Input vector must not be empty.");
  } else if (x.length() == 1) {
    return 0;
  }
  int max_val = *std::max_element(x.begin(), x.end());
  if (max_val < 0) {
    Rcpp::stop("Input vector must contain positive integers.");
  } else if (max_val == 0) {
    return 0; // Inapplicables only
  }
  
  const int nState = std::floor(std::log2(max_val)) + 1;
  if (nState > max_token_bits) {
    Rcpp::stop("maximum_length(): too many states (nState > 64) for uint64_t representation.");
  }
  
  const int nToken = (1 << nState) - 1;
  ActiveTokenLog token(nToken);
  Counts counts(token);
  int nInapp = 0;
  
  for (int xi : x) {
    if (xi == 0) {
      ++nInapp;
    } else if (xi > 0 && xi <= max_val) {
      counts.increment(xi - 1);
    } else {
      Rcpp::stop("Input out of range: " + std::to_string(xi));
    }
  }
  
  int regions = std::max(1, nInapp - 1);
  int steps = 0;
  
  // Step 2: Proceed only if more than one state observed
  if (counts.num_positive() <= 1) {
    return steps + std::max(0, std::min(counts.total(), regions) - 1);
  }
  
  
  int loopCount = 0;
  bool escape = false;
  
  while (true) {
    int amb = -1;
    token.forEachActive([&](int i) {
      if (counts.get(i) > 0) {
        amb = std::max(amb, Token::sum(i));
      }
    });
    
    // Break if no ambiguous tokens remain
    if (amb < 1) break;
    
    ++loopCount;
    if (loopCount > 10000) {
      Rcpp::stop("maximum_length() failed. Please report this bug.");
    }
    
    escape = false;
    
    // Outer for loop: descending union sizes from nState down to amb + 1
    // Objective: pair ...+++ with +++... to yield ++++++
    //            before considering ++.... for ++.+++
    for (int unionSize = nState; unionSize >= amb + 1; --unionSize) {
      token.forEachActive([&](int i) {
        assert(token.isActive(i));
        if (counts.get(i) == 0 || Token::sum(i) != amb) {
          return;
        }
        
        // Compute options: tokens that are active, counts > 0,
        // and non-intersecting with token i
        std::vector<int> optionIndices;
        std::vector<int> candidateIndices;
        token.forEachActive([&](int j) {
          if (j != i && counts.get(j) > 0 && !Token::intersect(i, j)) {
            assert(token.isActive(j)); // Validate snapshot active
            // We have an option for the future, though other options
            // potentially yield a larger union.
            optionIndices.push_back(j);
            if (Token::union_size(i, j) == unionSize) {
              // No options yield a larger union, so consider this candidate now
              candidateIndices.push_back(j);
            }
          }
        });
        
        if (optionIndices.empty()) {
          // No options: disable token i
          token.deactivate(i);
        } else if (!candidateIndices.empty()) {
          // Find candidate with maximum count
          int chosen = candidateIndices[0];
          int maxCount = counts.get(chosen);
          for (size_t idx = 1; idx < candidateIndices.size(); ++idx) {
            if (counts.get(candidateIndices[idx]) > maxCount) {
              chosen = candidateIndices[idx];
              maxCount = counts.get(chosen);
            }
          }
          
          // Update counts: decrement i and chosen, increment merged token
          counts.decrement(i);
          counts.decrement(chosen);
          const int product = Token::merge(i, chosen);
          counts.increment(product);
          ++steps;
          
          escape = true;   // Mark a successful merge
          return;           // Break out of inner i-loop
        }
      });
      
      if (escape) break;  // Break out of unionSize loop if a merge happened
    }
  }
  
  return steps + std::max(0, std::min(counts.total(), regions) - 1);
}
