#include <Rcpp.h>
#include <bitset>
#include <vector>
#include <algorithm>
#include <cmath>

constexpr int max_token_bits = 64;

struct Counts {
  std::unordered_map<int, int> data;
  
  // Get count at index
  int get(int index) const {
    auto it = data.find(index);
    return (it != data.end()) ? it->second : 0;
  }
  
  // Increment by one
  void increment(int index) {
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
  Counts counts;
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
  
  const int nState = std::floor(std::log2(max_val)) + 1;
  if (nState > max_token_bits) {
    Rcpp::stop("maximum_length(): too many states (nState > 64) for uint64_t representation.");
  }
  
  // Step 3: Generate tokens
  const int nToken = (1 << nState) - 1;
  std::vector<uint64_t> tokens(nToken);
  std::vector<int> tokenSums(nToken);
  
  for (int i = 0; i < nToken; ++i) {
    uint64_t token = static_cast<uint64_t>(i + 1);  // tokens 1..n map directly to bits
    tokens[i] = token;
    tokenSums[i] = __builtin_popcountll(token);
  }
  
  std::vector<bool> active(nToken, true);
  active[nToken - 1] = false;  // Final token (all bits) is ambiguity
  
  // Step 4: Precompute intersections and unions
  std::vector<std::vector<bool>> nonIntersect(nToken, std::vector<bool>(nToken));
  std::vector<std::vector<int>> unions(nToken, std::vector<int>(nToken));
  for (int i = 0; i < nToken; ++i) {
    for (int j = 0; j < nToken; ++j) {
      uint64_t a = tokens[i];
      uint64_t b = tokens[j];
      nonIntersect[i][j] = !(a & b);
      unions[i][j] = __builtin_popcountll(a | b);
    }
  }
  
  // Step 5: Define Merge()
  auto Merge = [&](const int a, const int b) -> int {
    uint64_t merged = tokens[a] | tokens[b];
    // Bit trick: token i has value i + 1
    int merged_index = static_cast<int>(merged) - 1;
    if (merged_index < 0 || merged_index >= nToken) {
      Rcpp::stop("Internal error: merged token index out of bounds.");
    }
    return merged_index;
  };
  
  int loopCount = 0;
  bool escape = false;
  
  while (true) {
    int amb = -1;
    for (int i = 0; i < nToken; ++i) {
      if (counts.get(i) > 0 && active[i]) {
        amb = std::max(amb, tokenSums[i]);
      }
    }
    
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
      for (int i = 0; i < nToken; ++i) {
        if (tokenSums[i] != amb || counts.get(i) == 0 || !active[i]) {
          continue;
        }
        // Compute options: tokens that are active, counts > 0,
        // and non-intersecting with token i
        std::vector<int> optionIndices;
        std::vector<int> candidateIndices;
        for (int j = 0; j < nToken; ++j) {
          if (active[j] &&
              counts.get(j) > 0 &&
              nonIntersect[i][j]) {
            // We have an option for the future, though other options
            // potentially yield a larger union.
            optionIndices.push_back(j);
            if (unions[i][j] == unionSize) {
              // No options yield a larger union, so consider this candidate now
              candidateIndices.push_back(j);
            }
          }
        }
        
        if (optionIndices.empty()) {
          // No options: disable token i
          active[i] = false;
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
          const int product = Merge(i, chosen);
          counts.increment(product);
          ++steps;
          
          escape = true;   // Mark a successful merge
          break;           // Break out of inner i-loop
        }
      }
      
      if (escape) break;  // Break out of unionSize loop if a merge happened
    }
  }
  
  return steps + std::max(0, std::min(counts.total(), regions) - 1);
}
