#include <Rcpp.h>
#include <bitset>
#include <vector>
#include <algorithm>
#include <cmath>

constexpr int max_token_bits = 64;

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
  std::vector<int> counts(max_val, 0);
  int nInapp = 0;
  
  for (int xi : x) {
    if (xi == 0) {
      ++nInapp;
    } else if (xi > 0 && xi <= max_val) {
      counts[xi - 1]++;
    } else {
      Rcpp::stop("Input out of range: " + std::to_string(xi));
    }
  }
  
  int regions = std::max(1, nInapp - 1);
  int steps = 0;
  
  // Step 2: Proceed only if more than one state observed
  int nPositiveStates = std::count_if(counts.begin(), counts.end(),
                                      [](int c) {return c > 0; });
  if (nPositiveStates <= 1) {
    return steps + std::max(0, std::min(std::accumulate(counts.begin(), counts.end(), 0), regions) - 1);
  }
  
  int nState = std::floor(std::log2(counts.size())) + 1;
  if (nState > max_token_bits) {
    Rcpp::stop("maximum_length(): too many states (nState > 64) for uint64_t representation.");
  }
  
  const int nToken = (1 << nState) - 1;
  counts.resize(nToken, 0);
  
  // Step 3: Generate tokens
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
      if (counts[i] > 0 && active[i]) {
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
        if (tokenSums[i] != amb || counts[i] == 0 || !active[i]) {
          continue;
        }
        // Compute options: tokens that are active, counts > 0,
        // and non-intersecting with token i
        std::vector<int> optionIndices;
        std::vector<int> candidateIndices;
        for (int j = 0; j < nToken; ++j) {
          if (active[j] &&
              counts[j] > 0 &&
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
          int maxCount = counts[chosen];
          for (size_t idx = 1; idx < candidateIndices.size(); ++idx) {
            if (counts[candidateIndices[idx]] > maxCount) {
              chosen = candidateIndices[idx];
              maxCount = counts[chosen];
            }
          }
          
          // Update counts: decrement i and chosen, increment merged token
          counts[i]--;
          counts[chosen]--;
          const int product = Merge(i, chosen);
          counts[product]++;
          ++steps;
          
          escape = true;   // Mark a successful merge
          break;           // Break out of inner i-loop
        }
      }
      
      if (escape) break;  // Break out of unionSize loop if a merge happened
    }
  }
  
  int remaining = std::accumulate(counts.begin(), counts.end(), 0);
  return steps + std::max(0, std::min(remaining, regions) - 1);
}
