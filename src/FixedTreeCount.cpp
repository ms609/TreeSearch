#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include <limits>

using namespace Rcpp;

// constant for log(0)
const double LOG_ZERO = -std::numeric_limits<double>::infinity();

// -----------------------------------------------------------------------------
// Helper: LogSumExp
// Computes log(exp(a) + exp(b)) safely
// -----------------------------------------------------------------------------
inline double logAdd(double logA, double logB) {
  if (logA == LOG_ZERO) return logB;
  if (logB == LOG_ZERO) return logA;
  if (logA > logB) {
    return logA + std::log1p(std::exp(logB - logA));
  }
  return logB + std::log1p(std::exp(logA - logB));
}

// -----------------------------------------------------------------------------
// Helper: Fitch Parsimony Operation
// Returns pair {state_mask, added_cost}
// -----------------------------------------------------------------------------
struct FitchResult {
  int state;
  int cost;
};

inline FitchResult fitchOp(int mask1, int mask2) {
  int intersection = mask1 & mask2;
  if (intersection > 0) {
    return {intersection, 0};
  } else {
    return {mask1 | mask2, 1};
  }
}

// -----------------------------------------------------------------------------
// Class to encapsulate the recursion and cache
// -----------------------------------------------------------------------------
class TreeCounter {
public:
  // Inputs
  const IntegerMatrix edge;
  const int nTip;
  const int maxSteps;
  const int nRows; // 2^nStates - 1 (though we allocate 2^n for easy indexing)
  
  // Computed Node Sizes (number of leaves under each node)
  std::vector<int> nodeSizes;
  
  // Memoization Cache
  // Key: Pair(NodeID, TokenVector)
  // Value: Flattened Matrix (rows=masks, cols=steps)
  std::map<std::pair<int, std::vector<int>>, std::vector<double>> cache;
  
  TreeCounter(IntegerMatrix edge_, int nTip_, int maxSteps_, int nStates_) 
    : edge(edge_), nTip(nTip_), maxSteps(maxSteps_), 
      nRows((1 << nStates_)) { // Allocate full power of 2 for direct indexing
    
    // Pre-calculate node sizes (number of tips descendant from node)
    // R nodes: 1..nTip are tips, nTip+1.. are internal
    int maxNode = nTip + edge.nrow(); // Approximation of max node index
    nodeSizes.resize(maxNode + 1, 0);
    
    // Initialize tips
    for(int i = 1; i <= nTip; ++i) nodeSizes[i] = 1;
    
    // Post-order traversal to sum sizes (assuming edges are standard phylo order)
    // Safest to do a quick recursive compute or iterate backwards if edges are post-order.
    // Here we use a recursive filler to be robust against edge ordering.
    computeNodeSizes(nTip + 1);
  }
  
  int computeNodeSizes(int node) {
    if (node <= nTip) return 1;
    if (nodeSizes[node] > 0) return nodeSizes[node]; // Already computed
    
    // Find children
    int size = 0;
    // Note: This linear scan is slow for massive trees, 
    // but standard R "edge" matrices allow O(N) if we assumed post-order.
    // For safety with generic inputs, we scan. 
    // (Optimization: Build an adjacency list in constructor if N is huge).
    for(int i = 0; i < edge.nrow(); ++i) {
      if (edge(i, 0) == node) {
        size += computeNodeSizes(edge(i, 1));
      }
    }
    nodeSizes[node] = size;
    return size;
  }
  
  // ---------------------------------------------------------------------------
  // Recursively generate splits of tokens
  // ---------------------------------------------------------------------------
  void generateSplits(size_t tokenTypeIdx, int currentLeftSize, int targetLeftSize,
                      const std::vector<int>& totalTokens,
                      std::vector<int>& currentSplit,
                      std::vector<std::vector<int>>& validSplits) {
    
    // Base case: processed all token types
    if (tokenTypeIdx == totalTokens.size()) {
      if (currentLeftSize == targetLeftSize) {
        validSplits.push_back(currentSplit);
      }
      return;
    }
    
    // Pruning: if we can't possibly fill the target, or already exceeded
    int remainingCapacity = targetLeftSize - currentLeftSize;
    if (remainingCapacity < 0) return;
    
    // Calculate max tokens of this type we can take
    int available = totalTokens[tokenTypeIdx];
    
    // Optimization: Don't iterate 0..available if we MUST take some to fill quota
    // Sum of all FUTURE tokens
    int futureTokens = 0;
    for(size_t i = tokenTypeIdx + 1; i < totalTokens.size(); ++i) futureTokens += totalTokens[i];
    
    int minToTake = std::max(0, remainingCapacity - futureTokens);
    int maxToTake = std::min(available, remainingCapacity);
    
    for (int count = minToTake; count <= maxToTake; ++count) {
      currentSplit[tokenTypeIdx] = count;
      generateSplits(tokenTypeIdx + 1, currentLeftSize + count, targetLeftSize,
                     totalTokens, currentSplit, validSplits);
    }
  }
  
  // ---------------------------------------------------------------------------
  // Core Recursive Function
  // ---------------------------------------------------------------------------
  std::vector<double> recurse(int node, const std::vector<int>& currentTokens) {
    
    // 1. Check Cache
    std::pair<int, std::vector<int>> key = {node, currentTokens};
    if (cache.count(key)) {
      return cache[key];
    }
    
    // Initialize Result Matrix (Log-Space 0 = -Inf)
    // Dimensions: nRows (masks) x (maxSteps + 1)
    // Flattened: index = mask + nRows * step
    std::vector<double> resMat(nRows * (maxSteps + 1), LOG_ZERO);
    
    // 2. Base Case: Leaf Node
    if (node <= nTip) {
      // Find which token is here (index where count is 1)
      int tokenIndex = -1;
      for(size_t i = 0; i < currentTokens.size(); ++i) {
        if (currentTokens[i] == 1) {
          tokenIndex = i;
          break;
        }
      }
      
      if (tokenIndex != -1) {
        // State mask is 2^(tokenIndex). Note: user code used 2^(idx-1), 
        // implies 0-based token index in power.
        // We will stick to 0-based index for tokens.
        // If tokenIndex is 0 (first token type), mask is 1.
        int stateMask = 1 << tokenIndex; 
        
        // At step 0, log probability is 0 (prob 1)
        // Index: row [stateMask], col [0]
        int idx = stateMask + nRows * 0; 
        resMat[idx] = 0.0; 
      }
      
      cache[key] = resMat;
      return resMat;
    }
    
    // 3. Recursive Step: Internal Node
    // Identify children
    int leftNode = -1, rightNode = -1;
    // Simple scan for children (assuming binary tree)
    for(int i = 0; i < edge.nrow(); ++i) {
      if (edge(i, 0) == node) {
        if (leftNode == -1) leftNode = edge(i, 1);
        else rightNode = edge(i, 1);
      }
    }
    
    int leftSize = nodeSizes[leftNode];
    
    // Generate Splits
    std::vector<std::vector<int>> splits;
    std::vector<int> buffer(currentTokens.size());
    generateSplits(0, 0, leftSize, currentTokens, buffer, splits);
    
    // Process Splits
    for (const auto& tokensL : splits) {
      std::vector<int> tokensR = currentTokens;
      for(size_t i = 0; i < tokensR.size(); ++i) tokensR[i] -= tokensL[i];
      
      std::vector<double> matL = recurse(leftNode, tokensL);
      std::vector<double> matR = recurse(rightNode, tokensR);
      
      // Convolve
      // Iterate over non-zero entries in L and R
      // Since we are in log space, "non-zero" means > LOG_ZERO
      
      // Optimization: Collect valid indices first to avoid O(N^2) loop over empties?
      // Given the matrix size is likely small (masks < 128, steps < 20),
      // dense iteration is acceptable and likely faster than sparse overhead 
      // unless sparse factor is huge.
      
      for (int stepL = 0; stepL <= maxSteps; ++stepL) {
        for (int maskL = 1; maskL < nRows; ++maskL) {
          double valL = matL[maskL + nRows * stepL];
          if (valL == LOG_ZERO) continue;
          
          for (int stepR = 0; stepR <= maxSteps - stepL; ++stepR) { // Pruning step sum
            for (int maskR = 1; maskR < nRows; ++maskR) {
              double valR = matR[maskR + nRows * stepR];
              if (valR == LOG_ZERO) continue;
              
              FitchResult fitch = fitchOp(maskL, maskR);
              int totalSteps = stepL + stepR + fitch.cost;
              
              if (totalSteps <= maxSteps) {
                int targetIdx = fitch.state + nRows * totalSteps;
                double combinedProb = valL + valR; // Log space mult
                resMat[targetIdx] = logAdd(resMat[targetIdx], combinedProb);
              }
            }
          }
        }
      }
    }
    
    cache[key] = resMat;
    return resMat;
  }
};

// -----------------------------------------------------------------------------
// Main Exported Function
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector FixedTreeCount(List tree, std::vector<int> tokens,
                             double steps = -1.0) {
  
  // 1. Setup
  IntegerMatrix edge = tree["edge"];
  int nTip = 0;
  
  // R "phylo" objects usually have a "tip.label" vector
  if (tree.containsElementNamed("tip.label")) {
    CharacterVector tips = tree["tip.label"];
    nTip = tips.size();
  } else {
    // Fallback if tip.label missing: max value in edge matrix
    // (This is rare for valid phylo objects)
    for(int i=0; i<edge.length(); ++i) if(edge[i] > nTip) nTip = edge[i];
    // Usually internal nodes are > nTip, so this heuristic is imperfect 
    // without standard phylo structure. Assuming standard inputs.
  }
  
  // Clean tokens (remove zeros and sort) - mimicking R logic
  std::vector<int> activeTokens;
  for(int t : tokens) {
    if (t > 0) activeTokens.push_back(t);
  }
  std::sort(activeTokens.begin(), activeTokens.end());
  
  int sumTokens = 0;
  for(int t : activeTokens) sumTokens += t;
  
  if (sumTokens != nTip) {
    stop("Number of leaves does not match total tokens.");
  }
  
  // Calculate Max Steps
  // R logic: min(steps, nTip - length(tokens) + 1)
  int calcMaxSteps = nTip - activeTokens.size() + 1;
  int limitSteps = (steps < 0) ? calcMaxSteps : (int)steps;
  int actualMaxSteps = std::min(limitSteps, calcMaxSteps);
  
  // 2. Instantiate Logic
  TreeCounter solver(edge, nTip, actualMaxSteps, activeTokens.size());
  
  // 3. Run (Root is typically nTip + 1)
  int rootNode = nTip + 1;
  std::vector<double> finalMat = solver.recurse(rootNode, activeTokens);
  
  // 4. Summarize Results
  // The result is colSums of finalMat. 
  // finalMat is flattened [mask + nRows * step]
  
  NumericVector results(actualMaxSteps + 1);
  
  for (int s = 0; s <= actualMaxSteps; ++s) {
    double stepTotal = LOG_ZERO;
    for (int m = 1; m < solver.nRows; ++m) {
      stepTotal = logAdd(stepTotal, finalMat[m + solver.nRows * s]);
    }
    results[s] = stepTotal;
  }
  
  // Assign names "0", "1", ...
  CharacterVector names(actualMaxSteps + 1);
  for(int i=0; i<=actualMaxSteps; ++i) names[i] = std::to_string(i);
  results.attr("names") = names;
  
  return results;
}
