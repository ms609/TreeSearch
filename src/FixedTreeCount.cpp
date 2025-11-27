#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <limits>
#include <cstring>

#include <TreeTools/renumber_tree.h> /* for preorder_edges_and_nodes */

using namespace Rcpp;

// constant for log(0)
const double LOG_ZERO = -std::numeric_limits<double>::infinity();

// -----------------------------------------------------------------------------
// Helper: LogSumExp
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
// Efficient cache key using a single 64-bit hash for small token vectors
// For larger vectors, store them but use a precomputed hash
// -----------------------------------------------------------------------------
struct CacheKey {
  int node;
  size_t tokenHash; // Precomputed hash of tokens
  const std::vector<int>* tokens; // Pointer to avoid copying
  
  CacheKey(int n, const std::vector<int>& t, size_t h) 
    : node(n), tokenHash(h), tokens(&t) {}
  
  bool operator==(const CacheKey& other) const {
    return node == other.node && 
      tokenHash == other.tokenHash &&
      *tokens == *other.tokens;
  }
};

struct CacheKeyHash {
  std::size_t operator()(const CacheKey& k) const {
    // Combine node and precomputed token hash
    size_t h = std::hash<int>()(k.node);
    h ^= k.tokenHash + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
  }
};

// -----------------------------------------------------------------------------
// Helper to compute token vector hash once
// -----------------------------------------------------------------------------
inline size_t hashTokens(const std::vector<int>& tokens) {
  size_t h = tokens.size();
  for(int v : tokens) {
    h ^= std::hash<int>()(v) + 0x9e3779b9 + (h << 6) + (h >> 2);
  }
  return h;
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
  const int nRows;
  const int matrixSize;
  
  // Computed Node Sizes
  std::vector<int> nodeSizes;
  
  // Precomputed child lookup
  std::vector<std::pair<int, int>> children;
  
  // Memoization Cache - now keyed by hash for faster lookup
  // Store actual token vectors separately to handle hash collisions
  std::unordered_map<CacheKey, double*, CacheKeyHash> cache;
  
  // Store token vectors to keep them alive for cache keys
  std::vector<std::vector<int>> storedTokens;
  
  // Memory pool for result matrices
  std::vector<double*> memoryPool;
  
  TreeCounter(IntegerMatrix edge_, int nTip_, int maxSteps_, int nStates_) 
    : edge(edge_), nTip(nTip_), maxSteps(maxSteps_), 
      nRows((1 << nStates_)),
                  matrixSize(nRows * (maxSteps_ + 1)) {
    
    int maxNode = nTip + edge.nrow();
    nodeSizes.resize(maxNode + 1, 0);
    children.resize(maxNode + 1, {-1, -1});
    
    for(int i = 1; i <= nTip; ++i) nodeSizes[i] = 1;
    
    // Precompute children
    for(int i = 0; i < edge.nrow(); ++i) {
      int parent = edge(i, 0);
      int child = edge(i, 1);
      if (children[parent].first == -1) {
        children[parent].first = child;
      } else {
        children[parent].second = child;
      }
    }
    
    computeNodeSizes(nTip + 1);
    
    cache.reserve(2000);
    storedTokens.reserve(2000);
  }
  
  ~TreeCounter() {
    for(auto ptr : memoryPool) {
      delete[] ptr;
    }
  }
  
  int computeNodeSizes(int node) {
    if (node <= nTip) return 1;
    if (nodeSizes[node] > 0) return nodeSizes[node];
    
    int size = 0;
    if (children[node].first != -1) {
      size += computeNodeSizes(children[node].first);
      if (children[node].second != -1) {
        size += computeNodeSizes(children[node].second);
      }
    }
    nodeSizes[node] = size;
    return size;
  }
  
  double* allocateMatrix() {
    double* mat = new double[matrixSize];
    memoryPool.push_back(mat);
    return mat;
  }
  
  inline void initMatrix(double* mat) {
    std::fill(mat, mat + matrixSize, LOG_ZERO);
  }
  
  // ---------------------------------------------------------------------------
  // Generate splits - optimized to reduce allocations
  // ---------------------------------------------------------------------------
  void generateSplits(size_t tokenTypeIdx, int currentLeftSize, int targetLeftSize,
                      const std::vector<int>& totalTokens,
                      std::vector<int>& currentSplit,
                      std::vector<std::vector<int>>& validSplits) {
    
    if (tokenTypeIdx == totalTokens.size()) {
      if (currentLeftSize == targetLeftSize) {
        validSplits.push_back(currentSplit);
      }
      return;
    }
    
    int remainingCapacity = targetLeftSize - currentLeftSize;
    if (remainingCapacity < 0) return;
    
    int available = totalTokens[tokenTypeIdx];
    
    int futureTokens = 0;
    for(size_t i = tokenTypeIdx + 1; i < totalTokens.size(); ++i) 
      futureTokens += totalTokens[i];
    
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
  double* recurse(int node, const std::vector<int>& currentTokens, size_t tokenHash) {
    
    // 1. Check Cache
    CacheKey key(node, currentTokens, tokenHash);
    auto it = cache.find(key);
    if (it != cache.end()) {
      return it->second;
    }
    
    // Allocate and initialize result matrix
    double* resMat = allocateMatrix();
    initMatrix(resMat);
    
    // 2. Base Case: Leaf Node
    if (node <= nTip) {
      int tokenIndex = -1;
      for(size_t i = 0; i < currentTokens.size(); ++i) {
        if (currentTokens[i] == 1) {
          tokenIndex = i;
          break;
        }
      }
      
      if (tokenIndex != -1) {
        int stateMask = 1 << tokenIndex;
        resMat[stateMask] = 0.0;
      }
      
      // Store token vector for cache key
      storedTokens.push_back(currentTokens);
      CacheKey storedKey(node, storedTokens.back(), tokenHash);
      cache.emplace(storedKey, resMat);
      return resMat;
    }
    
    // 3. Recursive Step: Internal Node
    int leftNode = children[node].first;
    int rightNode = children[node].second;
    int leftSize = nodeSizes[leftNode];
    
    // Generate Splits
    std::vector<std::vector<int>> splits;
    splits.reserve(100);
    std::vector<int> buffer(currentTokens.size());
    generateSplits(0, 0, leftSize, currentTokens, buffer, splits);
    
    // Preallocate vectors to avoid repeated allocation
    std::vector<int> tokensR(currentTokens.size());
    
    // Process Splits
    for (const auto& tokensL : splits) {
      // Calculate tokensR inline
      for(size_t i = 0; i < tokensR.size(); ++i) 
        tokensR[i] = currentTokens[i] - tokensL[i];
      
      // Precompute hashes to avoid recomputing in recursive calls
      size_t hashL = hashTokens(tokensL);
      size_t hashR = hashTokens(tokensR);
      
      double* matL = recurse(leftNode, tokensL, hashL);
      double* matR = recurse(rightNode, tokensR, hashR);
      
      // Convolve - unroll step loop for better performance
      for (int maskL = 1; maskL < nRows; ++maskL) {
        // Check if this mask has any non-zero values
        bool hasNonZero = false;
        for (int stepL = 0; stepL <= maxSteps; ++stepL) {
          if (matL[maskL + nRows * stepL] != LOG_ZERO) {
            hasNonZero = true;
            break;
          }
        }
        if (!hasNonZero) continue;
        
        for (int stepL = 0; stepL <= maxSteps; ++stepL) {
          double valL = matL[maskL + nRows * stepL];
          if (valL == LOG_ZERO) continue;
          
          int maxStepR = maxSteps - stepL;
          
          for (int maskR = 1; maskR < nRows; ++maskR) {
            // Precompute Fitch operation outside step loop
            FitchResult fitch = fitchOp(maskL, maskR);
            int baseCost = stepL + fitch.cost;
            
            if (baseCost > maxSteps) continue;
            
            for (int stepR = 0; stepR <= maxStepR; ++stepR) {
              double valR = matR[maskR + nRows * stepR];
              if (valR == LOG_ZERO) continue;
              
              int totalSteps = baseCost + stepR;
              
              if (totalSteps <= maxSteps) {
                int targetIdx = fitch.state + nRows * totalSteps;
                resMat[targetIdx] = logAdd(resMat[targetIdx], valL + valR);
              }
            }
          }
        }
      }
    }
    
    // Store token vector for cache key
    storedTokens.push_back(currentTokens);
    CacheKey storedKey(node, storedTokens.back(), tokenHash);
    cache.emplace(storedKey, resMat);
    return resMat;
  }
};

// -----------------------------------------------------------------------------
// Main Exported Function
// -----------------------------------------------------------------------------
//' Number of leaf labellings producing given tree lengths
//' 
//' @param tree A binary tree of class phylo
//' @param tokens Integer vector: Occurrences of each token
//' @param steps Integer vector: maximum number of steps to compute.
//' A negative value computes all possible step counts.
//' 
//' Given a tree, how many distinct leaf labellings produce a tree length of
//' $k$ under Fitch parsimony?  
//' The number of leaves exhibiting each character state is given by `tokens`.
//' 
//' 
//' @examples
//' tree <- TreeTools::BalancedTree(7)
//' tokens <- c(2, 3, 2) # e.g. 0 0 1 1 1 2 2
//' FixedTreeCount(tree, 2:4, tokens)
//' 
//' Note: setting `Inf` for steps will give all possible outcomes.
//' Setting a lower value will allow some recursions to terminate early,
//' potentially improving runtime - but probably not by much.
//' 
//' @returns `FixedTreeCount()` returns a vector of names `0`...`maxSteps`,
//' where each entry lists the natural logarithm of the number of distinct
//' labellings that produce that Fitch length.
//' 
//' @importFrom TreeTools CladeSizes
//' @export
// [[Rcpp::export]]
NumericVector FixedTreeCount(List tree, std::vector<int> tokens,
                             double steps = -1.0) {
  
  IntegerMatrix edge = tree["edge"];
  bool in_preorder = false;
  if (tree.containsElementNamed("order")) {
    const CharacterVector order = tree["order"];
    in_preorder = (order[0] == "Preorder");
  }
  
  if (!in_preorder) {
    edge = TreeTools::preorder_edges_and_nodes(edge(_, 0), edge(_, 1));
  }
  
  int nTip = 0;
  
  if (tree.containsElementNamed("tip.label")) {
    CharacterVector tips = tree["tip.label"];
    nTip = tips.size();
  } else {
    stop("Tree lacks tip.label");
  }
  
  std::vector<int> activeTokens;
  int sumTokens = 0;
  for(int t : tokens) {
    if (t > 0) {
      activeTokens.push_back(t);
      sumTokens += t;
    }
    if (t < 0) {
      stop("`tokens` may not be negative");
    }
  }
  std::sort(activeTokens.begin(), activeTokens.end());
  
  if (sumTokens != nTip) {
    stop("Number of leaves does not match total tokens.");
  }
  
  int calcMaxSteps = nTip - activeTokens.size() + 1;
  int limitSteps = (steps < 0) ? calcMaxSteps : (int)steps;
  int actualMaxSteps = std::min(limitSteps, calcMaxSteps);
  
  TreeCounter solver(edge, nTip, actualMaxSteps, activeTokens.size());
  
  int rootNode = nTip + 1;
  size_t rootHash = hashTokens(activeTokens);
  double* finalMat = solver.recurse(rootNode, activeTokens, rootHash);
  
  NumericVector results(actualMaxSteps + 1);
  
  for (int s = 0; s <= actualMaxSteps; ++s) {
    double stepTotal = LOG_ZERO;
    int offset = solver.nRows * s;
    for (int m = 1; m < solver.nRows; ++m) {
      stepTotal = logAdd(stepTotal, finalMat[m + offset]);
    }
    results[s] = stepTotal;
  }
  
  CharacterVector names(actualMaxSteps + 1);
  for(int i=0; i<=actualMaxSteps; ++i) 
    names[i] = std::to_string(i);
  results.attr("names") = names;
  
  return results;
}

//' @describeIn FixedTreeCount
//' @param tokenMatrix Integer matrix in which each column corresponds to a
//' character.
// [[Rcpp::export]]
NumericMatrix FixedTreeCountBatch(List tree, IntegerMatrix tokenMatrix,
                                  double steps = -1.0) {
  // Extract tree edge matrix
  IntegerMatrix edge = tree["edge"];
  bool in_preorder = false;
  if (tree.containsElementNamed("order")) {
    const CharacterVector order = tree["order"];
    in_preorder = (order[0] == "Preorder");
  }
  if (!in_preorder) {
    edge = TreeTools::preorder_edges_and_nodes(edge(_, 0), edge(_, 1));
  }
  
  // Get number of tips
  int nTip = 0;
  if (tree.containsElementNamed("tip.label")) {
    CharacterVector tips = tree["tip.label"];
    nTip = tips.size();
  } else {
    stop("Tree lacks tip.label");
  }
  
  // Dimensions
  int nStates = tokenMatrix.nrow();
  int nChars = tokenMatrix.ncol();
  
  // Compute max steps
  int calcMaxSteps = nTip - nStates + 1;
  int limitSteps = (steps < 0) ? calcMaxSteps : (int)steps;
  int actualMaxSteps = std::min(limitSteps, calcMaxSteps);
  
  // Initialize solver once
  TreeCounter solver(edge, nTip, actualMaxSteps, nStates);
  int rootNode = nTip + 1;
  
  // Prepare output matrix: rows = steps, cols = characters
  NumericMatrix results(actualMaxSteps + 1, nChars);
  
  // Loop over characters
  for (int c = 0; c < nChars; ++c) {
    std::vector<int> tokens(nStates);
    int sumTokens = 0;
    for (int r = 0; r < nStates; ++r) {
      int val = tokenMatrix(r, c);
      if (val < 0) stop("`tokens` may not be negative");
      tokens[r] = val;
      sumTokens += val;
    }
    if (sumTokens != nTip) {
      stop("Number of leaves does not match total tokens for character " + std::to_string(c + 1));
    }
    
    // Sort tokens for consistency
    std::sort(tokens.begin(), tokens.end());
    
    // Compute hash and recurse
    size_t rootHash = hashTokens(tokens);
    double* finalMat = solver.recurse(rootNode, tokens, rootHash);
    
    // Aggregate results for this character
    for (int s = 0; s <= actualMaxSteps; ++s) {
      double stepTotal = LOG_ZERO;
      int offset = solver.nRows * s;
      for (int m = 1; m < solver.nRows; ++m) {
        stepTotal = logAdd(stepTotal, finalMat[m + offset]);
      }
      results(s, c) = stepTotal; // log counts
    }
  }
  
  // Assign row names (steps)
  CharacterVector rowNames(actualMaxSteps + 1);
  for (int i = 0; i <= actualMaxSteps; ++i) {
    rowNames[i] = std::to_string(i);
  }
  results.attr("dimnames") = List::create(rowNames, R_NilValue);
  
  return results;
}
