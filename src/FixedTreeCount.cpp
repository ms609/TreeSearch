#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <limits>
#include <cstring>

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
// Custom hash function for pair<int, vector<int>> keys
// -----------------------------------------------------------------------------
struct PairVectorHash {
  std::size_t operator()(const std::pair<int, std::vector<int>>& p) const {
    std::size_t seed = std::hash<int>()(p.first);
    seed ^= p.second.size() + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    for(auto& i : p.second) {
      seed ^= std::hash<int>()(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

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
  const int matrixSize; // nRows * (maxSteps + 1)
  
  // Computed Node Sizes
  std::vector<int> nodeSizes;
  
  // Memoization Cache - using unordered_map for better performance
  std::unordered_map<std::pair<int, std::vector<int>>, double*, 
                     PairVectorHash> cache;
  
  // Memory pool for result matrices to avoid repeated allocations
  std::vector<double*> memoryPool;
  
  TreeCounter(IntegerMatrix edge_, int nTip_, int maxSteps_, int nStates_) 
    : edge(edge_), nTip(nTip_), maxSteps(maxSteps_), 
      nRows((1 << nStates_)),
                  matrixSize(nRows * (maxSteps_ + 1)) {
    
    int maxNode = nTip + edge.nrow();
    nodeSizes.resize(maxNode + 1, 0);
    
    for(int i = 1; i <= nTip; ++i) nodeSizes[i] = 1;
    
    computeNodeSizes(nTip + 1);
  }
  
  ~TreeCounter() {
    // Clean up memory pool
    for(auto ptr : memoryPool) {
      delete[] ptr;
    }
  }
  
  int computeNodeSizes(int node) {
    if (node <= nTip) return 1;
    if (nodeSizes[node] > 0) return nodeSizes[node];
    
    int size = 0;
    for(int i = 0; i < edge.nrow(); ++i) {
      if (edge(i, 0) == node) {
        size += computeNodeSizes(edge(i, 1));
      }
    }
    nodeSizes[node] = size;
    return size;
  }
  
  // Allocate a new matrix from pool
  double* allocateMatrix() {
    double* mat = new double[matrixSize];
    memoryPool.push_back(mat);
    return mat;
  }
  
  // Initialize matrix to LOG_ZERO
  inline void initMatrix(double* mat) {
    std::fill(mat, mat + matrixSize, LOG_ZERO);
  }
  
  // ---------------------------------------------------------------------------
  // Recursively generate splits of tokens
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
  // Core Recursive Function - now returns pointer instead of copying
  // ---------------------------------------------------------------------------
  double* recurse(int node, const std::vector<int>& currentTokens) {
    
    // 1. Check Cache
    std::pair<int, std::vector<int>> key = {node, currentTokens};
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
        int idx = stateMask + nRows * 0;
        resMat[idx] = 0.0;
      }
      
      cache[key] = resMat;
      return resMat;
    }
    
    // 3. Recursive Step: Internal Node
    int leftNode = -1, rightNode = -1;
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
    
    // Preallocate tokensR to avoid repeated allocations
    std::vector<int> tokensR(currentTokens.size());
    
    // Process Splits
    for (const auto& tokensL : splits) {
      // Calculate tokensR
      for(size_t i = 0; i < tokensR.size(); ++i) 
        tokensR[i] = currentTokens[i] - tokensL[i];
      
      double* matL = recurse(leftNode, tokensL);
      double* matR = recurse(rightNode, tokensR);
      
      // Convolve - optimized inner loops
      for (int stepL = 0; stepL <= maxSteps; ++stepL) {
        int offsetL = nRows * stepL;
        
        for (int maskL = 1; maskL < nRows; ++maskL) {
          double valL = matL[maskL + offsetL];
          if (valL == LOG_ZERO) continue;
          
          int maxStepR = maxSteps - stepL;
          for (int stepR = 0; stepR <= maxStepR; ++stepR) {
            int offsetR = nRows * stepR;
            
            for (int maskR = 1; maskR < nRows; ++maskR) {
              double valR = matR[maskR + offsetR];
              if (valR == LOG_ZERO) continue;
              
              FitchResult fitch = fitchOp(maskL, maskR);
              int totalSteps = stepL + stepR + fitch.cost;
              
              if (totalSteps <= maxSteps) {
                int targetIdx = fitch.state + nRows * totalSteps;
                double combinedProb = valL + valR;
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
  
  // 1. Setup
  IntegerMatrix edge = tree["edge"];
  int nTip = 0;
  
  if (tree.containsElementNamed("tip.label")) {
    CharacterVector tips = tree["tip.label"];
    nTip = tips.size();
  } else {
    for(int i=0; i<edge.length(); ++i) 
      if(edge[i] > nTip) nTip = edge[i];
  }
  
  // Clean tokens
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
  int calcMaxSteps = nTip - activeTokens.size() + 1;
  int limitSteps = (steps < 0) ? calcMaxSteps : (int)steps;
  int actualMaxSteps = std::min(limitSteps, calcMaxSteps);
  
  // 2. Instantiate Logic
  TreeCounter solver(edge, nTip, actualMaxSteps, activeTokens.size());
  
  // 3. Run
  int rootNode = nTip + 1;
  double* finalMat = solver.recurse(rootNode, activeTokens);
  
  // 4. Summarize Results
  NumericVector results(actualMaxSteps + 1);
  
  for (int s = 0; s <= actualMaxSteps; ++s) {
    double stepTotal = LOG_ZERO;
    int offset = solver.nRows * s;
    for (int m = 1; m < solver.nRows; ++m) {
      stepTotal = logAdd(stepTotal, finalMat[m + offset]);
    }
    results[s] = stepTotal;
  }
  
  // Assign names
  CharacterVector names(actualMaxSteps + 1);
  for(int i=0; i<=actualMaxSteps; ++i) 
    names[i] = std::to_string(i);
  results.attr("names") = names;
  
  return results;
}
 