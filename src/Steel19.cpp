#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <map>

using namespace Rcpp;

// Constants
const int MAX_R = 8; // Hard limit for safety (2^8 = 256 states)
const double NEG_INF = -std::numeric_limits<double>::infinity();

// --- User-Provided LSE Accumulator ---
struct LSEAccumulator {
  long double maxv;
  long double acc; // accumulates sum of exp(v - maxv)
  bool empty;
  
  LSEAccumulator() : maxv(-std::numeric_limits<long double>::infinity()), acc(0.0L), empty(true) {}
  
  inline void add(long double v) {
    // ignore non-finite terms
    if (!std::isfinite((double)v)) return;
    
    if (empty) {
      maxv = v;
      acc = 1.0L;        // exp(v - maxv) = exp(0) = 1
      empty = false;
    } else {
      if (v > maxv) {
        // bring old accumulator to the scale of the new max
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
  
  inline bool is_empty() const { return empty; }
};

// --- Internal Data Structures ---

// Represents the distribution at a node:
// dp[mask][cost] = log_probability
// We use a flattened vector for cache locality + a list of active masks to skip zeros.
struct NodeProfile {
  // Table dimensions: (1 << r) rows, (k + 1) cols
  // Indexing: flat_table[mask * (max_k + 1) + cost]
  std::vector<double> flat_table; 
  std::vector<int> active_masks; // Only indices with finite probs
  int max_k;
  int num_states_r;
  int num_masks;
  
  NodeProfile(int r, int k) : max_k(k), num_states_r(r) {
    num_masks = 1 << r;
    flat_table.assign(num_masks * (max_k + 1), NEG_INF);
  }
  
  // Helper to access the table
  inline double get(int mask, int cost) const {
    return flat_table[mask * (max_k + 1) + cost];
  }
  
  // Helper to set (only used during initialization really)
  inline void set(int mask, int cost, double val) {
    if (flat_table[mask * (max_k + 1) + cost] == NEG_INF) {
      // If this is the first time we touch this mask, add to active list
      bool mask_exists = false;
      // Check if mask is already active (simple linear scan is fine for small r)
      // Note: For r=6, max 64 masks. 
      for(int m : active_masks) if(m == mask) { mask_exists = true; break; }
      if(!mask_exists) active_masks.push_back(mask);
    }
    flat_table[mask * (max_k + 1) + cost] = val;
  }
};

// --- Core Recursive Logic ---

// Computes the profile for a parent node based on two children
// This implements the Fitch operation logic in Log-Probability space
NodeProfile combine_nodes(const NodeProfile& left, const NodeProfile& right, int max_k, int r) {
  NodeProfile parent(r, max_k);
  
  // We need to accumulate probabilities into the parent table.
  // Since we are adding probabilities (log-sum-exp), we need LSE accumulators.
  // However, creating a full grid of LSEAccumulators is expensive.
  // We will iterate children, compute target indices, and push to a temp map or direct update?
  // Better: Use a flat vector of LSEAccumulators for the parent.
  std::vector<LSEAccumulator> accumulators(parent.num_masks * (max_k + 1));
  std::vector<bool> accum_touched(parent.num_masks, false); // Track which masks got hit
  
  // Iterate over ACTIVE masks of children only (Sparse Optimization)
  for (int mask_L : left.active_masks) {
    for (int mask_R : right.active_masks) {
      
      // Fitch Operation Logic
      int intersection = mask_L & mask_R;
      int mask_new;
      int cost_inc;
      
      if (intersection != 0) {
        mask_new = intersection;
        cost_inc = 0;
      } else {
        mask_new = mask_L | mask_R;
        cost_inc = 1;
      }
      
      // Convolution of Costs
      // We need: cost_L + cost_R + cost_inc <= max_k
      // Iterate valid costs for Left
      for (int cL = 0; cL <= max_k; ++cL) {
        double pL = left.get(mask_L, cL);
        if (!std::isfinite(pL)) continue;
        
        // Iterate valid costs for Right
        // cR <= max_k - cL - cost_inc
        int max_cR = max_k - cL - cost_inc;
        
        for (int cR = 0; cR <= max_cR; ++cR) {
          double pR = right.get(mask_R, cR);
          if (!std::isfinite(pR)) continue;
          
          int total_cost = cL + cR + cost_inc;
          
          // Add log_prob to accumulator
          // P(new) += P(L) * P(R)  =>  logP(new) = LSE(logP(new), logP(L) + logP(R))
          accumulators[mask_new * (max_k + 1) + total_cost].add(pL + pR);
          accum_touched[mask_new] = true;
        }
      }
    }
  }
  
  // Finalize parent profile
  for (int m = 0; m < parent.num_masks; ++m) {
    if (accum_touched[m]) {
      parent.active_masks.push_back(m);
      for (int c = 0; c <= max_k; ++c) {
        parent.flat_table[m * (max_k + 1) + c] = accumulators[m * (max_k + 1) + c].result();
      }
    }
  }
  
  return parent;
}


// --- Main Exported Function ---

//' Exact Distribution of Parsimony Score on a Tree (I.I.D. Model)
//' 
//' Computes the log-probability that a character generated by independent 
//' resampling of leaf states has a parsimony length L(T) = k.
//' 
//' @param tree A phylo object (list with edge matrix).
//' @param state_probs A numeric vector of probabilities for each state (0..r-1).
//'        Must sum to 1. Length determines r.
//' @param steps The maximum cost k to compute (returns 0..k).
//' 
//' @return A named numeric vector of log-probabilities for steps 0 to k.
//' @export
// [[Rcpp::export]]
NumericVector active_parsimony_dist(List tree, NumericVector state_probs, int steps) {
  
  // 1. Validate Inputs
  int r = state_probs.size();
  if (r > MAX_R) {
    stop("Number of states (r) exceeds maximum supported limit of 8.");
  }
  if (std::abs(sum(state_probs) - 1.0) > 1e-6) {
    warning("State probabilities do not sum to 1.0; results may be unnormalized.");
  }
  
  IntegerMatrix edge = tree["edge"];
  int n_edge = edge.nrow();
  
  // Determine N_tips and N_nodes
  // In ape, tips are 1..Ntip. Internal are Ntip+1..
  // We scan edge matrix to be safe.
  int max_id = 0;
  for(int i=0; i < n_edge * 2; ++i) {
    if(edge[i] > max_id) max_id = edge[i];
  }
  
  // Build Adjacency (Children map)
  // children[u] = {child1, child2}
  // Note: R indices are 1-based. We will shift to 0-based for C++ logic if needed, 
  // but keeping 1-based for node IDs is less confusing relative to R.
  std::vector<std::vector<int>> children(max_id + 1);
  std::vector<bool> is_child(max_id + 1, false);
  
  for (int i = 0; i < n_edge; ++i) {
    int u = edge(i, 0); // Parent
    int v = edge(i, 1); // Child
    children[u].push_back(v);
    is_child[v] = true;
  }
  
  // Identify Root and Tips
  int root = -1;
  for (int i = 1; i <= max_id; ++i) {
    if (!is_child[i] && !children[i].empty()) {
      root = i; 
      break;
    }
    if (children[i].empty()) {
      // It's a tip (assuming connected tree)
    }
  }
  
  if (root == -1) stop("Could not identify root node.");
  
  // 2. Post-Order Traversal Strategy
  // We need to process children before parents.
  // A simple DFS post-order.
  std::vector<int> post_order;
  std::vector<int> stack;
  stack.push_back(root);
  
  // Iterative Post-Order Generation
  // We'll use a 2-stack approach or simple recursive helper. Recursive is fine for tree depth < 10k.
  // Let's use an explicit stack for safety.
  // Actually, ape "edge" matrix is usually already sorted or easily sortable.
  // But let's do a topological sort to be robust.
  std::vector<int> visit_stack;
  visit_stack.push_back(root);
  while(!visit_stack.empty()) {
    int u = visit_stack.back();
    visit_stack.pop_back();
    post_order.push_back(u); // Pre-order push
    
    // Push children
    for(int child : children[u]) {
      visit_stack.push_back(child);
    }
  }
  // Reverse to get post-order (children first)
  std::reverse(post_order.begin(), post_order.end());
  
  // 3. DP Execution
  // We store calculated profiles in a map or vector. 
  // Since node IDs are dense 1..max_id, vector is fine.
  // We use pointers to NodeProfile to avoid copying huge vectors, 
  // but since we process bottom-up, we can construct/move.
  
  // Store *active* profiles. Once a parent is computed, children can be discarded.
  // However, `max_id` isn't huge, keeping them is easier for implementation unless RAM is tight.
  // Let's use a vector of unique_ptrs or just objects.
  std::vector<NodeProfile> profiles; 
  profiles.reserve(max_id + 1);
  
  // Initialize with dummy
  for(int i=0; i<=max_id; ++i) profiles.emplace_back(r, steps);
  
  for (int u : post_order) {
    if (children[u].empty()) {
      // --- LEAF NODE ---
      // Initialize based on I.I.D. probabilities
      // For a random character, the leaf has probability state_probs[s] of being in state s.
      // Fitch state is exactly {s}, i.e., mask (1 << s). Cost is 0.
      
      for (int s = 0; s < r; ++s) {
        double log_prob = std::log(state_probs[s]);
        int mask = (1 << s);
        profiles[u].set(mask, 0, log_prob);
      }
      
    } else {
      // --- INTERNAL NODE ---
      // Assume binary tree for Fitch (standard).
      if (children[u].size() != 2) {
        stop("Tree must be binary (internal nodes must have 2 children). Node %d has %d children.", u, children[u].size());
      }
      
      int child1 = children[u][0];
      int child2 = children[u][1];
      
      profiles[u] = combine_nodes(profiles[child1], profiles[child2], steps, r);
      
      // Optional: Free memory of children to save RAM? 
      // profiles[child1] = NodeProfile(0,0); 
    }
  }
  
  // 4. Aggregate Results at Root
  // Result is sum over all Fitch States for each length
  NodeProfile& root_prof = profiles[root];
  
  NumericVector results(steps + 1);
  
  for (int k = 0; k <= steps; ++k) {
    LSEAccumulator final_acc;
    for (int mask : root_prof.active_masks) {
      final_acc.add(root_prof.get(mask, k));
    }
    results[k] = final_acc.result();
  }
  
  // Assign names 0..k
  CharacterVector names(steps + 1);
  for(int i=0; i<=steps; ++i) names[i] = std::to_string(i);
  results.names() = names;
  
  return results;
 }
