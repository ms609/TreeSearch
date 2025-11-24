#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <string>

using namespace Rcpp;

// Constants
const int MAX_R = 8; // Hard limit for safety (2^8 = 256 states)
const double NEG_INF = -std::numeric_limits<double>::infinity();

// ============================================================================
// PART 1: Helpers & Structures
// ============================================================================

// --- LSE Accumulator for Log-Space (Used in Distribution) ---
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
    return (double)(maxv + logl(acc));
  }
};

// --- Moments Accumulator (Linear Space) ---
// Used to aggregate conditional expectations for a specific Fitch mask
struct MomentAccumulator {
  double prob_sum;
  double w_len_sum;    // Weighted sum of lengths: Sum(P * L)
  double w_len_sq_sum; // Weighted sum of squared lengths: Sum(P * L^2)
  
  MomentAccumulator() : prob_sum(0.0), w_len_sum(0.0), w_len_sq_sum(0.0) {}
  
  // Add a contribution from a child-pair combination
  // p_pair: Probability of this specific Left/Right state combo
  // e1_pair: E[L] for this combo
  // e2_pair: E[L^2] for this combo
  void add(double p_pair, double e1_pair, double e2_pair) {
    prob_sum += p_pair;
    w_len_sum += p_pair * e1_pair;
    w_len_sq_sum += p_pair * e2_pair;
  }
};

// --- Tree Traversal Helper ---
// Returns post-order node indices and children map
struct TreeTopo {
  int root;
  int max_node;
  std::vector<int> post_order;
  std::vector<std::vector<int>> children;
};

TreeTopo process_tree(const IntegerMatrix& edge) {
  TreeTopo topo;
  int n_edge = edge.nrow();
  topo.max_node = 0;
  
  for(int i=0; i < n_edge * 2; ++i) {
    if(edge[i] > topo.max_node) topo.max_node = edge[i];
  }
  
  topo.children.resize(topo.max_node + 1);
  std::vector<bool> is_child(topo.max_node + 1, false);
  
  for (int i = 0; i < n_edge; ++i) {
    int u = edge(i, 0); 
    int v = edge(i, 1); 
    topo.children[u].push_back(v);
    is_child[v] = true;
  }
  
  topo.root = -1;
  for (int i = 1; i <= topo.max_node; ++i) {
    if (!is_child[i] && !topo.children[i].empty()) {
      topo.root = i; break;
    }
  }
  if (topo.root == -1) stop("Invalid tree: could not find root.");
  
  std::vector<int> stack;
  stack.push_back(topo.root);
  while(!stack.empty()) {
    int u = stack.back();
    stack.pop_back();
    topo.post_order.push_back(u); 
    for(int child : topo.children[u]) stack.push_back(child);
  }
  std::reverse(topo.post_order.begin(), topo.post_order.end());
  
  return topo;
}

// ============================================================================
// PART 2: Full Distribution (Log-Space DP)
// ============================================================================

struct DistNodeProfile {
  std::vector<double> flat_table; 
  std::vector<int> active_masks;
  int max_k, num_masks;
  
  DistNodeProfile(int r, int k) : max_k(k) {
    num_masks = 1 << r;
    flat_table.assign(num_masks * (max_k + 1), NEG_INF);
  }
  
  inline double get(int mask, int cost) const { return flat_table[mask * (max_k + 1) + cost]; }
  
  inline void set(int mask, int cost, double val) {
    if (flat_table[mask * (max_k + 1) + cost] == NEG_INF) {
      // Check if mask exists in active list (linear scan ok for small r)
      bool exists = false;
      for(int m : active_masks) if(m == mask) { exists = true; break; }
      if(!exists) active_masks.push_back(mask);
    }
    flat_table[mask * (max_k + 1) + cost] = val;
  }
};

DistNodeProfile combine_dist(const DistNodeProfile& L, const DistNodeProfile& R, int max_k, int r) {
  DistNodeProfile P(r, max_k);
  std::vector<LSEAccumulator> accs(P.num_masks * (max_k + 1));
  std::vector<bool> touched(P.num_masks, false);
  
  for (int mL : L.active_masks) {
    for (int mR : R.active_masks) {
      int intersection = mL & mR;
      int mNew = (intersection != 0) ? intersection : (mL | mR);
      int cost_inc = (intersection == 0) ? 1 : 0;
      
      for (int cL = 0; cL <= max_k; ++cL) {
        double pL = L.get(mL, cL);
        if (!std::isfinite(pL)) continue;
        
        int max_cR = max_k - cL - cost_inc;
        for (int cR = 0; cR <= max_cR; ++cR) {
          double pR = R.get(mR, cR);
          if (!std::isfinite(pR)) continue;
          
          accs[mNew * (max_k + 1) + (cL + cR + cost_inc)].add(pL + pR);
          touched[mNew] = true;
        }
      }
    }
  }
  
  for (int m = 0; m < P.num_masks; ++m) {
    if (touched[m]) {
      P.active_masks.push_back(m);
      for (int c = 0; c <= max_k; ++c) P.flat_table[m * (max_k + 1) + c] = accs[m * (max_k + 1) + c].result();
    }
  }
  return P;
}

//' Exact Distribution of Parsimony Score on a Tree (I.I.D. Model)
//'
//' `active_parsimony_dist()` computes the log-probability that a character generated by independent 
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
  int r = state_probs.size();
  if (r > MAX_R) stop("r exceeds max limit");
  TreeTopo topo = process_tree(tree["edge"]);
  
  std::vector<DistNodeProfile> profiles;
  profiles.reserve(topo.max_node + 1);
  for(int i=0; i<=topo.max_node; ++i) profiles.emplace_back(r, steps);
  
  for (int u : topo.post_order) {
    if (topo.children[u].empty()) { // Leaf
      for (int s = 0; s < r; ++s) {
        profiles[u].set(1 << s, 0, std::log(state_probs[s]));
      }
    } else { // Internal
      profiles[u] = combine_dist(profiles[topo.children[u][0]], profiles[topo.children[u][1]], steps, r);
    }
  }
  
  NumericVector res(steps + 1);
  for (int k = 0; k <= steps; ++k) {
    LSEAccumulator final;
    for (int m : profiles[topo.root].active_masks) final.add(profiles[topo.root].get(m, k));
    res[k] = final.result();
  }
  CharacterVector names(steps + 1);
  for(int i=0; i<=steps; ++i) names[i] = std::to_string(i);
  res.names() = names;
  return res;
}

// ============================================================================
// PART 3: Moments (Expectation & Variance) (Linear Space)
// ============================================================================

struct MomentNodeProfile {
  std::vector<double> prob; // P(S=mask)
  std::vector<double> e1;   // E[L | S=mask]
  std::vector<double> e2;   // E[L^2 | S=mask]
  std::vector<int> active_masks;
  int num_masks;
  
  MomentNodeProfile(int r) {
    num_masks = 1 << r;
    prob.assign(num_masks, 0.0);
    e1.assign(num_masks, 0.0);
    e2.assign(num_masks, 0.0);
  }
  
  void set_leaf(int mask, double p) {
    active_masks.push_back(mask);
    prob[mask] = p;
    e1[mask] = 0.0;
    e2[mask] = 0.0;
  }
};

MomentNodeProfile combine_moments(const MomentNodeProfile& L, const MomentNodeProfile& R, int r) {
  MomentNodeProfile P(r);
  std::vector<MomentAccumulator> accs(P.num_masks);
  std::vector<bool> touched(P.num_masks, false);
  
  for (int mL : L.active_masks) {
    for (int mR : R.active_masks) {
      double p_joint = L.prob[mL] * R.prob[mR];
      if (p_joint <= 0.0) continue; // Skip unlikely branches
      
      int intersection = mL & mR;
      int mNew = (intersection != 0) ? intersection : (mL | mR);
      double cost_inc = (intersection == 0) ? 1.0 : 0.0;
      
      // Values from children
      double E_L = L.e1[mL];
      double E_R = R.e1[mR];
      double E2_L = L.e2[mL];
      double E2_R = R.e2[mR];
      
      // Formulae for sums of independent random variables (+ constant cost_inc)
      // New Length Y = L + R + c
      // E[Y] = E[L] + E[R] + c
      double E_new = E_L + E_R + cost_inc;
      
      // E[Y^2] = E[(L + R + c)^2]
      //        = E[L^2 + R^2 + c^2 + 2LR + 2Lc + 2Rc]
      //        = E[L^2] + E[R^2] + c^2 + 2E[L]E[R] + 2c(E[L] + E[R])
      double E2_new = E2_L + E2_R + (cost_inc * cost_inc) 
        + 2.0 * E_L * E_R 
      + 2.0 * cost_inc * (E_L + E_R);
      
      accs[mNew].add(p_joint, E_new, E2_new);
      touched[mNew] = true;
    }
  }
  
  // Normalize Accumulators to get Conditional Expectations
  for(int m=0; m < P.num_masks; ++m) {
    if(touched[m]) {
      P.active_masks.push_back(m);
      P.prob[m] = accs[m].prob_sum;
      if (P.prob[m] > 0) {
        P.e1[m] = accs[m].w_len_sum / P.prob[m];
        P.e2[m] = accs[m].w_len_sq_sum / P.prob[m];
      }
    }
  }
  return P;
}

//' Expectation and Variance of Parsimony Score
//'
//' `parsimony_moments()` efficiently computes the Expected value (E) and Variance (V) of the parsimony 
//' length L(T) under the I.I.D. model.
//'
//' @rdname active_parsimony_dist
//' @return `parsimony_moments` returns a list with components `expectation` and `variance`.
//' @export
// [[Rcpp::export]]
List parsimony_moments(List tree, NumericVector state_probs) {
  int r = state_probs.size();
  if (r > MAX_R) stop("r exceeds max limit");
  if (std::abs(sum(state_probs) - 1.0) > 1e-6) warning("Probs do not sum to 1");
  
  TreeTopo topo = process_tree(tree["edge"]);
  
  std::vector<MomentNodeProfile> profiles;
  profiles.reserve(topo.max_node + 1);
  for(int i=0; i<=topo.max_node; ++i) profiles.emplace_back(r);
  
  for (int u : topo.post_order) {
    if (topo.children[u].empty()) {
      for (int s = 0; s < r; ++s) {
        profiles[u].set_leaf(1 << s, state_probs[s]);
      }
    } else {
      profiles[u] = combine_moments(profiles[topo.children[u][0]], profiles[topo.children[u][1]], r);
    }
  }
  
  // Aggregate at Root
  double grand_E1 = 0.0;
  double grand_E2 = 0.0;
  double total_prob = 0.0;
  
  MomentNodeProfile& root = profiles[topo.root];
  for (int m : root.active_masks) {
    grand_E1 += root.prob[m] * root.e1[m];
    grand_E2 += root.prob[m] * root.e2[m];
    total_prob += root.prob[m];
  }
  
  // Normalize (in case total_prob drifted slightly from 1.0 due to float math, though unlikely)
  grand_E1 /= total_prob;
  grand_E2 /= total_prob;
  
  double variance = grand_E2 - (grand_E1 * grand_E1);
  
  return List::create(
    Named("expectation") = grand_E1,
    Named("variance") = variance
  );
}
