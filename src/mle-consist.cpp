// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <algorithm>
#include <random>
#include <cmath>
#include <cstring> // for memset
using namespace Rcpp;

// ---------- Prepare tree metadata (call once from R) ---------------------
// Input: edge matrix (parent, child) 1-based
// Output: list with nNodes, root, postorder (IntegerVector), children_flat, child_offsets, child_counts

// [[Rcpp::export]]
List mlci_prepare_tree(const IntegerMatrix edge) {
  int nEdge = edge.nrow();
  int nNodes = 0;
  for (int i = 0; i < nEdge; ++i) {
    nNodes = std::max(nNodes, (int)edge(i,0));
    nNodes = std::max(nNodes, (int)edge(i,1));
  }
  // Build children lists (temporary)
  std::vector< std::vector<int> > children(nNodes + 1);
  std::vector<int> isChild(nNodes + 1, 0);
  for (int i = 0; i < nEdge; ++i) {
    int p = edge(i,0), c = edge(i,1);
    children[p].push_back(c);
    isChild[c] = 1;
  }
  // find root
  int root = -1;
  for (int v = 1; v <= nNodes; ++v) if (!isChild[v]) { root = v; break; }
  if (root == -1) stop("could not determine root");
  
  // postorder using iterative stack
  std::vector<int> postorder;
  postorder.reserve(nNodes);
  {
    std::vector<int> stack;
    stack.reserve(nNodes);
    std::vector<int> state(nNodes + 1, 0);
    stack.push_back(root);
    while (!stack.empty()) {
      int v = stack.back();
      if (state[v] == 0) {
        state[v] = 1;
        for (int c : children[v]) if (state[c] == 0) stack.push_back(c);
      } else {
        stack.pop_back();
        if (state[v] == 1) { state[v] = 2; postorder.push_back(v); }
      }
    }
  }
  
  // Flatten children into single array and offsets
  std::vector<int> child_counts(nNodes + 1, 0);
  size_t total_children = 0;
  for (int v = 1; v <= nNodes; ++v) {
    child_counts[v] = (int)children[v].size();
    total_children += children[v].size();
  }
  std::vector<int> child_offsets(nNodes + 1, 0);
  int offset = 0;
  for (int v = 1; v <= nNodes; ++v) {
    child_offsets[v] = offset;
    offset += child_counts[v];
  }
  std::vector<int> children_flat; children_flat.resize(total_children);
  // fill flat array
  for (int v = 1; v <= nNodes; ++v) {
    int pos = child_offsets[v];
    for (size_t i = 0; i < children[v].size(); ++i) children_flat[pos + i] = children[v][i];
  }
  
  // Convert vectors to Rcpp objects for returning
  IntegerVector r_post(postorder.size());
  for (size_t i = 0; i < postorder.size(); ++i) r_post[i] = postorder[i];
  IntegerVector r_children_flat(children_flat.size());
  for (size_t i = 0; i < children_flat.size(); ++i) r_children_flat[i] = children_flat[i];
  IntegerVector r_child_offsets(child_offsets.size());
  for (size_t i = 0; i < child_offsets.size(); ++i) r_child_offsets[i] = child_offsets[i];
  IntegerVector r_child_counts(child_counts.size());
  for (size_t i = 0; i < child_counts.size(); ++i) r_child_counts[i] = child_counts[i];
  
  return List::create(
    Named("nNodes") = nNodes,
    Named("root") = root,
    Named("postorder") = r_post,
    Named("children_flat") = r_children_flat,
    Named("child_offsets") = r_child_offsets,
    Named("child_counts") = r_child_counts
  );
}

// ---------- logLik using prepared tree (flat buffers) --------------------
static double logLik_prepared_internal(
    const IntegerMatrix &tip_states,        // Ntip x nPatterns
    const NumericVector &weights, int k,
    const IntegerVector &postorder, const IntegerVector &children_flat,
    const IntegerVector &child_offsets, const IntegerVector &child_counts,
    int root, int Ntip, int nNodes,
    double t
) {
  int nPatterns = tip_states.ncol();
  const double e = std::exp(- (double)k * t);
  const double invk = 1.0 / (double)k;
  
  // cond buffer: (nNodes+1) * k, index (node * k + s)
  std::vector<double> cond_buf((size_t)(nNodes + 1) * (size_t)k);
  // zero with memset (fast)
  std::memset(cond_buf.data(), 0, cond_buf.size() * sizeof(double));
  
  // temp accum vector reused
  std::vector<double> accum((size_t)k);
  
  double totalLogLik = 0.0;
  
  for (int p = 0; p < nPatterns; ++p) {
    // zero only the buffer (we already zeroed above, but we need to zero per pattern)
    std::memset(cond_buf.data(), 0, cond_buf.size() * sizeof(double));
    // set tips: tips are nodes 1..Ntip
    for (int tip = 0; tip < Ntip; ++tip) {
      int state = tip_states(tip, p);
      if (state < 0 || state >= k) stop("tip state out of range");
      cond_buf[(size_t)(tip + 1) * k + state] = 1.0;
    }
    // postorder combine
    for (int idx = 0; idx < postorder.size(); ++idx) {
      int node = postorder[idx];
      if (node <= Ntip) continue;
      // reset accum to 1
      for (int s = 0; s < k; ++s) accum[s] = 1.0;
      int offs = child_offsets[node];
      int cnt  = child_counts[node];
      for (int ci = 0; ci < cnt; ++ci) {
        int child = children_flat[offs + ci];
        // compute meanChild quickly
        double sumChild = 0.0;
        double *child_ptr = &cond_buf[(size_t)child * k];
        for (int s = 0; s < k; ++s) sumChild += child_ptr[s];
        double meanChild = sumChild * invk;
        for (int s = 0; s < k; ++s) {
          double pts = meanChild + e * (child_ptr[s] - meanChild);
          accum[s] *= pts;
        }
      }
      // write accum into cond_buf for node
      double *node_ptr = &cond_buf[(size_t)node * k];
      for (int s = 0; s < k; ++s) node_ptr[s] = accum[s];
    }
    // root site likelihood
    double sumRoot = 0.0;
    double *root_ptr = &cond_buf[(size_t)root * k];
    for (int s = 0; s < k; ++s) sumRoot += root_ptr[s];
    double siteLik = sumRoot * invk;
    if (siteLik <= 0.0) {
      totalLogLik += weights[p] * log(std::numeric_limits<double>::min());
    } else {
      totalLogLik += weights[p] * log(siteLik);
    }
  }
  return totalLogLik;
}

// R-exported wrapper for convenience
// [[Rcpp::export]]
double logLik_equal_t_prepared(const IntegerMatrix tip_states, const NumericVector weights, int k,
                               const List &tree_prep, int Ntip, double t) {
  IntegerVector postorder = tree_prep["postorder"];
  IntegerVector children_flat = tree_prep["children_flat"];
  IntegerVector child_offsets = tree_prep["child_offsets"];
  IntegerVector child_counts = tree_prep["child_counts"];
  int root = as<int>(tree_prep["root"]);
  int nNodes = as<int>(tree_prep["nNodes"]);
  return logLik_prepared_internal(tip_states, weights, k, postorder, children_flat, child_offsets, child_counts, root, Ntip, nNodes, t);
}

// ---------- MLE (prepared) using golden-section calling logLik_prepared_internal ---
static double mle_t_prepared_internal(const IntegerMatrix &tip_states, const NumericVector &weights, int k,
                                      const List &tree_prep, int Ntip, double lower, double upper, double tol) {
  if (lower <= 0) lower = 1e-12;
  if (upper <= lower) stop("upper must be > lower");
  IntegerVector postorder = tree_prep["postorder"];
  IntegerVector children_flat = tree_prep["children_flat"];
  IntegerVector child_offsets = tree_prep["child_offsets"];
  IntegerVector child_counts = tree_prep["child_counts"];
  int root = as<int>(tree_prep["root"]);
  int nNodes = as<int>(tree_prep["nNodes"]);
  
  const double gr = (sqrt(5.0) + 1.0) / 2.0;
  double a = lower, b = upper;
  double c = b - (b - a) / gr;
  double d = a + (b - a) / gr;
  double fc = logLik_prepared_internal(tip_states, weights, k, postorder, children_flat, child_offsets, child_counts, root, Ntip, nNodes, c);
  double fd = logLik_prepared_internal(tip_states, weights, k, postorder, children_flat, child_offsets, child_counts, root, Ntip, nNodes, d);
  int iter = 0, maxiter = 200;
  while (fabs(b - a) > tol * (fabs(c) + fabs(d)) && iter < maxiter) {
    if (fc > fd) {
      b = d; d = c; fd = fc; c = b - (b - a) / gr; fc = logLik_prepared_internal(tip_states, weights, k, postorder, children_flat, child_offsets, child_counts, root, Ntip, nNodes, c);
    } else {
      a = c; c = d; fc = fd; d = a + (b - a) / gr; fd = logLik_prepared_internal(tip_states, weights, k, postorder, children_flat, child_offsets, child_counts, root, Ntip, nNodes, d);
    }
    iter++;
  }
  double t_hat = (b + a) / 2.0;
  return t_hat;
}

// [[Rcpp::export]]
List mle_t_prepared(const IntegerMatrix tip_states, const NumericVector weights, int k,
                    const List &tree_prep, int Ntip, double lower = 1e-8, double upper = 10.0, double tol = 1e-6) {
  double t_hat = mle_t_prepared_internal(tip_states, weights, k, tree_prep, Ntip, lower, upper, tol);
  // compute logLik at t_hat
  IntegerVector postorder = tree_prep["postorder"];
  IntegerVector children_flat = tree_prep["children_flat"];
  IntegerVector child_offsets = tree_prep["child_offsets"];
  IntegerVector child_counts = tree_prep["child_counts"];
  int root = as<int>(tree_prep["root"]);
  int nNodes = as<int>(tree_prep["nNodes"]);
  double logL = logLik_prepared_internal(tip_states, weights, k, postorder, children_flat, child_offsets, child_counts, root, Ntip, nNodes, t_hat);
  return List::create(Named("t_hat") = t_hat, Named("logLik") = logL);
}

// ---------- Resampler using prepared tree ---------------------------------
// [[Rcpp::export]]
List mlci_resample(
    const List &tree_prep,
    const IntegerMatrix edge, // kept for compatibility when calling mle_t_prepared
    int Ntip,
    const IntegerVector &obs_states, // length Ntip, integers 0..k-1
    const IntegerVector &best_states, // length Ntip, 0/1
    
    double precision = 1e-3,
    int maxResample = 10000,
    double lower = 1e-8,
    double upper = 10.0,
    double tol = 1e-6
) {
  int nEdge = edge.nrow();
  // observed tip matrix
  IntegerMatrix tip_obs(Ntip, 1);
  for (int i = 0; i < Ntip; ++i) tip_obs(i,0) = obs_states[i];
  NumericVector w1(1); w1[0] = 1.0;
  // infer k_obs
  int maxv = -1;
  for (int i=0;i<Ntip;++i) if (obs_states[i] > maxv) maxv = obs_states[i];
  int k_obs = maxv + 1;
  // observed MLE
  double t_hat_obs = mle_t_prepared_internal(tip_obs, w1, k_obs, tree_prep, Ntip, lower, upper, tol);
  double score = t_hat_obs * nEdge;
  // best
  IntegerMatrix tip_best(Ntip, 1);
  for (int i = 0; i < Ntip; ++i) tip_best(i,0) = best_states[i];
  double t_hat_best = mle_t_prepared_internal(tip_best, w1, 2, tree_prep, Ntip, lower, upper, tol);
  double bestScore = t_hat_best * nEdge;
  
  // setup permutation vector (copy obs_states)
  std::vector<int> perm((size_t)Ntip);
  for (int i=0;i<Ntip;++i) perm[i] = obs_states[i];
  
  long long n = 0;
  double mean_rm = 0.0, M2 = 0.0;
  
  // allocate integer matrix for permuted tip states reused
  IntegerMatrix tip_r(Ntip, 1);
  NumericVector wr(1); wr[0] = 1.0;
  
  GetRNGstate();
  // Fisher-Yates in-place shuffle using R's RNG
  auto shuffle_inplace = [&](std::vector<int> &v) {
    for (int i = (int)v.size() - 1; i > 0; --i) {
      double u = unif_rand();
      int j = (int)(u * (i + 1));
      int tmp = v[i];
      v[i] = v[j];
      v[j] = tmp;
    }
  };
  
  double rmScore = 0.0;
  double rmMean = NA_REAL, rmSE = NA_REAL, mlci = NA_REAL, mlciSE = NA_REAL;
  while (true) {
    shuffle_inplace(perm);
    for (int i=0;i<Ntip;++i) tip_r(i,0) = perm[i];
    // infer k_r quickly
    int maxvv = -1;
    for (int i=0;i<Ntip;++i) if (perm[i] > maxvv) maxvv = perm[i];
    int k_r = maxvv + 1;
    if (k_r < 1) k_r = 1;
    
    double t_hat_r = mle_t_prepared_internal(tip_r, wr, k_r, tree_prep, Ntip, lower, upper, tol);
    rmScore = t_hat_r * nEdge;
    
    n++;
    double delta = rmScore - mean_rm;
    mean_rm += delta / (double)n;
    double delta2 = rmScore - mean_rm;
    M2 += delta * delta2;
    
    if (n >= 2) {
      double var = M2 / (double)(n - 1);
      double se_mean = std::sqrt(var) / std::sqrt((double)n);
      double denom = (bestScore - mean_rm);
      if (denom == 0.0) mlciSE = R_PosInf;
      else {
        double gprime = (score - bestScore) / (denom * denom);
        mlciSE = std::fabs(gprime) * se_mean;
      }
      mlci = (bestScore - mean_rm == 0.0) ? NA_REAL : (score - mean_rm) / denom;
      rmMean = mean_rm;
      rmSE = se_mean;
      if (mlciSE <= precision) break;
      if (n >= maxResample) break;
    }
  }
  PutRNGstate();
  
  if (n < 2) {
    rmMean = mean_rm;
    rmSE = NA_REAL;
    mlci = (bestScore - mean_rm == 0.0) ? NA_REAL : (score - mean_rm) / (bestScore - mean_rm);
    mlciSE = NA_REAL;
  }
  
  return List::create(
    Named("t_hat_obs") = t_hat_obs,
    Named("score") = score,
    Named("t_hat_best") = t_hat_best,
    Named("bestScore") = bestScore,
    Named("rmMean") = rmMean,
    Named("rmSE") = rmSE,
    Named("mlci") = mlci,
    Named("mlciSE") = mlciSE,
    Named("nResample") = (int)n
  );
}


// [[Rcpp::export]]
DataFrame MLCI_rcpp(IntegerMatrix edge,
                    IntegerMatrix tipStates,
                    IntegerVector bestSplitInt,
                    double precision = 1e-2,
                    int maxResample = 10000) {
  
  int nEdge = edge.nrow();
  int nTip  = tipStates.nrow();
  int nChar = tipStates.ncol();
  
  // Prepare tree once
  List treePrep = mlci_prepare_tree(edge);
  
  // Result vectors
  NumericVector t_hat(nChar), score(nChar), rmMean(nChar),
  rmSE(nChar), bestScore(nChar), mlci(nChar);
  IntegerVector nResampleVec(nChar);
  
  for(int j=0; j<nChar; j++) {
    IntegerVector col = tipStates(_, j);
    int k = *std::max_element(col.begin(), col.end()) + 1;
    
    // Observed character MLE
    IntegerMatrix tip_col(nTip, 1);
    for (int i = 0; i < nTip; ++i) tip_col(i, 0) = col[i];
    double tObs = mle_t_prepared_internal(tip_col, NumericVector::create(1.0), k, treePrep, nTip, 1e-8, 10.0, 1e-6);
    double scoreObs = tObs * nEdge;
    
    // Best split MLE
    int kBest = std::max(2, *std::max_element(bestSplitInt.begin(), bestSplitInt.end()) + 1);
    IntegerMatrix tip_best(nTip, 1);
    for (int i = 0; i < nTip; ++i) tip_best(i, 0) = bestSplitInt[i];
    double tBest = mle_t_prepared_internal(tip_best, NumericVector::create(1.0), kBest, treePrep, nTip, 1e-8, 10.0, 1e-6);
    double scoreBest = tBest * nEdge;
    
    // Resampling
    std::vector<int> tokenVec(col.begin(), col.end());
    std::vector<double> rmScores;
    int nResample = 0;
    std::random_device rd;
    std::mt19937 g(rd());
    IntegerMatrix tip_r(nTip, 1);
    
    while(true) {
      std::shuffle(tokenVec.begin(), tokenVec.end(), g);
      for(int i=0;i<nTip;++i) tip_r(i,0) = tokenVec[i];
      
      int kPerm = *std::max_element(tokenVec.begin(), tokenVec.end()) + 1;
      double tPerm = mle_t_prepared_internal(tip_r, NumericVector::create(1.0), kPerm, treePrep, nTip, 1e-8, 10.0, 1e-6);
      double rmScore = tPerm * nEdge;
      
      rmScores.push_back(rmScore);
      nResample++;
      
      if(nResample >= 2) {
        double mean_rm = std::accumulate(rmScores.begin(), rmScores.end(), 0.0) / nResample;
        double var = 0.0;
        for(double x : rmScores) var += (x - mean_rm) * (x - mean_rm);
        var /= (nResample - 1);
        double se_rm = std::sqrt(var) / std::sqrt(nResample);
        if(se_rm <= precision || nResample >= maxResample) break;
      }
    }
    
    // Compute RM mean/SE
    double rm_mean = std::accumulate(rmScores.begin(), rmScores.end(), 0.0) / rmScores.size();
    double rm_se = (rmScores.size() >= 2) ? std::sqrt(std::accumulate(rmScores.begin(), rmScores.end(), 0.0,
                                  [rm_mean](double acc, double x){ return acc + (x - rm_mean)*(x - rm_mean); }) / rmScores.size()) : NA_REAL;
    
    // Fill result vectors
    t_hat[j] = tObs;
    score[j] = scoreObs;
    rmMean[j] = rm_mean;
    rmSE[j] = rm_se;
    bestScore[j] = scoreBest;
    mlci[j] = (scoreBest - rm_mean == 0.0) ? NA_REAL : (scoreObs - rm_mean) / (scoreBest - rm_mean);
    nResampleVec[j] = nResample;
  }
  
  return DataFrame::create(
    Named("t_hat") = t_hat,
    Named("score") = score,
    Named("rmMean") = rmMean,
    Named("rmSE") = rmSE,
    Named("bestScore") = bestScore,
    Named("mlci") = mlci,
    Named("nResample") = nResampleVec
  );
}
