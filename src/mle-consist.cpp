// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <R.h>        // for GetRNGstate / PutRNGstate / unif_rand
#include <Rmath.h>    // sometimes helps; safe to include
using namespace Rcpp;

// ---------- low-level ER log-likelihood (same as earlier) --------------
// P(t) closed form used. tip_states: Ntip x nPatterns, values 0..k-1

double logLik_equal_t_internal(const IntegerMatrix &edge, int Ntip,
                               const IntegerMatrix &tip_states,
                               const NumericVector &weights, int k, double t) {
  int nEdge = edge.nrow();
  int nNodes = 0;
  for (int i = 0; i < nEdge; ++i) {
    nNodes = std::max(nNodes, (int)edge(i,0));
    nNodes = std::max(nNodes, (int)edge(i,1));
  }
  int nPatterns = tip_states.ncol();
  
  if (k <= 0) stop("k must be > 0");
  if (Ntip <= 0) stop("Ntip must be > 0");
  if (tip_states.nrow() != Ntip) stop("tip_states must have Ntip rows");
  if ((int)weights.size() != nPatterns) stop("weights length must equal number of patterns");
  
  std::vector< std::vector<int> > children(nNodes + 1);
  std::vector<int> isChild(nNodes + 1, 0);
  for (int i = 0; i < nEdge; ++i) {
    int parent = edge(i,0);
    int child  = edge(i,1);
    if (parent < 1 || child < 1 || parent > nNodes || child > nNodes) stop("edge matrix contains invalid node indices");
    children[parent].push_back(child);
    isChild[child] = 1;
  }
  
  int root = -1;
  for (int v = 1; v <= nNodes; ++v) if (!isChild[v]) { root = v; break; }
  if (root == -1) stop("could not determine root from edge matrix");
  
  // postorder
  std::vector<int> postorder; postorder.reserve(nNodes);
  {
    std::vector<int> stack; stack.reserve(nNodes);
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
  
  // allocate cond buffers once
  std::vector< std::vector<double> > cond(nNodes + 1, std::vector<double>(k));
  const double e = std::exp(- (double)k * t);
  const double invk = 1.0 / (double)k;
  double totalLogLik = 0.0;
  
  for (int p = 0; p < nPatterns; ++p) {
    // zero
    for (int i = 1; i <= nNodes; ++i) std::fill(cond[i].begin(), cond[i].end(), 0.0);
    // set tips
    for (int tip = 0; tip < Ntip; ++tip) {
      int state = tip_states(tip, p);
      if (state < 0 || state >= k) stop("tip state out of range");
      cond[ tip + 1 ][ state ] = 1.0;
    }
    // postorder combine
    for (int idx = 0; idx < (int)postorder.size(); ++idx) {
      int node = postorder[idx];
      if (node <= Ntip) continue;
      std::vector<double> accum(k, 1.0);
      for (int child : children[node]) {
        double sumChild = 0.0;
        for (int s = 0; s < k; ++s) sumChild += cond[child][s];
        double meanChild = sumChild * invk;
        for (int s = 0; s < k; ++s) {
          double pts = meanChild + e * (cond[child][s] - meanChild);
          accum[s] *= pts;
        }
      }
      for (int s = 0; s < k; ++s) cond[node][s] = accum[s];
    }
    // root
    double sumRoot = 0.0;
    for (int s = 0; s < k; ++s) sumRoot += cond[root][s];
    double siteLik = sumRoot * invk;
    if (siteLik <= 0) {
      totalLogLik += weights[p] * log(std::numeric_limits<double>::min());
    } else {
      totalLogLik += weights[p] * log(siteLik);
    }
  }
  return totalLogLik;
}

// Exported convenience wrapper, identical behaviour as earlier function name
// [[Rcpp::export]]
double logLik_equal_t(IntegerMatrix edge, int Ntip, IntegerMatrix tip_states, NumericVector weights, int k, double t) {
  return logLik_equal_t_internal(edge, Ntip, tip_states, weights, k, t);
}

// ---------- internal golden-section maximiser for logLik (maximises logLik) -----------
static double mle_t_internal(IntegerMatrix &edge, int Ntip, IntegerMatrix &tip_states, NumericVector &weights, int k, double lower, double upper, double tol) {
  if (lower <= 0) lower = 1e-12;
  if (upper <= lower) stop("upper must be > lower");
  const double gr = (sqrt(5.0) + 1.0) / 2.0;
  double a = lower, b = upper;
  double c = b - (b - a) / gr;
  double d = a + (b - a) / gr;
  double fc = logLik_equal_t_internal(edge, Ntip, tip_states, weights, k, c);
  double fd = logLik_equal_t_internal(edge, Ntip, tip_states, weights, k, d);
  int iter = 0, maxiter = 200;
  while (fabs(b - a) > tol * (fabs(c) + fabs(d)) && iter < maxiter) {
    if (fc > fd) {
      b = d; d = c; fd = fc; c = b - (b - a) / gr; fc = logLik_equal_t_internal(edge, Ntip, tip_states, weights, k, c);
    } else {
      a = c; c = d; fc = fd; d = a + (b - a) / gr; fd = logLik_equal_t_internal(edge, Ntip, tip_states, weights, k, d);
    }
    iter++;
  }
  double t_hat = (b + a) / 2.0;
  return t_hat;
}

// exported small wrapper to match earlier name (optional)
// [[Rcpp::export]]
List mle_t(IntegerMatrix edge, int Ntip, IntegerMatrix tip_states, NumericVector weights, int k, double lower = 1e-8, double upper = 10.0, double tol = 1e-6) {
  double t_hat = mle_t_internal(edge, Ntip, tip_states, weights, k, lower, upper, tol);
  double logL = logLik_equal_t_internal(edge, Ntip, tip_states, weights, k, t_hat);
  return List::create(Named("t_hat") = t_hat, Named("logLik") = logL);
}

// ---------- C++ resampling MLCI routine -------------------------------------
// obs_states: IntegerVector length Ntip with values 0..k-1 (observed tokens coded as integers)
// best_split: IntegerVector length Ntip with values 0 or 1 (TRUE state -> 1, FALSE -> 0)
// precision: desired SE on MLCI, maxResample: cap
//
// Returns a list with t_hat_obs, score, t_hat_best, bestScore,
// rmMean (mean rmScore), rmSE (SE of rmMean), mlci, mlciSE, nResample

// [[Rcpp::export]]
List mlci_resample(IntegerMatrix edge, int Ntip,
                   IntegerVector obs_states, IntegerVector best_split,
                   double precision = 1e-3, int maxResample = 10000,
                   double lower = 1e-8, double upper = 10.0, double tol = 1e-6) {
  
  constexpr int minResample = 3;
  
  if ((int)obs_states.size() != Ntip) stop("obs_states must have length Ntip");
  if ((int)best_split.size() != Ntip) stop("best_split must have length Ntip");
  
  int nEdge = edge.nrow();
  
  // infer k from obs_states unique values
  std::vector<int> present; // map values -> seen
  int maxval = -1;
  for (int i = 0; i < Ntip; ++i) if (obs_states[i] > maxval) maxval = obs_states[i];
  if (maxval < 0) stop("obs_states must be non-negative integers");
  int possible_k = maxval + 1;
  // build a boolean seen array
  std::vector<char> seen(possible_k, 0);
  int k = 0;
  for (int i = 0; i < Ntip; ++i) {
    int v = obs_states[i];
    if (v < 0) stop("obs_states must be non-negative integers");
    if (v >= (int)seen.size()) seen.resize(v+1, 0);
    if (!seen[v]) { seen[v] = 1; ++k; }
  }
  
  if (k <= 0) stop("no states found");
  
  // Build tip_states matrix Ntip x 1 for observed
  IntegerMatrix tip_obs(Ntip, 1);
  for (int i = 0; i < Ntip; ++i) tip_obs(i,0) = obs_states[i];
  NumericVector w1(1); w1[0] = 1.0;
  
  // observed MLE
  double t_hat_obs = mle_t_internal(edge, Ntip, tip_obs, w1, k, lower, upper, tol);
  double score = t_hat_obs * nEdge;
  
  // best character: create tip_states with 0/1 states per best_split.
  IntegerMatrix tip_best(Ntip, 1);
  for (int i = 0; i < Ntip; ++i) tip_best(i,0) = best_split[i] ? 1 : 0;
  NumericVector w2(1); w2[0] = 1.0;
  double t_hat_best = mle_t_internal(edge, Ntip, tip_best, w2, 2, lower, upper, tol);
  double bestScore = t_hat_best * nEdge;
  
  // Pre-allocated buffers for use inside loop
  IntegerVector perm(Ntip);
  for (int i = 0; i < Ntip; ++i) perm[i] = obs_states[i];
  
  // Welford accumulators
  long long n = 0;
  double mean = 0.0;
  double M2 = 0.0;
  
  // prepare a tip_states matrix reused
  IntegerMatrix tip_r(Ntip, 1);
  NumericVector wr(1); wr[0] = 1.0;
  
  // RNG state
  GetRNGstate();
  
  // Fisher-Yates shuffle function inline
  auto shuffle_inplace = [&](IntegerVector &v) {
    int m = v.size();
    for (int i = m - 1; i > 0; --i) {
      double u = unif_rand(); // in [0,1)
      int j = (int)(u * (i + 1)); // 0..i
      if (j != i) {
        int tmp = v[i];
        v[i] = v[j];
        v[j] = tmp;
      }
    }
  };
  
  double rmScore = 0.0;
  double rmMean = NA_REAL;
  double rmSE = NA_REAL;
  double mlci = NA_REAL;
  double mlciSE = NA_REAL;
  
  // If denom zero edge-case: bestScore == score, then mlci is NaN/infinite; handle later.
  // Loop
  while (true) {
    // shuffle perm (it contains integer states)
    shuffle_inplace(perm);
    // build tip_r from perm
    for (int i = 0; i < Ntip; ++i) tip_r(i,0) = perm[i];
    // compute mle t for this resampled character
    // infer k_r for this permutation (some states may vanish) -> compute unique count quickly
    int maxv = -1;
    for (int i=0;i<Ntip;++i) if (perm[i] > maxv) maxv = perm[i];
    std::vector<char> seen_r(maxv+1, 0);
    int k_r = 0;
    for (int i=0;i<Ntip;++i) {
      int v = perm[i];
      if (v >= (int)seen_r.size()) seen_r.resize(v+1,0);
      if (!seen_r[v]) { seen_r[v] = 1; ++k_r; }
    }
    if (k_r <= 0) k_r = 1; // unlikely
    
    double t_hat_r = mle_t_internal(edge, Ntip, tip_r, wr, k_r, lower, upper, tol);
    rmScore = t_hat_r * nEdge;
    
    // update Welford
    n += 1;
    double delta = rmScore - mean;
    mean += delta / (double)n;
    double delta2 = rmScore - mean;
    M2 += delta * delta2;
    
    if (n >= minResample) {
      double var = M2 / (double)(n - 1);           // sample variance of rmScores
      double se_mean = std::sqrt(var) / std::sqrt((double)n); // SE of mean
      // propagate to SE of mlci via delta method
      double denom = (bestScore - mean);
      if (denom == 0.0) {
        mlciSE = R_PosInf;
      } else {
        double gprime = (score - bestScore) / (denom * denom);
        mlciSE = std::fabs(gprime) * se_mean;
      }
      // current mlci estimate
      mlci = (score - mean) / (bestScore - mean);
      rmMean = mean;
      rmSE = se_mean;
      // stopping condition:
      if (mlciSE <= precision) break;
      if (n >= maxResample) break;
    } // else need at least 2 samples to compute SE
  }
  
  PutRNGstate();
  
  // If only 1 sample: compute current rmMean = mean, rmSE = NA, mlciSE = NA
  if (n < 2) {
    rmMean = mean;
    rmSE = NA_REAL;
    mlci = (bestScore - mean == 0.0) ? NA_REAL : (score - mean) / (bestScore - mean);
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
