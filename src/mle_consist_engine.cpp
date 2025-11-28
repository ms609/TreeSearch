// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <cstring> // memset
using namespace Rcpp;

// C++ engine for ER-model likelihood with equal branch length t
// Preallocates buffers once; supports k up to maxK passed at construction.

struct ERLikelihoodEngine {
  int nEdge;
  int nNodes;
  int Ntip;
  int root;
  int maxK;
  
  // flattened tree representation
  std::vector<int> postorder;      // postorder nodes
  std::vector<int> children_flat;  // flattened children list
  std::vector<int> child_offsets;  // offsets into children_flat (1..nNodes)
  std::vector<int> child_counts;   // counts per node
  
  // single contiguous conditional likelihood buffer: (nNodes+1) * maxK
  std::vector<double> cond_buf;
  
  ERLikelihoodEngine(const IntegerMatrix &edge, int Ntip_, int maxK_) :
    Ntip(Ntip_), maxK(maxK_)
  {
    nEdge = edge.nrow();
    // find nNodes
    nNodes = 0;
    for (int i = 0; i < nEdge; ++i) {
      nNodes = std::max(nNodes, (int)edge(i,0));
      nNodes = std::max(nNodes, (int)edge(i,1));
    }
    
    // temporary children lists
    std::vector< std::vector<int> > children(nNodes + 1);
    std::vector<int> isChild(nNodes + 1, 0);
    for (int i = 0; i < nEdge; ++i) {
      int p = edge(i,0), c = edge(i,1);
      children[p].push_back(c);
      isChild[c] = 1;
    }
    
    // root
    root = -1;
    for (int v = 1; v <= nNodes; ++v) if (!isChild[v]) { root = v; break; }
    if (root == -1) Rf_error("could not determine root");
    
    // postorder
    {
      postorder.reserve(nNodes);
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
    
    // flatten children
    child_counts.assign(nNodes + 1, 0);
    for (int v = 1; v <= nNodes; ++v) child_counts[v] = (int)children[v].size();
    
    child_offsets.assign(nNodes + 1, 0);
    int offs = 0;
    for (int v = 1; v <= nNodes; ++v) {
      child_offsets[v] = offs;
      offs += child_counts[v];
    }
    children_flat.assign(offs, 0);
    for (int v = 1; v <= nNodes; ++v) {
      int pos = child_offsets[v];
      for (size_t i = 0; i < children[v].size(); ++i) children_flat[pos + i] = children[v][i];
    }
    
    // allocate cond_buf once (nNodes+1) * maxK
    if (maxK <= 0) Rf_error("maxK must be > 0");
    cond_buf.assign((size_t)(nNodes + 1) * (size_t)maxK, 0.0);
  }
  
  // core log-likelihood using preallocated cond_buf; does not allocate memory
  // tip_states: IntegerMatrix Ntip x nPatterns with values 0..k-1
  double logLik(const IntegerMatrix &tip_states, const NumericVector &weights, int k, double t) {
    if (k <= 0) Rf_error("k must be > 0");
    if (k > maxK) Rf_error("k > maxK for this engine");
    if (tip_states.nrow() != Ntip) Rf_error("tip_states must have Ntip rows");
    int nPatterns = tip_states.ncol();
    
    const double e = std::exp(- (double)k * t);
    const double invk = 1.0 / (double)k;
    
    // pointers for speed
    double *cond = cond_buf.data();
    int stride = maxK; // use stride = maxK so nodes align; only first k entries used
    
    // buffer for accumulating per-node product
    std::vector<double> accum((size_t)k);
    
    double totalLogLik = 0.0;
    
    for (int p = 0; p < nPatterns; ++p) {
      // zero cond_buf (only the used portion)
      std::memset(cond, 0, (size_t)(nNodes + 1) * (size_t)maxK * sizeof(double));
      // set tips: tips are nodes 1..Ntip
      for (int tip = 0; tip < Ntip; ++tip) {
        int state = tip_states(tip, p);
        if (state < 0 || state >= k) Rf_error("tip state out of range");
        cond[(tip + 1) * stride + state] = 1.0;
      }
      // postorder combine
      for (size_t idx = 0; idx < postorder.size(); ++idx) {
        int node = postorder[idx];
        if (node <= Ntip) continue;
        // init accum to 1
        for (int s = 0; s < k; ++s) accum[s] = 1.0;
        int offs = child_offsets[node];
        int cnt  = child_counts[node];
        for (int ci = 0; ci < cnt; ++ci) {
          int child = children_flat[offs + ci];
          double sumChild = 0.0;
          double *child_ptr = &cond[child * stride];
          for (int s = 0; s < k; ++s) sumChild += child_ptr[s];
          double meanChild = sumChild * invk;
          for (int s = 0; s < k; ++s) {
            double pts = meanChild + e * (child_ptr[s] - meanChild);
            accum[s] *= pts;
          }
        }
        double *node_ptr = &cond[node * stride];
        for (int s = 0; s < k; ++s) node_ptr[s] = accum[s];
      }
      // root likelihood
      double sumRoot = 0.0;
      double *root_ptr = &cond[root * stride];
      for (int s = 0; s < k; ++s) sumRoot += root_ptr[s];
      double siteLik = sumRoot * invk;
      if (siteLik <= 0.0) {
        totalLogLik += weights[p] * std::log(std::numeric_limits<double>::min());
      } else {
        totalLogLik += weights[p] * std::log(siteLik);
      }
    }
    return totalLogLik;
  }
  
  // golden-section MLE that calls logLik (no allocations inside logLik)
  double mle_t_internal(const IntegerMatrix &tip_states, const NumericVector &weights, int k, double lower, double upper, double tol) {
    if (lower <= 0) lower = 1e-12;
    if (upper <= lower) Rf_error("upper must be > lower");
    const double gr = (std::sqrt(5.0) + 1.0) / 2.0;
    double a = lower, b = upper;
    double c = b - (b - a) / gr;
    double d = a + (b - a) / gr;
    double fc = logLik(tip_states, weights, k, c);
    double fd = logLik(tip_states, weights, k, d);
    int iter = 0, maxiter = 200;
    while (std::fabs(b - a) > tol * (std::fabs(c) + std::fabs(d)) && iter < maxiter) {
      if (fc > fd) {
        b = d; d = c; fd = fc; c = b - (b - a) / gr; fc = logLik(tip_states, weights, k, c);
      } else {
        a = c; c = d; fc = fd; d = a + (b - a) / gr; fd = logLik(tip_states, weights, k, d);
      }
      ++iter;
    }
    double t_hat = (b + a) / 2.0;
    return t_hat;
  }
};

// XPtr factory for engine
// [[Rcpp::export]]
SEXP mlci_make_engine(const IntegerMatrix &edge, int Ntip, int maxK = 10) {
  ERLikelihoodEngine *eng = new ERLikelihoodEngine(edge, Ntip, maxK);
  XPtr<ERLikelihoodEngine> p(eng, true);
  return p;
}

// engine helpers
// [[Rcpp::export]]
double cpp_engine_logLik(SEXP eng_xptr, const IntegerMatrix &tip_states, const NumericVector &weights, int k, double t) {
  XPtr<ERLikelihoodEngine> eng(eng_xptr);
  return eng->logLik(tip_states, weights, k, t);
}

// [[Rcpp::export]]
List cpp_engine_mle(SEXP eng_xptr, const IntegerMatrix &tip_states, const NumericVector &weights, int k, double lower = 1e-8, double upper = 10.0, double tol = 1e-6) {
  XPtr<ERLikelihoodEngine> eng(eng_xptr);
  double t_hat = eng->mle_t_internal(tip_states, weights, k, lower, upper, tol);
  double logL  = eng->logLik(tip_states, weights, k, t_hat);
  return List::create(Named("t_hat") = t_hat, Named("logLik") = logL);
}

// resampling MLCI using the engine (all in C++)
// [[Rcpp::export]]
List mlci_resample_engine(SEXP eng_xptr,
                          const IntegerVector &obs_states,
                          const IntegerVector &best_states,
                          double precision = 1e-3, int maxResample = 10000,
                          double lower = 1e-8, double upper = 10.0, double tol = 1e-6) {
  XPtr<ERLikelihoodEngine> eng(eng_xptr);
  int Ntip = eng->Ntip;
  int nEdge = eng->nEdge;
  
  if ((int)obs_states.size() != Ntip) Rf_error("obs_states length mismatch");
  if ((int)best_states.size() != Ntip) Rf_error("best_states length mismatch");
  
  // observed tip matrix (Ntip x 1)
  IntegerMatrix tip_obs(Ntip, 1);
  for (int i = 0; i < Ntip; ++i) tip_obs(i,0) = obs_states[i];
  NumericVector w1(1); w1[0] = 1.0;
  
  // infer k_obs
  int maxv = -1;
  for (int i = 0; i < Ntip; ++i) if (obs_states[i] > maxv) maxv = obs_states[i];
  int k_obs = maxv + 1;
  if (k_obs <= 0) Rf_error("k_obs invalid");
  
  // observed MLE and score
  double t_hat_obs = eng->mle_t_internal(tip_obs, w1, k_obs, lower, upper, tol);
  double score = t_hat_obs * nEdge;
  
  // best character MLE
  IntegerMatrix tip_best(Ntip, 1);
  for (int i = 0; i < Ntip; ++i) tip_best(i,0) = best_states[i];
  double t_hat_best = eng->mle_t_internal(tip_best, w1, 2, lower, upper, tol);
  double bestScore = t_hat_best * nEdge;
  
  // prepare permutation vector
  std::vector<int> perm((size_t)Ntip);
  for (int i = 0; i < Ntip; ++i) perm[i] = obs_states[i];
  
  long long n = 0;
  double mean_rm = 0.0, M2 = 0.0;
  IntegerMatrix tip_r(Ntip, 1);
  NumericVector wr(1); wr[0] = 1.0;
  
  GetRNGstate();
  auto shuffle_inplace = [&](std::vector<int> &v) {
    for (int i = (int)v.size() - 1; i > 0; --i) {
      double u = unif_rand();
      int j = (int)(u * (i + 1));
      int tmp = v[i]; v[i] = v[j]; v[j] = tmp;
    }
  };
  
  double rmScore = 0.0;
  double rmMean = NA_REAL, rmSE = NA_REAL, mlci = NA_REAL, mlciSE = NA_REAL;
  
  while (true) {
    shuffle_inplace(perm);
    for (int i = 0; i < Ntip; ++i) tip_r(i,0) = perm[i];
    // infer k_r
    int maxvv = -1;
    for (int i = 0; i < Ntip; ++i) if (perm[i] > maxvv) maxvv = perm[i];
    int k_r = maxvv + 1;
    if (k_r < 1) k_r = 1;
    
    double t_hat_r = eng->mle_t_internal(tip_r, wr, k_r, lower, upper, tol);
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
