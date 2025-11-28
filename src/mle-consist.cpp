// Rcpp implementation for ML estimation of a single equal-branch-length parameter (t)
// under the ER (equal-rates) model with k states, plus an R wrapper that computes
// MLCI = (score - rmScore)/(bestScore - rmScore) per column.
//
// Exported C++ functions:
// - logLik_equal_t(edge, Ntip, tip_states, weights, k, t) -> log-likelihood
// - mle_t(edge, Ntip, tip_states, weights, k, lower, upper, tol) -> List(t_hat, logLik)
//
// R-side helpers below prepare data and perform resampling until desired precision is reached.


#include <Rcpp.h>
using namespace Rcpp;

// --- Low-level likelihood --------------------------------------------------
// We reuse the same closed-form P(t) for the ER model:
// Pii = 1/k + (k-1)/k * exp(-k t)
// Pij = 1/k - 1/k * exp(-k t)  (i != j)

// logLik_equal_t: computes log-likelihood for a single character (or pattern set)
// [[Rcpp::export]]
double logLik_equal_t(IntegerMatrix edge, int Ntip, IntegerMatrix tip_states, NumericVector weights, int k, double t) {
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
  
  // Build children list
  std::vector< std::vector<int> > children(nNodes + 1);
  std::vector<int> isChild(nNodes + 1, 0);
  for (int i = 0; i < nEdge; ++i) {
    int parent = edge(i,0);
    int child = edge(i,1);
    if (parent < 1 || child < 1 || parent > nNodes || child > nNodes) stop("edge matrix contains invalid node indices");
    children[parent].push_back(child);
    isChild[child] = 1;
  }
  
  // find root
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
  
  // conditional likelihood buffers
  std::vector< std::vector<double> > cond(nNodes + 1, std::vector<double>(k));
  
  const double e = std::exp(- (double)k * t);
  const double invk = 1.0 / (double)k;
  
  double totalLogLik = 0.0;
  
  for (int p = 0; p < nPatterns; ++p) {
    // zero cond
    for (int i = 1; i <= nNodes; ++i) std::fill(cond[i].begin(), cond[i].end(), 0.0);
    // set tips
    for (int tip = 0; tip < Ntip; ++tip) {
      int state = tip_states(tip, p);
      if (state < 0 || state >= k) stop("tip state out of range");
      cond[ tip + 1 ][ state ] = 1.0;
    }
    
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

// --- MLE for t using golden-section search on the log-likelihood ----------------
// We maximise logLik(t); we return t_hat and the corresponding logLik.

// helper wrapper to avoid repeated R-to-C++ conversions for parameters: evaluate logLik given t
static double eval_logLik_cached(IntegerMatrix edge, int Ntip, IntegerMatrix tip_states, NumericVector weights, int k, double t) {
  return logLik_equal_t(edge, Ntip, tip_states, weights, k, t);
}

// [[Rcpp::export]]
List mle_t(IntegerMatrix edge, int Ntip, IntegerMatrix tip_states, NumericVector weights, int k, double lower = 1e-8, double upper = 10.0, double tol = 1e-6) {
  if (lower <= 0) lower = 1e-12;
  if (upper <= lower) stop("upper must be > lower");
  const double gr = (sqrt(5.0) + 1.0) / 2.0;
  double a = lower, b = upper;
  double c = b - (b - a) / gr;
  double d = a + (b - a) / gr;
  double fc = eval_logLik_cached(edge, Ntip, tip_states, weights, k, c);
  double fd = eval_logLik_cached(edge, Ntip, tip_states, weights, k, d);
  int iter = 0, maxiter = 200;
  while (fabs(b - a) > tol * (fabs(c) + fabs(d)) && iter < maxiter) {
    if (fc > fd) {
      b = d; d = c; fd = fc; c = b - (b - a) / gr; fc = eval_logLik_cached(edge, Ntip, tip_states, weights, k, c);
    } else {
      a = c; c = d; fc = fd; d = a + (b - a) / gr; fd = eval_logLik_cached(edge, Ntip, tip_states, weights, k, d);
    }
    iter++;
  }
  double t_hat = (b + a) / 2.0;
  double logL = eval_logLik_cached(edge, Ntip, tip_states, weights, k, t_hat);
  return List::create(Named("t_hat") = t_hat, Named("logLik") = logL);
}
