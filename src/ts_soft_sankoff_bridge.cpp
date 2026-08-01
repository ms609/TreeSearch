// R bridge for the soft-Sankoff prototype kernel.
//
// EXPLORATION CODE.  See dev/plans/2026-08-01-soft-sankoff-temperature-dial.md.
// Kept in its own translation unit so that ts_rcpp.cpp, which serves the live
// search path, is untouched by this branch.
//
// Tip costs are taken as explicit per-character (n_tip x n_states) cost
// matrices rather than as integer state indices.  That mirrors the verified R
// reference (`TipCosts()` in tests/testthat/helper-soft-sankoff.R) exactly, and
// side-steps the 0-based-tip-state trap that ts_sankoff_test() carries: an
// out-of-range index there does not error, it leaves every state at infinity
// and the whole score returns Inf.

#include <Rcpp.h>
#include <algorithm>
#include "ts_soft_sankoff.h"

using namespace Rcpp;

// Build left/right child arrays and an internal-node postorder from an ape
// edge matrix, validating that the tree really is rooted and binary.
//
// This is a deliberately stricter sibling of build_topo_from_edge() in
// ts_rcpp.cpp, which fills left/right from the first two edges it sees and
// silently drops any third child.  The R reference handles arbitrary arity, so
// a polytomy passed to a permissive bridge would give a silent mismatch rather
// than an error.
static void soft_build_topo(
    const int* edge_parent, const int* edge_child, int n_edge,
    int n_tip,
    std::vector<int>& left_out, std::vector<int>& right_out,
    std::vector<int>& postorder_out)
{
  const int n_internal = n_tip - 1;
  const int n_node = n_tip + n_internal;

  if (n_edge != 2 * n_internal) {
    stop("edge matrix has %d edges; a rooted binary tree on %d tips has %d",
         n_edge, n_tip, 2 * n_internal);
  }

  left_out.assign(n_internal, -1);
  right_out.assign(n_internal, -1);
  std::vector<int> n_kids(n_internal, 0);
  std::vector<bool> is_child(n_node, false);

  for (int i = 0; i < n_edge; ++i) {
    const int p = edge_parent[i] - 1;
    const int c = edge_child[i] - 1;
    if (p < n_tip || p >= n_node) {
      stop("edge %d has parent %d, outside the internal-node range %d..%d",
           i + 1, p + 1, n_tip + 1, n_node);
    }
    if (c < 0 || c >= n_node) {
      stop("edge %d has child %d, outside the node range 1..%d",
           i + 1, c + 1, n_node);
    }
    if (is_child[c]) {
      stop("node %d appears as a child of more than one edge", c + 1);
    }
    is_child[c] = true;

    const int pi = p - n_tip;
    if (n_kids[pi] == 0)      left_out[pi]  = c;
    else if (n_kids[pi] == 1) right_out[pi] = c;
    else {
      stop("node %d has more than two children; this kernel is binary-only",
           p + 1);
    }
    ++n_kids[pi];
  }

  for (int ni = 0; ni < n_internal; ++ni) {
    if (n_kids[ni] != 2) {
      stop("internal node %d has %d children, expected 2",
           ni + n_tip + 1, n_kids[ni]);
    }
  }

  // The root is the one node that is nobody's child.
  int root = -1;
  for (int v = 0; v < n_node; ++v) {
    if (!is_child[v]) {
      if (root >= 0) stop("edge matrix describes more than one root");
      root = v;
    }
  }
  if (root < n_tip) stop("the root of the edge matrix is a tip");

  // Two-stack postorder over internal nodes, leaves-to-root.
  postorder_out.clear();
  postorder_out.reserve(n_internal);
  std::vector<int> stk;
  stk.push_back(root);
  while (!stk.empty()) {
    const int nd = stk.back(); stk.pop_back();
    if (nd >= n_tip) {
      postorder_out.push_back(nd);
      const int ni = nd - n_tip;
      if (left_out[ni]  >= 0) stk.push_back(left_out[ni]);
      if (right_out[ni] >= 0) stk.push_back(right_out[ni]);
    }
  }
  if (static_cast<int>(postorder_out.size()) != n_internal) {
    stop("edge matrix is disconnected: %d of %d internal nodes reachable "
         "from the root", static_cast<int>(postorder_out.size()), n_internal);
  }
  std::reverse(postorder_out.begin(), postorder_out.end());
}

//' Soft-Sankoff score (exploration prototype)
//'
//' Sankoff's dynamic program with `min` replaced by a soft-minimum at
//' temperature `temperature`.  Not part of any exported scoring path.
//'
//' @param edge Two-column integer edge matrix, `ape` convention, rooted binary.
//' @param n_tip Number of tips.
//' @param tip_costs List of length `n_chars`; element `ch` is an
//'   `n_tip x n_states[ch]` numeric matrix of per-state tip costs (0 for an
//'   observed state, `Inf` otherwise).
//' @param cost_matrices List of length `n_chars`; element `ch` is an
//'   `n_states x n_states` matrix, `cost[from, to]`, applied to every edge.
//'   Ignored for a character that supplies `branch_costs`.
//' @param temperature Soft-min temperature; `0` gives hard weighted parsimony.
//' @param root_costs Optional list of length `n_chars` of additive per-state
//'   root costs (supply `-log(pi)` for root frequencies `pi`), or `NULL`.
//' @param branch_costs Optional list of length `n_chars`; element `ch` is
//'   either `NULL` or a list indexed by *child node* (`1..n_node`) of
//'   per-branch cost matrices, matching the R reference's `cost[[child]]`
//'   convention.  The root's element is never read.
//' @param n_rep Number of times to repeat the whole scoring pass.  Only the
//'   last pass's result is returned; the repetition exists so that timing
//'   measures the kernel rather than R-side marshalling.
//'
//' @return List with `score` (total) and `per_char`.
//' @keywords internal
// [[Rcpp::export]]
List ts_soft_sankoff_test(
    IntegerMatrix edge,
    int n_tip,
    List tip_costs,
    List cost_matrices,
    double temperature,
    Nullable<List> root_costs = R_NilValue,
    Nullable<List> branch_costs = R_NilValue,
    int n_rep = 1)
{
  const int n_edge = edge.nrow();
  const int n_chars = tip_costs.size();
  const int n_internal = n_tip - 1;
  const int n_node = n_tip + n_internal;
  const double INF = std::numeric_limits<double>::infinity();

  if (edge.ncol() != 2) stop("edge must have two columns");
  if (n_tip < 1) stop("n_tip must be at least 1");
  if (n_chars < 1) stop("tip_costs must hold at least one character");
  if (cost_matrices.size() != n_chars) {
    stop("cost_matrices has %d elements, tip_costs has %d",
         cost_matrices.size(), n_chars);
  }
  if (n_rep < 1) stop("n_rep must be at least 1");
  if (!(temperature >= 0.0)) {
    stop("temperature must be non-negative (0 = hard parsimony)");
  }

  std::vector<int> left_v, right_v, postorder;
  soft_build_topo(&edge(0, 0), &edge(0, 1), n_edge, n_tip,
                  left_v, right_v, postorder);

  const bool have_root = root_costs.isNotNull();
  const bool have_branch = branch_costs.isNotNull();
  List root_list, branch_list;
  if (have_root) {
    root_list = List(root_costs.get());
    if (root_list.size() != n_chars) {
      stop("root_costs has %d elements, expected %d", root_list.size(), n_chars);
    }
  }
  if (have_branch) {
    branch_list = List(branch_costs.get());
    if (branch_list.size() != n_chars) {
      stop("branch_costs has %d elements, expected %d",
           branch_list.size(), n_chars);
    }
  }

  // ---- Build SoftSankoffData ----
  ts::SoftSankoffData sd;
  sd.n_tips = n_tip;
  sd.n_chars = n_chars;
  sd.max_states = 0;
  sd.chars.resize(n_chars);

  std::vector<int> ns_vec(n_chars);
  for (int ch = 0; ch < n_chars; ++ch) {
    NumericMatrix tc = as<NumericMatrix>(tip_costs[ch]);
    if (tc.nrow() != n_tip) {
      stop("tip_costs[[%d]] has %d rows, expected n_tip = %d",
           ch + 1, tc.nrow(), n_tip);
    }
    const int ns = tc.ncol();
    if (ns < 1) stop("tip_costs[[%d]] has no state columns", ch + 1);
    ns_vec[ch] = ns;
    sd.chars[ch].n_states = ns;
    if (ns > sd.max_states) sd.max_states = ns;

    NumericMatrix cm = as<NumericMatrix>(cost_matrices[ch]);
    if (cm.nrow() != ns || cm.ncol() != ns) {
      stop("cost_matrices[[%d]] is %d x %d but character %d has %d states",
           ch + 1, cm.nrow(), cm.ncol(), ch + 1, ns);
    }
    sd.chars[ch].cost_matrix.resize(static_cast<size_t>(ns) * ns);
    for (int r = 0; r < ns; ++r) {
      for (int c = 0; c < ns; ++c) {
        sd.chars[ch].cost_matrix[static_cast<size_t>(r) * ns + c] = cm(r, c);
      }
    }

    if (have_root && !Rf_isNull(root_list[ch])) {
      NumericVector rc = as<NumericVector>(root_list[ch]);
      if (rc.size() != ns) {
        stop("root_costs[[%d]] has length %d, expected %d",
             ch + 1, rc.size(), ns);
      }
      sd.chars[ch].root_cost.assign(rc.begin(), rc.end());
    }

    if (have_branch && !Rf_isNull(branch_list[ch])) {
      List bl = as<List>(branch_list[ch]);
      if (bl.size() != n_node) {
        stop("branch_costs[[%d]] has %d elements; it must be indexed by child "
             "node, so length n_node = %d", ch + 1, bl.size(), n_node);
      }
      const size_t cm_size = static_cast<size_t>(ns) * ns;
      // Fill with infinity, so that a child node whose matrix was omitted
      // produces Inf rather than silently scoring as zero cost.
      sd.chars[ch].branch_costs.assign(
          static_cast<size_t>(n_node) * cm_size, INF);
      for (int v = 0; v < n_node; ++v) {
        if (Rf_isNull(bl[v])) continue;      // the root legitimately has none
        NumericMatrix bm = as<NumericMatrix>(bl[v]);
        if (bm.nrow() != ns || bm.ncol() != ns) {
          stop("branch_costs[[%d]][[%d]] is %d x %d, expected %d x %d",
               ch + 1, v + 1, bm.nrow(), bm.ncol(), ns, ns);
        }
        double* dst = sd.chars[ch].branch_costs.data()
            + static_cast<size_t>(v) * cm_size;
        for (int r = 0; r < ns; ++r) {
          for (int c = 0; c < ns; ++c) {
            dst[static_cast<size_t>(r) * ns + c] = bm(r, c);
          }
        }
      }
    }
  }

  // ---- Tip costs, strided to match ts::SoftSankoffData ----
  const int stride = sd.stride();
  sd.tip_costs.assign(static_cast<size_t>(n_tip) * stride, INF);
  for (int ch = 0; ch < n_chars; ++ch) {
    NumericMatrix tc = as<NumericMatrix>(tip_costs[ch]);
    for (int t = 0; t < n_tip; ++t) {
      double* dst = sd.tip_costs.data()
          + static_cast<size_t>(t) * stride
          + static_cast<size_t>(ch) * sd.max_states;
      for (int s = 0; s < ns_vec[ch]; ++s) {
        dst[s] = tc(t, s);
      }
    }
  }

  // ---- Score ----
  double total = 0.0;
  for (int rep = 0; rep < n_rep; ++rep) {
    total = ts::soft_sankoff_score(
        left_v.data(), right_v.data(),
        postorder.data(), n_internal, n_tip, sd, temperature);
  }

  NumericVector per_char(n_chars);
  for (int ch = 0; ch < n_chars; ++ch) {
    const double* ch_tip = sd.tip_costs.data() + ch * sd.max_states;
    per_char[ch] = ts::soft_sankoff_score_char(
        left_v.data(), right_v.data(),
        postorder.data(), n_internal, n_tip,
        sd.chars[ch], ch_tip, stride, temperature);
  }

  return List::create(
    Named("score") = total,
    Named("per_char") = per_char);
}
