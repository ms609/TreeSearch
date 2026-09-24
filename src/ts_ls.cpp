#include "ts_ls.h"
#include "ts_rng.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

namespace ts {

// ---------------------------------------------------------------------------
//  Target data
// ---------------------------------------------------------------------------

LSData build_ls_data(const double* dmat, int n_tip, const double* wmat) {
  LSData ls;
  ls.n_tip = n_tip;
  ls.n_pairs = n_tip * (n_tip - 1) / 2;
  ls.target.resize(ls.n_pairs);
  ls.sqrt_weight.assign(ls.n_pairs, 1.0);
  ls.weighted = (wmat != nullptr);

  int p = 0;
  for (int i = 0; i < n_tip; ++i) {
    for (int j = i + 1; j < n_tip; ++j, ++p) {
      // D is symmetric, so column-major and row-major indexing coincide.
      ls.target[p] = dmat[i + static_cast<size_t>(j) * n_tip];
      if (wmat) {
        double w = wmat[i + static_cast<size_t>(j) * n_tip];
        ls.sqrt_weight[p] = (w > 0.0) ? std::sqrt(w) : 0.0;
      }
    }
  }
  return ls;
}

// ---------------------------------------------------------------------------
//  Design (normal equations) assembled from the tree's unrooted branches
// ---------------------------------------------------------------------------
//
// We never materialise the full design matrix.  For each tip pair we find the
// branches on the path between the two tips (those whose split separates them)
// and accumulate the weighted normal-equation contributions
//   A[a][b] += w   for a,b on the path     (A = X^T W X)
//   g[a]    += w * D     for a on the path  (g = X^T W y)
//   yy      += w * D^2                      (= y^T W y)
// where w = weight of the pair.  Branch lengths then solve A v = g (OLS) or
// the non-negative variant (NNLS), and RSS = yy - 2 v.g + v.A.v.

struct Design {
  int n_branch = 0;
  int wps = 0;
  std::vector<double> A;            // n_branch * n_branch (row-major, symmetric)
  std::vector<double> g;            // n_branch
  double yy = 0.0;
  std::vector<int> branch_node;     // [n_branch] child node id of each branch
  std::vector<uint64_t> child_bits; // [n_branch * wps] tip set on the child side
};

static void build_design(const TreeState& tree, const LSData& ls, Design& d) {
  const int n_tip = tree.n_tip;
  const int wps = (n_tip + 63) / 64;
  d.wps = wps;

  // Per-node subtree tip bitsets (postorder union).
  std::vector<uint64_t> node_bits(static_cast<size_t>(tree.n_node) * wps, 0ULL);
  for (int t = 0; t < n_tip; ++t) {
    node_bits[static_cast<size_t>(t) * wps + (t / 64)] = 1ULL << (t % 64);
  }
  for (int node : tree.postorder) {
    int ni = node - n_tip;
    const uint64_t* lb = &node_bits[static_cast<size_t>(tree.left[ni]) * wps];
    const uint64_t* rb = &node_bits[static_cast<size_t>(tree.right[ni]) * wps];
    uint64_t* db = &node_bits[static_cast<size_t>(node) * wps];
    for (int w = 0; w < wps; ++w) db[w] = lb[w] | rb[w];
  }

  // Enumerate unrooted branches: every node except the root and the root's
  // right child (whose split duplicates the root's left child in the unrooted
  // tree).  Yields n_tip pendant + (n_tip - 3) internal = 2*n_tip - 3 columns.
  const int root = n_tip;
  const int root_right = tree.right[0];
  d.branch_node.clear();
  d.child_bits.clear();
  d.branch_node.reserve(2 * n_tip - 3);
  for (int node = 0; node < tree.n_node; ++node) {
    if (node == root || node == root_right) continue;
    d.branch_node.push_back(node);
    const uint64_t* nb = &node_bits[static_cast<size_t>(node) * wps];
    d.child_bits.insert(d.child_bits.end(), nb, nb + wps);
  }
  const int B = static_cast<int>(d.branch_node.size());
  d.n_branch = B;

  d.A.assign(static_cast<size_t>(B) * B, 0.0);
  d.g.assign(B, 0.0);
  d.yy = 0.0;

  // Accumulate over tip pairs.
  std::vector<int> sep;
  sep.reserve(B);
  int p = 0;
  for (int i = 0; i < n_tip; ++i) {
    const int iw = i / 64, ib = i % 64;
    for (int j = i + 1; j < n_tip; ++j, ++p) {
      const double sw = ls.sqrt_weight[p];
      if (sw == 0.0) continue;            // zero-weight pair drops out
      const int jw = j / 64, jb = j % 64;
      const double w = sw * sw;
      const double dij = ls.target[p];

      // Branches separating i and j: split contains exactly one of them.
      sep.clear();
      for (int b = 0; b < B; ++b) {
        const uint64_t* cb = &d.child_bits[static_cast<size_t>(b) * wps];
        bool has_i = (cb[iw] >> ib) & 1ULL;
        bool has_j = (cb[jw] >> jb) & 1ULL;
        if (has_i != has_j) sep.push_back(b);
      }

      d.yy += w * dij * dij;
      const double wd = w * dij;
      for (int a : sep) {
        d.g[a] += wd;
        double* Aa = &d.A[static_cast<size_t>(a) * B];
        for (int b : sep) Aa[b] += w;
      }
    }
  }
}

// ---------------------------------------------------------------------------
//  Linear algebra
// ---------------------------------------------------------------------------

// Cholesky solve of the SPD system A x = b, where A (n*n, row-major) is given
// by the index list `idx` into a larger symmetric matrix `Afull` (size NxN).
// Solves the idx-restricted subsystem.  Returns false if A is not PD.
static bool chol_solve_sub(const double* Afull, int N,
                           const int* idx, int n,
                           const double* bfull, double* x) {
  if (n == 0) return true;
  std::vector<double> L(static_cast<size_t>(n) * n, 0.0);
  // Scale-aware pivot floor.
  double maxdiag = 0.0;
  for (int i = 0; i < n; ++i) {
    double aii = Afull[static_cast<size_t>(idx[i]) * N + idx[i]];
    maxdiag = std::max(maxdiag, aii);
  }
  const double floor_piv = (maxdiag > 0.0 ? maxdiag : 1.0) * 1e-12;

  for (int j = 0; j < n; ++j) {
    double sum = Afull[static_cast<size_t>(idx[j]) * N + idx[j]];
    for (int k = 0; k < j; ++k) sum -= L[j * n + k] * L[j * n + k];
    if (sum <= floor_piv) return false;
    double Ljj = std::sqrt(sum);
    L[j * n + j] = Ljj;
    for (int i = j + 1; i < n; ++i) {
      double s = Afull[static_cast<size_t>(idx[i]) * N + idx[j]];
      for (int k = 0; k < j; ++k) s -= L[i * n + k] * L[j * n + k];
      L[i * n + j] = s / Ljj;
    }
  }
  // Forward solve L z = b
  std::vector<double> z(n);
  for (int i = 0; i < n; ++i) {
    double s = bfull[idx[i]];
    for (int k = 0; k < i; ++k) s -= L[i * n + k] * z[k];
    z[i] = s / L[i * n + i];
  }
  // Back solve L^T x = z
  for (int i = n - 1; i >= 0; --i) {
    double s = z[i];
    for (int k = i + 1; k < n; ++k) s -= L[k * n + i] * x[k];
    x[i] = s / L[i * n + i];
  }
  return true;
}

// Ordinary least squares: solve A v = g for all branches.
static bool solve_ols(const Design& d, std::vector<double>& v) {
  const int B = d.n_branch;
  v.assign(B, 0.0);
  std::vector<int> idx(B);
  for (int i = 0; i < B; ++i) idx[i] = i;
  return chol_solve_sub(d.A.data(), B, idx.data(), B, d.g.data(), v.data());
}

// Non-negative least squares (Lawson & Hanson active-set), operating on the
// precomputed normal equations (A, g).  Matches phangorn::nnls.tree().
static bool solve_nnls(const Design& d, std::vector<double>& v) {
  const int B = d.n_branch;
  v.assign(B, 0.0);
  if (B == 0) return true;

  std::vector<char> passive(B, 0);
  std::vector<double> w(B), s(B, 0.0);
  std::vector<int> idx;
  idx.reserve(B);

  double maxdiag = 0.0;
  for (int i = 0; i < B; ++i)
    maxdiag = std::max(maxdiag, d.A[static_cast<size_t>(i) * B + i]);
  const double scale = (maxdiag > 0.0 ? maxdiag : 1.0);
  const double tol = scale * 1e-11;

  const int max_outer = 3 * B + 10;
  for (int outer = 0; outer < max_outer; ++outer) {
    // gradient w = g - A v
    for (int a = 0; a < B; ++a) {
      double Av = 0.0;
      const double* Aa = &d.A[static_cast<size_t>(a) * B];
      for (int b = 0; b < B; ++b) Av += Aa[b] * v[b];
      w[a] = d.g[a] - Av;
    }
    // pick most-violated active coordinate
    int t = -1; double wmax = tol;
    for (int a = 0; a < B; ++a) {
      if (!passive[a] && w[a] > wmax) { wmax = w[a]; t = a; }
    }
    if (t < 0) break;
    passive[t] = 1;

    const int max_inner = 3 * B + 10;
    for (int inner = 0; inner < max_inner; ++inner) {
      idx.clear();
      for (int a = 0; a < B; ++a) if (passive[a]) idx.push_back(a);
      std::vector<double> zk(idx.size(), 0.0);
      if (!chol_solve_sub(d.A.data(), B, idx.data(),
                          static_cast<int>(idx.size()), d.g.data(),
                          zk.data())) {
        return false;   // passive submatrix singular — degenerate topology
      }
      std::fill(s.begin(), s.end(), 0.0);
      bool all_pos = true;
      for (size_t r = 0; r < idx.size(); ++r) {
        s[idx[r]] = zk[r];
        if (zk[r] <= tol) all_pos = false;
      }
      if (all_pos) { v = s; break; }

      // blocking step: largest alpha keeping v >= 0
      double alpha = std::numeric_limits<double>::infinity();
      for (int a : idx) {
        if (s[a] <= tol) {
          double denom = v[a] - s[a];
          if (denom > 0.0) alpha = std::min(alpha, v[a] / denom);
        }
      }
      if (!std::isfinite(alpha)) alpha = 0.0;
      for (int a = 0; a < B; ++a) v[a] += alpha * (s[a] - v[a]);
      for (int a = 0; a < B; ++a) {
        if (passive[a] && v[a] <= tol) { passive[a] = 0; v[a] = 0.0; }
      }
    }
  }
  return true;
}

// RSS = yy - 2 v.g + v.A.v  (clamped at 0 to absorb round-off).
static double residual_ss(const Design& d, const std::vector<double>& v) {
  const int B = d.n_branch;
  double vg = 0.0, vAv = 0.0;
  for (int a = 0; a < B; ++a) {
    vg += v[a] * d.g[a];
    double Av = 0.0;
    const double* Aa = &d.A[static_cast<size_t>(a) * B];
    for (int b = 0; b < B; ++b) Av += Aa[b] * v[b];
    vAv += v[a] * Av;
  }
  double rss = d.yy - 2.0 * vg + vAv;
  return rss > 0.0 ? rss : 0.0;
}

// ---------------------------------------------------------------------------
//  Public fitting / scoring entry points
// ---------------------------------------------------------------------------

LSFit ls_fit(const TreeState& tree, const LSData& ls, LSMethod method) {
  LSFit fit;
  Design d;
  build_design(tree, ls, d);
  fit.n_branch = d.n_branch;
  fit.branch_node = d.branch_node;

  std::vector<double> v;
  bool ok = (method == LSMethod::NNLS) ? solve_nnls(d, v) : solve_ols(d, v);
  if (!ok) {
    // Rank-deficient normal equations (e.g. enough zero-weight pairs to leave a
    // branch unidentifiable).  Signal failure, but still return a fully-sized,
    // finite length vector so callers never index past the end of
    // branch_length.  RSS is +Inf to mark the fit as unusable.
    fit.ok = false;
    fit.branch_length.assign(fit.n_branch, 0.0);
    fit.rss = std::numeric_limits<double>::infinity();
    return fit;
  }

  fit.ok = true;
  fit.branch_length = v;
  fit.rss = residual_ss(d, v);
  return fit;
}

double ls_score(const TreeState& tree, const LSData& ls, LSMethod method) {
  Design d;
  build_design(tree, ls, d);
  std::vector<double> v;
  bool ok = (method == LSMethod::NNLS) ? solve_nnls(d, v) : solve_ols(d, v);
  if (!ok) return std::numeric_limits<double>::infinity();
  return residual_ss(d, v);
}

// ---------------------------------------------------------------------------
//  Topology search (full LS rescore per candidate)
// ---------------------------------------------------------------------------

LSSearchResult ls_nni_search(TreeState& tree, const LSData& ls,
                             LSMethod method, int max_hits) {
  LSSearchResult res;
  double best = ls_score(tree, ls, method);
  res.rss = best;
  if (tree.n_tip < 4) return res;   // single unrooted topology

  std::vector<int> edges = tree.nni_edges();
  std::mt19937 rng = ts::make_rng();
  const double eps = 1e-9;
  int hits = 1;

  bool improved = true;
  while (improved) {
    improved = false;
    std::shuffle(edges.begin(), edges.end(), rng);

    for (int c : edges) {
      for (int which = 0; which < 2; ++which) {
        auto undo = tree.nni_apply(c, which);
        tree.build_postorder();
        ++res.n_iterations;
        double cand = ls_score(tree, ls, method);

        bool accept = false;
        if (cand < best - eps) {
          best = cand; hits = 1; accept = true;
        } else if (cand < best + eps && hits <= max_hits) {
          // equal score: accept up to max_hits to explore plateaus
          ++hits; accept = true;
        }

        if (accept) {
          res.rss = best;
          ++res.n_moves;
          improved = true;
          break;       // first-improvement: restart edge scan
        }
        tree.nni_undo(undo);
        tree.build_postorder();
      }
      if (improved) break;
      if (ts::check_interrupt()) { improved = false; goto done; }
    }
  }
done:
  tree.build_postorder();
  res.rss = ls_score(tree, ls, method);
  return res;
}

// SPR search: clip every movable subtree, try every regraft edge, full rescore.
LSSearchResult ls_spr_search(TreeState& tree, const LSData& ls,
                             LSMethod method, int max_hits) {
  LSSearchResult res;
  double best = ls_score(tree, ls, method);
  res.rss = best;
  if (tree.n_tip < 4) return res;

  std::mt19937 rng = ts::make_rng();
  const double eps = 1e-9;
  int hits = 1;

  std::vector<int> clip_candidates;
  for (int node = 0; node < tree.n_node; ++node) {
    if (node == tree.n_tip) continue;             // root
    clip_candidates.push_back(node);
  }

  std::vector<std::pair<int,int>> destinations;

  bool improved = true;
  while (improved) {
    improved = false;
    std::shuffle(clip_candidates.begin(), clip_candidates.end(), rng);

    for (int clip_node : clip_candidates) {
      if (tree.parent[clip_node] == tree.n_tip) continue;  // child of root

      tree.spr_clip(clip_node);
      tree.build_postorder();

      int ns = tree.clip_state.clip_sibling;
      int nz = tree.clip_state.clip_grandpar;

      // Collect destination edges in the divided main tree.
      destinations.clear();
      {
        std::vector<int> stack;
        stack.push_back(tree.n_tip);
        while (!stack.empty()) {
          int node = stack.back(); stack.pop_back();
          if (node < tree.n_tip) continue;
          int ni = node - tree.n_tip;
          int lc = tree.left[ni], rc = tree.right[ni];
          destinations.push_back({node, lc});
          destinations.push_back({node, rc});
          stack.push_back(lc);
          stack.push_back(rc);
        }
      }

      bool accepted = false;
      for (auto& [above, below] : destinations) {
        if (above == nz && below == ns) continue;   // original position
        tree.spr_regraft(above, below);
        tree.build_postorder();
        ++res.n_iterations;
        double cand = ls_score(tree, ls, method);

        bool accept = false;
        if (cand < best - eps) { best = cand; hits = 1; accept = true; }
        else if (cand < best + eps && hits <= max_hits) { ++hits; accept = true; }

        if (accept) {
          res.rss = best;
          ++res.n_moves;
          accepted = true;
          improved = true;
          break;
        }
        tree.spr_unregraft(above, below);
      }

      if (!accepted) {
        tree.spr_unclip();
        tree.build_postorder();
      }
      if (improved) break;
      if (ts::check_interrupt()) { improved = false; goto spr_done; }
    }
  }
spr_done:
  tree.build_postorder();
  res.rss = ls_score(tree, ls, method);
  return res;
}

} // namespace ts
