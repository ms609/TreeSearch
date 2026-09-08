#include "ts_hsj.h"
#include "ts_fitch.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <initializer_list>
#include <limits>

namespace ts {

std::vector<int> partition_weights(
    const int* index_r, int n_orig_chars,
    const int* weight_r, int n_patterns,
    const std::vector<int>& hierarchy_chars)
{
  std::vector<int> adjusted(weight_r, weight_r + n_patterns);

  for (int c : hierarchy_chars) {
    if (c < 0 || c >= n_orig_chars) continue;
    int pat = index_r[c];
    if (pat >= 0 && pat < n_patterns && adjusted[pat] > 0) {
      --adjusted[pat];
    }
  }

  return adjusted;
}

// A traversal of the tree rooted at tip 0, used for the secondary labelling
// only (T-374).
//
// The HSJ score is defined as a minimum over internal-node labelings of a sum
// of SYMMETRIC dissimilarities over the branches of an UNROOTED tree (Hopkins
// & St John 2021, p.3, p.6), so it may not depend on where the tree happens to
// be rooted.  The a(n)/p(n) DP below satisfies that already -- its branch
// costs are symmetric and it minimizes over the root's own state.  The
// secondary labelling did not: fitch_label_char() resolves ambiguous nodes
// with a DELTRAN-style uppass whose direction, and with subtree support counts
// whose subtrees, are both properties of the INPUT rooting.  Two rootings of
// one topology therefore disagreed about d(u, v) and so about the score.
//
// Rooting the labelling pass at tip 0 makes it a function of the unrooted
// topology and the data alone.  Tip indices come from the dataset, not from
// the rooting, so this is canonical.  Note this is NOT the reporting-boundary
// canonicalisation PR #278 applied to XFORM: this sits inside the kernel, so
// the objective the SEARCH optimizes is itself rooting-invariant, and
// MaximizeParsimony()'s reported score agrees with TreeLength() of the trees
// it returns by construction rather than by re-scoring at the boundary.
//
// Neighbour lists are sorted by node index so the traversal order depends only
// on the numbering, never on the incoming parent/child orientation.
struct CanonOrder {
  std::vector<int> post;     // postorder; canonical root (tip 0) last
  std::vector<int> kids;     // children, flattened
  std::vector<int> kidOff;   // kids[kidOff[n] .. kidOff[n] + kidNum[n])
  std::vector<int> kidNum;
};

static CanonOrder build_canon_order(const TreeState& tree) {
  const int n_tip = tree.n_tip;
  const int n_node = tree.n_node;

  // Undirected adjacency.  Every node has degree <= 3: the kernel's root is a
  // degree-2 subdivision point of an unrooted edge, and Fitch passes such a
  // node through transparently, so its position cannot bias the labelling.
  std::vector<int> adj(static_cast<size_t>(n_node) * 3, -1);
  std::vector<int> deg(n_node, 0);
  auto link = [&](int u, int v) {
    if (deg[u] < 3) adj[static_cast<size_t>(u) * 3 + deg[u]++] = v;
    if (deg[v] < 3) adj[static_cast<size_t>(v) * 3 + deg[v]++] = u;
  };
  for (int node : tree.postorder) {
    int ni = node - n_tip;
    link(node, tree.left[ni]);
    link(node, tree.right[ni]);
  }
  for (int n = 0; n < n_node; ++n) {
    std::sort(adj.begin() + static_cast<size_t>(n) * 3,
              adj.begin() + static_cast<size_t>(n) * 3 + deg[n]);
  }

  CanonOrder co;
  co.kidOff.assign(n_node, 0);
  co.kidNum.assign(n_node, 0);
  co.post.reserve(n_node);
  co.kids.reserve(n_node);

  // Iterative DFS from tip 0, emitting a preorder we then reverse.
  std::vector<int> pre;
  pre.reserve(n_node);
  std::vector<char> seen(n_node, 0);
  std::vector<int> stack;
  stack.push_back(0);
  seen[0] = 1;
  while (!stack.empty()) {
    int n = stack.back();
    stack.pop_back();
    pre.push_back(n);
    co.kidOff[n] = static_cast<int>(co.kids.size());
    for (int k = 0; k < deg[n]; ++k) {
      int nb = adj[static_cast<size_t>(n) * 3 + k];
      if (nb < 0 || seen[nb]) continue;
      seen[nb] = 1;
      co.kids.push_back(nb);
      ++co.kidNum[n];
      stack.push_back(nb);
    }
  }
  co.post.assign(pre.rbegin(), pre.rend());
  return co;
}

// Fitch downpass + uppass for a single character represented as integer state
// labels.  After the downpass, state sets at internal nodes are ambiguous
// (intersection or union of children).  The uppass resolves each node to a
// single state so that parent–child mismatches can be detected for HSJ
// secondary dissimilarity.
// Returns number of Fitch steps (union operations in the downpass).
//
// tip_labels holds 0-based TOKEN (allLevels/contrast-row) indices, not state
// indices (T-375): a token like "?" is its own row in the contrast matrix,
// generally with several columns set, so treating the token index itself as
// a bit position (as this function formerly did) is a category error --
// `label > 30` can never fire on a valid token index, so "?" silently scored
// as one concrete, arbitrary state instead of the wildcard it denotes.
// token_states[label] gives the actual bitmask of states the token is
// compatible with (populated by build_dataset() from the contrast matrix),
// which is what state_sets must hold to make the Fitch downpass/uppass below
// correct for ambiguous tokens.
//
// `pri_free[t]` marks the tips at which this secondary carries no constraint:
// those whose controlling primary CANNOT code the structure present, so the
// character does not exist there and its "-" is not a state it takes (T-374).
// score_hierarchy_block() computes it and documents why the test is that
// strict one rather than "may be absent".  Admitting "-" as an ordinary
// concrete state -- as this function formerly did, and as the comment here
// formerly asserted was deliberate -- let the uppass propagate it INWARDS and
// resolve a node in the middle of the PRESENT region to it, where it is
// disjoint from every present neighbour in every secondary at once, and
// score_hierarchy_block() charged that branch d = m, the full alpha, for a
// node that by construction has no inapplicable secondaries.  That over-charge
// is wrong under any rooting (the paper's d counts "nonmatching secondary
// characters", p.5, among characters that APPLY), and because whether it fired
// depended on the DELTRAN direction it was also the dominant source of
// T-374's rooting-dependence.
static int fitch_label_char(
    const TreeState& tree,
    const std::vector<int>& tip_labels,
    int char_idx,
    int n_orig_chars,
    const std::vector<uint32_t>& token_states,
    int n_levels,
    const std::vector<char>& pri_free,
    const CanonOrder& co,
    std::vector<uint32_t>& state_sets)
{
  int n_tip = tree.n_tip;
  int n_node = tree.n_node;

  // The applicable domain: the states this character is observed in at tips
  // where it actually applies.  Wildcarding to this rather than to all
  // n_levels bits keeps the tie-break arrays below as small as they were.
  uint32_t domain = 0;
  for (int t = 0; t < n_tip; ++t) {
    if (!pri_free[t]) {
      domain |= token_states[tip_labels[t * n_orig_chars + char_idx]];
    }
  }
  // The character applies nowhere: it constrains nothing.  Give every node one
  // shared state so no branch can ever register a mismatch.
  if (domain == 0) domain = 1u;

  uint32_t used_mask = 0;
  std::vector<uint32_t> observed(n_tip);
  for (int t = 0; t < n_tip; ++t) {
    int label = tip_labels[t * n_orig_chars + char_idx];
    observed[t] = pri_free[t] ? domain : token_states[label];
    state_sets[t] = observed[t];
    used_mask |= observed[t];
  }

  // --- Downpass, in the canonical (tip-0-rooted) postorder ---
  // Generalized over arity: the kernel's own root becomes an ordinary degree-2
  // node here, which Fitch passes through unchanged, and the canonical root is
  // tip 0, which has one child AND an observation of its own to honour.
  int steps = 0;
  for (int node : co.post) {
    int nk = co.kidNum[node];
    if (nk == 0) continue;                       // canonical leaf: tip state
    const int* kid = &co.kids[co.kidOff[node]];
    bool have = false;
    uint32_t inter = 0, uni = 0;
    if (node < n_tip) {                          // the canonical root
      inter = uni = observed[node];
      have = true;
    }
    for (int k = 0; k < nk; ++k) {
      uint32_t s = state_sets[kid[k]];
      if (!have) { inter = uni = s; have = true; }
      else { inter &= s; uni |= s; }
    }
    if (inter != 0) {
      state_sets[node] = inter;
    } else {
      state_sets[node] = uni;
      ++steps;
    }
  }

  // --- Order-invariant tie-break support -------------------------------
  // The uppass below must resolve ambiguous nodes to a single state so that
  // parent-child mismatches (the HSJ secondary dissimilarity) can be counted.
  // Resolving by bit index (the old "lowest set bit") would make the result
  // depend on the arbitrary phyDat `levels` ordering, since which state
  // occupies the lowest bit is determined by `levels`.  Instead we resolve
  // toward the state with the most support in the node's own subtree,
  // breaking ties by the smallest supporting tip index.  Both keys are
  // properties of the *states* and the tree, not of the bit encoding, so the
  // resolution — and hence the mismatch count — is invariant to level
  // ordering.  This still yields a valid most-parsimonious reconstruction, so
  // the dissimilarity stays non-zero (the concern that motivated adding the
  // uppass in the first place).
  //
  // tb_cnt[node * K + s]    = # tips in subtree(node) carrying concrete state s
  // tb_mintip[node * K + s] = smallest tip index in subtree(node) with state s
  //
  // K only needs to cover states this CHARACTER actually uses -- the downpass
  // above only ever intersects/unions existing tip bits, so no bit outside
  // `used_mask` (accumulated while initializing state_sets) can appear at any
  // node. Bounding K by `used_mask`'s highest set bit (rather than the
  // dataset-global n_levels) keeps these arrays as small as the old
  // per-character sizing did, even when other characters in the dataset use
  // many more states than this one.
  int K = 0;
  for (uint32_t m = used_mask; m; m >>= 1) ++K;
  if (K == 0) K = 1;
  if (K > n_levels) K = n_levels;  // defensive: never exceed the real state space
  const int INF_TIP = std::numeric_limits<int>::max();
  std::vector<int> tb_cnt(static_cast<size_t>(n_node) * K, 0);
  std::vector<int> tb_mintip(static_cast<size_t>(n_node) * K, INF_TIP);
  for (int t = 0; t < n_tip; ++t) {
    uint32_t set = observed[t];
    // Only a tip observed as exactly one concrete state contributes tie-break
    // support -- an ambiguous tip (e.g. "?", or one the primary codes absent)
    // must not bias which state the uppass prefers.
    if (set != 0 && (set & (set - 1)) == 0) {
      int s = ctz64(set);
      tb_cnt[static_cast<size_t>(t) * K + s] = 1;
      tb_mintip[static_cast<size_t>(t) * K + s] = t;
    }
  }
  // Accumulate over CANONICAL subtrees, so the support counts are a property of
  // the unrooted tree rather than of the incoming rooting (T-374).
  for (int node : co.post) {
    int nk = co.kidNum[node];
    if (nk == 0) continue;
    const int* kid = &co.kids[co.kidOff[node]];
    size_t nb = static_cast<size_t>(node) * K;
    for (int k = 0; k < nk; ++k) {
      size_t cb = static_cast<size_t>(kid[k]) * K;
      for (int s = 0; s < K; ++s) {
        tb_cnt[nb + s] += tb_cnt[cb + s];
        tb_mintip[nb + s] = std::min(tb_mintip[nb + s], tb_mintip[cb + s]);
      }
    }
  }

  // Resolve `state_sets[node]` to a single state, preferring the best-supported
  // state (max subtree count; ties broken by smallest supporting tip index —
  // a strict order, since distinct states never share a supporting tip).  When
  // no state in the set has concrete support (the whole subtree is ambiguous
  // for this character), fall back to the lowest set bit: every node then
  // inherits it and no mismatch is affected, so the choice is score-neutral.
  auto pick_state = [&](int node) -> uint32_t {
    uint32_t set = state_sets[node];
    size_t base = static_cast<size_t>(node) * K;
    int best = -1;
    for (int s = 0; s < K; ++s) {
      if (!(set & (1u << s)) || tb_cnt[base + s] == 0) continue;
      if (best < 0 ||
          tb_cnt[base + s] > tb_cnt[base + best] ||
          (tb_cnt[base + s] == tb_cnt[base + best] &&
           tb_mintip[base + s] < tb_mintip[base + best])) {
        best = s;
      }
    }
    if (best < 0) return set & (~set + 1);   // no support: lowest bit (neutral)
    return 1u << best;
  };

  // --- Uppass: resolve each node to a single state, canonical preorder ---
  // The canonical root is tip 0, whose label is observed, not inferred, so it
  // resolves within its own observation rather than within the downpass set.
  const int root = co.post.back();
  state_sets[root] &= observed[root];
  if (state_sets[root] == 0) state_sets[root] = observed[root];
  state_sets[root] = pick_state(root);

  for (int i = static_cast<int>(co.post.size()) - 1; i >= 0; --i) {
    int node = co.post[i];
    int nk = co.kidNum[node];
    // A canonical leaf has no children to resolve, and forming
    // `&co.kids[co.kidOff[node]]` for one can dereference co.kids.end():
    // kidOff is written as the CURRENT kids.size() when the DFS pops a node,
    // so every childless node popped after the final push_back carries the
    // end offset -- always the last node popped, and usually several more
    // (>1 in 843 of 900 random 2-24 tip trees).  That is the OOB read
    // -D_GLIBCXX_ASSERTIONS aborts on.  The two loops above already skip on
    // nk == 0; this one did not (agent-issues/TreeSearch#51).
    if (nk == 0) continue;
    const int* kid = &co.kids[co.kidOff[node]];
    // Resolve each child: prefer parent's (already-resolved) state if it lies
    // in the child's set (DELTRAN-style); otherwise pick order-invariantly.
    for (int k = 0; k < nk; ++k) {
      int child = kid[k];
      if (state_sets[child] & state_sets[node]) {
        state_sets[child] = state_sets[node]; // inherit parent's state
      } else {
        state_sets[child] = pick_state(child);
      }
    }
  }

  return steps;
}

// Check if two state sets have disjoint states (= mismatch for HSJ).
static inline bool states_mismatch(uint32_t a, uint32_t b) {
  return (a & b) == 0;
}

// Score one hierarchy block via the HSJ a(n)/p(n) dynamic programming.
//
// For each internal node n with children c1, c2:
//   a(n) = min cost assuming primary ABSENT at n
//   p(n) = min cost assuming primary PRESENT at n
//
// Branch cost(parent_state, child_state):
//   absent → absent:   0
//   absent → present:  1  (gain)
//   present → absent:  1  (loss)
//   present → present: alpha * d(parent, child) / m
//     where d = number of secondaries with disjoint state sets
static double score_hierarchy_block(
    const TreeState& tree,
    const HierarchyBlock& block,
    double alpha,
    const std::vector<int>& tip_labels,
    int n_orig_chars,
    int inapp_state,
    const std::vector<uint32_t>& token_states,
    int n_levels)
{
  const int n_tip = tree.n_tip;
  const int n_node = tree.n_node;
  const int m = block.n_secondaries;
  const double INF = std::numeric_limits<double>::infinity();

  // Step 1: Determine primary feasibility at each tip via SET MEMBERSHIP.
  // `label` is a TOKEN (allLevels/contrast-row) index, but block.absent_state
  // and inapp_state are STATE (levels) indices -- comparing them directly via
  // `==` (the former code) was T-375/T-376's bug: it happened to work only
  // when `levels` and `allLevels` coincide, which is common but not
  // guaranteed, and it always mis-scored "?" (whose token index is never
  // itself a valid comparison target for either state index). token_states
  // translates the token into the state-space bitmask it actually denotes, so
  // both sides of the test are in the same index space.
  //
  // The structure is absent when the primary's state set includes the
  // explicit "absent" state (block.absent_state, e.g. "0") OR the
  // inapplicable state ("-"). This mirrors the x-transform recoding
  // (recode_hierarchy.R), which treats `pri == "0" || pri == "-"` as absent,
  // and is required for nested hierarchies where a controlling primary may
  // itself be inapplicable. a(leaf) = 0 iff the token's state set meets the
  // absent-coding states; p(leaf) = 0 iff it meets the present-coding states
  // (anything else). A fully ambiguous token (e.g. "?") meets both, so
  // a = p = 0 there: genuinely unconstrained, not "present" (the old code's
  // effective classification of "?").
  const uint32_t inapp_bit = (inapp_state >= 0) ? (1u << inapp_state) : 0u;
  const uint32_t absent_bits = (1u << block.absent_state) | inapp_bit;

  // Tips where the primary CANNOT code the structure present are exactly the
  // tips at which the secondaries cannot apply, so they must not constrain the
  // secondary reconstruction -- see fitch_label_char() (T-374).
  //
  // The test is deliberately strict ("cannot be present"), not the laxer "may
  // be absent".  Under the lax test a tip whose primary is "?" would be freed
  // too, discarding an OBSERVED secondary -- but observing a secondary is
  // itself evidence the structure is present, and that is real information the
  // alpha term should keep.  It would also be silently self-erasing: `domain`
  // in fitch_label_char() unions only the non-free tips, so a block with no
  // unambiguously present primary would collapse the whole alpha term to zero.
  // ValidateHierarchy whitelists "?" primaries, so that is reachable data.
  // This matches recode_hierarchy.R's `tipStates == -2L` / `tipSecKnown` path,
  // added under T-379 for exactly this reason: an ambiguous primary must not
  // free the secondaries that WERE observed.
  std::vector<char> pri_free(n_tip, 0);
  for (int t = 0; t < n_tip; ++t) {
    uint32_t set = token_states[tip_labels[t * n_orig_chars + block.primary_char]];
    pri_free[t] = ((set & ~absent_bits) == 0) ? 1 : 0;
  }

  // Step 2: Run Fitch downpass/uppass on each secondary character, over a
  // traversal rooted canonically at tip 0 so the labelling -- and hence
  // d(u, v) -- depends only on the unrooted topology (T-374).  Built once and
  // shared by every secondary in the block.
  std::vector<uint32_t> sec_states(m * n_node, 0);
  if (m > 0) {
    const CanonOrder co = build_canon_order(tree);
    std::vector<uint32_t> buf(n_node);
    for (int j = 0; j < m; ++j) {
      fitch_label_char(tree, tip_labels, block.secondary_chars[j],
                       n_orig_chars, token_states, n_levels, pri_free, co, buf);
      for (int nd = 0; nd < n_node; ++nd) {
        sec_states[j * n_node + nd] = buf[nd];
      }
    }
  }

  // Step 3: a(n)/p(n) DP
  // a[node], p[node]
  std::vector<double> a(n_node, 0.0);
  std::vector<double> p(n_node, 0.0);

  // Initialize leaves
  for (int t = 0; t < n_tip; ++t) {
    int label = tip_labels[t * n_orig_chars + block.primary_char];
    uint32_t set = token_states[label];
    a[t] = (set & absent_bits) ? 0.0 : INF;
    p[t] = (set & ~absent_bits) ? 0.0 : INF;
  }

  // Helper: count secondary mismatches between two nodes
  auto count_mismatches = [&](int node1, int node2) -> int {
    int d = 0;
    for (int j = 0; j < m; ++j) {
      if (states_mismatch(sec_states[j * n_node + node1],
                          sec_states[j * n_node + node2])) {
        ++d;
      }
    }
    return d;
  };

  // Postorder traversal
  for (int i = 0; i < static_cast<int>(tree.postorder.size()); ++i) {
    int node = tree.postorder[i];
    int ni = node - n_tip;
    int c1 = tree.left[ni];
    int c2 = tree.right[ni];

    // Compute branch costs for each combination of parent/child states
    // parent absent:
    double bc_aa_c1 = 0.0;  // absent→absent
    double bc_ap_c1 = 1.0;  // absent→present (gain)
    double bc_aa_c2 = 0.0;
    double bc_ap_c2 = 1.0;

    // parent present:
    double bc_pa_c1 = 1.0;  // present→absent (loss)
    double bc_pa_c2 = 1.0;
    double bc_pp_c1 = 0.0;  // present→present (secondary dissimilarity)
    double bc_pp_c2 = 0.0;

    if (m > 0 && alpha > 0.0) {
      int d1 = count_mismatches(node, c1);
      int d2 = count_mismatches(node, c2);
      bc_pp_c1 = alpha * d1 / m;
      bc_pp_c2 = alpha * d2 / m;
    }

    // a(n): parent is absent
    double best_a = INF;
    // Try all 4 combinations of child states
    double cost_aa = (a[c1] + bc_aa_c1) + (a[c2] + bc_aa_c2);  // both absent
    double cost_ap = (a[c1] + bc_aa_c1) + (p[c2] + bc_ap_c2);  // c1 abs, c2 pres
    double cost_pa = (p[c1] + bc_ap_c1) + (a[c2] + bc_aa_c2);  // c1 pres, c2 abs
    double cost_pp = (p[c1] + bc_ap_c1) + (p[c2] + bc_ap_c2);  // both present
    best_a = std::min({cost_aa, cost_ap, cost_pa, cost_pp});
    a[node] = best_a;

    // p(n): parent is present
    double best_p = INF;
    cost_aa = (a[c1] + bc_pa_c1) + (a[c2] + bc_pa_c2);  // both absent
    cost_ap = (a[c1] + bc_pa_c1) + (p[c2] + bc_pp_c2);  // c1 abs, c2 pres
    cost_pa = (p[c1] + bc_pp_c1) + (a[c2] + bc_pa_c2);  // c1 pres, c2 abs
    cost_pp = (p[c1] + bc_pp_c1) + (p[c2] + bc_pp_c2);  // both present
    best_p = std::min({cost_aa, cost_ap, cost_pa, cost_pp});
    p[node] = best_p;
  }

  // Root score = min(a[root], p[root])
  int root = tree.postorder.back();
  return std::min(a[root], p[root]);
}

double hsj_score(
    TreeState& tree,
    const DataSet& ds,
    const std::vector<HierarchyBlock>& hierarchy_blocks,
    double alpha,
    const std::vector<int>& tip_labels,
    int n_orig_chars)
{
  // Score non-hierarchy characters via standard Fitch.
  // Use fitch_score_ew() directly to avoid infinite recursion when
  // score_tree() dispatches to hsj_score() for HSJ mode.
  double fitch_total = fitch_score_ew(tree, ds);

  // Score each hierarchy block via HSJ DP
  double hsj_total = 0.0;
  for (const auto& block : hierarchy_blocks) {
    hsj_total += score_hierarchy_block(
        tree, block, alpha, tip_labels, n_orig_chars, ds.inapp_state,
        ds.token_states, ds.n_levels);
  }

  return fitch_total + hsj_total;
}

double hsj_score(TreeState& tree, const DataSet& ds)
{
  return hsj_score(tree, ds, ds.hierarchy_blocks, ds.hsj_alpha,
                   ds.tip_labels, ds.n_orig_chars);
}

} // namespace ts
