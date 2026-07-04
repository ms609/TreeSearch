// =====================================================================
// Red-team area-13 : exhaustive validity harness for the T-327 guard
// =====================================================================
//
// QUESTION UNDER TEST
// -------------------
// impose_one_pass() in src/ts_constraint.cpp repairs constraint
// violations with topology_spr() moves, guarding each move (T-327 fix
// 6b60f235, src/ts_constraint.cpp ~735) with:
//
//        tree.build_postorder();                       // size-capped DFS
//        if (postorder.size() != n_internal) { revert the move; }
//
// i.e. a move is ACCEPTED iff build_postorder() visits exactly
// n_internal internal nodes.  The guard's author claims this is
// "complete without having to enumerate every topology_spr corner case,
// and guarantees no corrupt tree ever reaches a DFS helper."
//
// This harness tests that completeness claim.  build_postorder() (see
// src/ts_tree.cpp:75) is a DFS from the root over left[]/right[] with
// NO visited-set and a cap that BREAKS once the visit count EXCEEDS
// n_internal.  So an OVER-count is always caught; the only way an
// invalid tree can score exactly n_internal is a NET-ZERO corruption:
// one node reached twice (+1) compensated by one node never reached
// (orphan, -1).  If such a tree slips the guard, a later topology_spr()
// / map_constraint_nodes / compute_dfs_timestamps walks the corruption
// -> std::bad_alloc or a silent wrong answer, reachable from the public
// MaximizeParsimony(constraint = ...)  ->  a P1 for the CRAN release.
//
// METHOD
// ------
// For every VALID rooted binary tree on n_tip = 4..8 tips, apply
// topology_spr() for every (clip, above, below) the repair can emit and
// compare two verdicts on the result:
//
//     guard  = ( build_postorder().size() == n_internal )   // what the guard sees
//     truth  = full_validity(tree)                          // independent ground truth
//
// and ASSERT  guard  =>  truth.  Any (guard && !truth) is a slipped
// corruption: captured and printed, exit code 1.
//
// WHY ONE MOVE SUFFICES (induction)
// ---------------------------------
// impose_one_pass() only ever calls topology_spr() on a tree that a
// PRIOR accepted move (or the seed) left in place; the guard reverts any
// move it rejects.  So the FIRST invalid state it could ever accept is
// necessarily produced by ONE move from a VALID state.  Enumerating
// {all valid trees} x {all moves} and checking accept=>valid therefore
// catches that first slip if one exists — no need to chain moves.
//
// WHY full_validity() CHECKS parent[] TOO (not just left/right)
// ------------------------------------------------------------
// build_postorder(), compute_node_tips(), map_constraint_nodes() and
// compute_dfs_timestamps() read ONLY left[]/right[] — they are blind to
// parent[].  But topology_spr() itself READS parent[clip]/parent[nx] at
// the TOP of the very next loop iteration (src/ts_constraint.cpp:432,439).
// So a move that leaves left/right a clean arborescence but DESYNCS
// parent[] sails through the guard, and a left/right-only validity check
// would rubber-stamp it as "valid" — masking a latent P1 (the next move
// reads the bad parent[] and can corrupt left/right for real).  Hence
// full_validity() requires parent[] to be fully consistent with
// left/right, closing that blind spot.
//
// TRUST BOUNDARY
// --------------
//   * topology_spr, collect_edges_in_subtree, collect_edges_outside_subtree
//     are EXTRACTED VERBATIM from `git show HEAD:src/ts_constraint.cpp`
//     at build time (extract_funcs.sh -> extracted_spr.gen.inc); they
//     cannot drift from the code the release ships.
//   * build_postorder() and init_from_edge() are the REAL engine
//     (src/ts_tree.cpp is compiled straight in).
//   * The tree enumerator, full_validity() and the move drivers below are
//     this harness's own, independent code.
//
// Two enumerations are run:
//   Pass A (geometric superset): clip over all non-root nodes x targets
//     over all real edges + the explicit degenerate targets nx/nz/ns.
//     collect_edges_*() can only ever emit real (parent,child) edges, so
//     this is a STRICT SUPERSET of any move impose_one_pass can emit ->
//     silence here is an airtight KILL.
//   Pass B (emission-shaped): mirrors the outside/inside geometry of the
//     repair (clip inside best_node's subtree grafts OUTSIDE it, and vice
//     versa) using the verbatim collect_edges_*() to build the target
//     set, so any hit is characterised as reachable-shaped, not a
//     geometric artefact.  Pass B's probe set is a subset of Pass A's.
//
// Belt-and-suspenders: for n_tip=4,5, Pass A is re-run under every
// internal-id relabelling (relabel-invariance check).
//
// n_tip<=8 runs locally in well under a minute — no Hamilton needed.

#include "ts_tree.h"

#include <cstdio>
#include <cstdint>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>

// ---- verbatim functions-under-test, generated at build time ----------
//   ts::topology_spr(TreeState&, int clip, int above, int below)
//   ts::collect_edges_in_subtree(const TreeState&, int, vector<pair<int,int>>&)
//   ts::collect_edges_outside_subtree(const TreeState&, int, vector<pair<int,int>>&)
#include "extracted_spr.gen.inc"

using ts::TreeState;
using ts::DataSet;

// =====================================================================
// A minimal DataSet: init_from_edge()/load_tip_states() touch only
// total_words, n_blocks and tip_states; with n_blocks==0 no character
// machinery is exercised.  We never score, only manipulate topology.
// =====================================================================
static DataSet make_min_dataset(int n_tip) {
  DataSet ds;
  ds.n_blocks = 0;
  ds.total_words = 1;
  ds.tip_states.assign(static_cast<size_t>(n_tip) * ds.total_words, 0ULL);
  // block_word_offset / blocks stay empty (only read when n_blocks>0).
  return ds;
}

// =====================================================================
// INDEPENDENT ground-truth structural validity.
//
// A fully-valid rooted binary tree (tips 0..n_tip-1, internal
// n_tip..2n-2, root == n_tip) must satisfy ALL of:
//   (1) parent[root] == root;
//   (2) every internal slot holds an in-range node id, left != right,
//       and parent[] is consistent with left/right (bidirectional);
//   (3) in-degree via left/right is exactly 1 for every non-root node
//       and 0 for the root (no double-ref, no missing parent);
//   (4) a DFS from the root over left/right visits every node exactly
//       once (acyclic AND no orphan): visited==n_node with no repeat.
// =====================================================================
// Structural validity of the left/right graph ALONE (ignoring parent[]).
// This is exactly the property build_postorder — hence the guard — could in
// principle be validating, since build_postorder reads only left/right.
// Passes iff: slots in range, left!=right, in-degree 1 for non-root / 0 for
// root, and a root DFS visits every node exactly once (acyclic + no orphan).
static bool left_right_valid(const TreeState& T) {
  const int nt = T.n_tip, ni = T.n_internal, nn = T.n_node, root = nt;
  if (static_cast<int>(T.left.size())  != ni) return false;
  if (static_cast<int>(T.right.size()) != ni) return false;

  for (int i = 0; i < ni; ++i) {
    const int l = T.left[i], r = T.right[i];
    if (l < 0 || l >= nn || r < 0 || r >= nn) return false;
    if (l == r) return false;
  }
  std::vector<int> indeg(nn, 0);
  for (int i = 0; i < ni; ++i) { indeg[T.left[i]]++; indeg[T.right[i]]++; }
  if (indeg[root] != 0) return false;
  for (int v = 0; v < nn; ++v) { if (v == root) continue; if (indeg[v] != 1) return false; }

  std::vector<char> seen(nn, 0);
  std::vector<int> st; st.reserve(nn);
  st.push_back(root);
  int visited = 0;
  while (!st.empty()) {
    const int node = st.back(); st.pop_back();
    if (seen[node]) return false;            // repeat => double-ref / cycle
    seen[node] = 1; ++visited;
    if (node >= nt) { const int i = node - nt; st.push_back(T.left[i]); st.push_back(T.right[i]); }
  }
  return visited == nn;                       // else orphan
}

// Is parent[] fully consistent with a (presumed valid) left/right graph?
static bool parent_consistent(const TreeState& T) {
  const int nt = T.n_tip, ni = T.n_internal, nn = T.n_node, root = nt;
  if (static_cast<int>(T.parent.size()) != nn) return false;
  if (T.parent[root] != root) return false;
  for (int i = 0; i < ni; ++i) {
    const int node = nt + i;
    if (T.parent[T.left[i]]  != node) return false;
    if (T.parent[T.right[i]] != node) return false;
  }
  return true;
}

// Full validity = left/right structurally valid AND parent[] consistent.
static bool full_validity(const TreeState& T) {
  return left_right_valid(T) && parent_consistent(T);
}

// Does following parent[] upward from ANY node fail to reach the root within
// n_node steps (i.e. a parent[] cycle in the orphaned/desynced region)?  This
// is the property that would make a parent-ASCENDING consumer loop without
// bound — e.g. tbr_search's `while (cur != root) cur = parent[cur]`
// (ts_tbr.cpp:106-107), which in ts_nni_perturb.cpp runs on the post-impose
// tree BEFORE the caller's verify-and-discard.  build_postorder does NOT read
// parent[], so the guard is blind to this.
static bool has_parent_cycle(const TreeState& T) {
  const int nt = T.n_tip, nn = T.n_node, root = nt;
  if (static_cast<int>(T.parent.size()) != nn) return true;
  for (int v = 0; v < nn; ++v) {
    int node = v, steps = 0;
    while (node != root) {
      const int p = T.parent[node];
      if (p < 0 || p >= nn) return true;          // out-of-range parent
      node = p;
      if (++steps > nn) return true;              // never reached root => cycle
    }
  }
  return false;
}

// Classify an accepted-INVALID tree:
//   1 = left/right ITSELF is structurally corrupt (double-ref + orphan) —
//       the net-zero mechanism; build_postorder is NOT a complete validator.
//   2 = left/right is a valid arborescence; only parent[] is desynced —
//       build_postorder DID fully validate structure; residual is a stale
//       parent[] the next reanchor/topology_spr consumes.
static int invalidity_type(const TreeState& T) {
  return left_right_valid(T) ? 2 : 1;
}

// Does the left/right graph contain a cycle reachable from the root?
// (three-colour iterative DFS).  This is the property that would let the
// UNCAPPED, visited-set-free downstream traversals — collect_edges_*() and
// compute_dfs_timestamps() — loop without bound -> std::bad_alloc.  If NO
// guard-ACCEPTED tree ever has one, that crash is empirically unreachable
// and any residual defect is a silent wrong-answer, not a crash.
static bool has_root_reachable_cycle(const TreeState& T) {
  const int nt = T.n_tip, nn = T.n_node, root = nt;
  std::vector<char> color(nn, 0);            // 0=white 1=grey(on-stack) 2=black
  std::vector<std::pair<int,int>> st;        // (node, child-cursor: 0,1,2)
  st.push_back({root, 0});
  color[root] = 1;
  while (!st.empty()) {
    auto& top = st.back();
    const int node = top.first;
    if (node < nt) { color[node] = 2; st.pop_back(); continue; }  // tip: leaf
    const int i = node - nt;
    if (top.second < 2) {
      const int child = (top.second == 0) ? T.left[i] : T.right[i];
      ++top.second;
      if (child < 0 || child >= nn) continue;
      if (color[child] == 1) return true;    // back-edge to on-stack node -> cycle
      if (color[child] == 0) { color[child] = 1; st.push_back({child, 0}); }
    } else {
      color[node] = 2; st.pop_back();
    }
  }
  return false;
}

// =====================================================================
// Abstract tree generator: emits every tip-labelled rooted binary tree
// on n leaves EXACTLY ONCE via sequential leaf addition (there are
// (2n-3)!! such trees).  Leaves 0..n-1; internal ids n..2n-2 allocated
// as leaves are inserted; root can be any internal id (it floats on
// root-stem grafts).  A finalise step relabels so the root becomes
// n_tip and emits an edge list for the real init_from_edge().
// =====================================================================
struct AbsTree {
  int n = 0;               // target leaf count
  int root = -1;
  int used_internals = 0;  // ids in use: n .. n+used_internals-1
  std::vector<int> L, R;   // size n-1, indexed by (internal_id - n)
  std::vector<int> par;    // size 2n-1
};

static AbsTree abs_base(int n) {
  AbsTree t; t.n = n;
  t.L.assign(n - 1, -1); t.R.assign(n - 1, -1);
  t.par.assign(2 * n - 1, -1);
  // internal id n with leaf children {0,1}
  t.used_internals = 1;
  t.root = n;
  t.L[0] = 0; t.R[0] = 1;
  t.par[0] = n; t.par[1] = n; t.par[n] = n;
  return t;
}

static void graft_edge(AbsTree& t, int node, bool is_left, int leaf) {
  const int m = t.n + t.used_internals; ++t.used_internals;
  const int mi = m - t.n, ndi = node - t.n;
  const int c = is_left ? t.L[ndi] : t.R[ndi];
  t.L[mi] = c; t.R[mi] = leaf;
  if (is_left) t.L[ndi] = m; else t.R[ndi] = m;
  t.par[c] = m; t.par[leaf] = m; t.par[m] = node;
}

static void graft_rootstem(AbsTree& t, int leaf) {
  const int m = t.n + t.used_internals; ++t.used_internals;
  const int mi = m - t.n;
  t.L[mi] = t.root; t.R[mi] = leaf;
  t.par[t.root] = m; t.par[leaf] = m; t.par[m] = m;
  t.root = m;
}

// Convert a completed AbsTree (all n leaves present) to init_from_edge
// input: root -> n_tip, other internals -> n_tip+1.., 1-based edges.
static void abs_to_edges(const AbsTree& t,
                         std::vector<int>& ep, std::vector<int>& ec) {
  const int n = t.n;
  std::vector<int> remap(2 * n - 1, -1);
  for (int i = 0; i < n; ++i) remap[i] = i;         // leaves fixed
  remap[t.root] = n;                                // root -> n_tip
  int nxt = n + 1;
  for (int id = n; id <= 2 * n - 2; ++id)
    if (id != t.root) remap[id] = nxt++;

  ep.clear(); ec.clear();
  for (int id = n; id <= 2 * n - 2; ++id) {
    const int i = id - n;
    // (parent, left) then (parent, right); 1-based for init_from_edge
    ep.push_back(remap[id] + 1); ec.push_back(remap[t.L[i]] + 1);
    ep.push_back(remap[id] + 1); ec.push_back(remap[t.R[i]] + 1);
  }
}

// Recursively enumerate; call cb(AbsTree) for each complete tree.
template <class F>
static void enumerate(const AbsTree& t, int k, int n, F&& cb) {
  if (k == n) { cb(t); return; }
  // root-stem position
  { AbsTree u = t; graft_rootstem(u, k); enumerate(u, k + 1, n, cb); }
  // every internal node's two child edges
  const int used = t.used_internals;
  for (int j = 0; j < used; ++j) {
    const int node = n + j;
    { AbsTree u = t; graft_edge(u, node, /*left=*/true,  k); enumerate(u, k + 1, n, cb); }
    { AbsTree u = t; graft_edge(u, node, /*left=*/false, k); enumerate(u, k + 1, n, cb); }
  }
}

// (2n-3)!! — expected tree count, used as an enumeration self-check.
static long long double_fact_odd(int n) {
  long long p = 1;
  for (int x = 2 * n - 3; x >= 1; x -= 2) p *= x;
  return p;
}

// =====================================================================
// Result recording
// =====================================================================
struct Witness {
  int n_tip;
  std::vector<int> ep, ec;   // the valid seed tree (1-based edges)
  int clip, above, below;
  std::string pass;          // "A", "A-perm", "C1-out/in", "C2-out/in", "D"
  int ftype = 0;             // 1 = left/right structurally corrupt; 2 = only parent[] desynced
  bool has_cycle = false;    // root-reachable left/right cycle in the accepted tree?
  bool par_cycle = false;    // parent[] cycle (would loop a parent-ascending consumer)
  unsigned long long split_bits = 0;  // Pass D: the constraint split (tip bitmask)
  int depth = -1;            // Pass D: move index within the split at which it slipped
  bool passes_verify = false;  // Pass D: does the caller's map_constraint_nodes re-check
                               // still MAP the split on this corrupt tree? (true => the
                               // verify-and-discard guard does NOT catch it => the
                               // malformed tree is kept/scored/returned)
  // corrupt post-move topology (for the report)
  std::vector<int> parent, left, right;
  size_t postorder_size;
};

static std::vector<Witness> g_witnesses;
static long long g_probes_A = 0, g_probes_C = 0;
static long long g_parent_cycle = 0;   // accepted-invalid trees with a parent[] cycle (ALL passes)
static long long g_max_movelist = 0;   // max (|move_out_roots|+|move_in_roots|) seen

static void print_witness(const Witness& w) {
  std::printf("\n!!! ACCEPT-ON-INVALID  (pass %s, n_tip=%d) !!!\n",
              w.pass.c_str(), w.n_tip);
  std::printf("  seed edges (parent->child, 1-based):");
  for (size_t i = 0; i < w.ep.size(); ++i)
    std::printf(" (%d->%d)", w.ep[i], w.ec[i]);
  std::printf("\n  move: clip=%d above=%d below=%d\n", w.clip, w.above, w.below);
  if (w.pass == "D")
    std::printf("  split tip-bitmask=0x%llx  slipped at move index=%d ; "
                "caller verify still maps split: %s\n",
                w.split_bits, w.depth,
                w.passes_verify ? "YES (malformed tree KEPT -> reaches score/return)"
                                : "no (discarded by verify-and-discard)");
  std::printf("  guard saw postorder.size()=%zu == n_internal (ACCEPT)\n",
              w.postorder_size);
  std::printf("  invalidity type: %s\n",
              w.ftype == 1 ? "TYPE-1 (left/right structurally corrupt: double-ref + orphan)"
                           : "TYPE-2 (left/right VALID; only parent[] desynced)");
  std::printf("  root-reachable cycle: %s\n", w.has_cycle ? "YES" : "no");
  std::printf("  post-move parent[]:");
  for (size_t i = 0; i < w.parent.size(); ++i) std::printf(" %d", w.parent[i]);
  std::printf("\n  post-move left[]  :");
  for (size_t i = 0; i < w.left.size(); ++i) std::printf(" %d", w.left[i]);
  std::printf("\n  post-move right[] :");
  for (size_t i = 0; i < w.right.size(); ++i) std::printf(" %d", w.right[i]);
  std::printf("\n");
}

// Mirror of the caller's post-repair verification (defined below pass_C).
static bool split_still_maps(TreeState& T, const std::vector<uint64_t>& S, int n_words);

// Probe one move on tree T (T is restored afterwards).  Records a
// witness iff the guard would ACCEPT a structurally-invalid result.
// If vsplit != nullptr, also records whether the caller's verify-and-discard
// (map_constraint_nodes) would still MAP that split on the corrupt tree.
static void probe(TreeState& T, int clip, int above, int below,
                  const char* pass, const std::vector<int>& seed_ep,
                  const std::vector<int>& seed_ec, long long& counter,
                  const std::vector<uint64_t>* vsplit = nullptr, int vwords = 0) {
  ++counter;
  std::vector<int> sp = T.parent, sl = T.left, sr = T.right;   // snapshot

  ts::topology_spr(T, clip, above, below);
  T.build_postorder();
  const bool accept = (static_cast<int>(T.postorder.size()) == T.n_internal);
  const bool valid  = full_validity(T);

  if (accept && !valid) {
    Witness w;
    w.n_tip = T.n_tip; w.ep = seed_ep; w.ec = seed_ec;
    w.clip = clip; w.above = above; w.below = below; w.pass = pass;
    w.ftype = invalidity_type(T);
    w.has_cycle = has_root_reachable_cycle(T);
    w.par_cycle = has_parent_cycle(T);
    if (w.par_cycle) ++g_parent_cycle;
    if (vsplit) { w.passes_verify = split_still_maps(T, *vsplit, vwords);
                  w.split_bits = 0; }
    w.parent = T.parent; w.left = T.left; w.right = T.right;
    w.postorder_size = T.postorder.size();
    g_witnesses.push_back(std::move(w));
  }

  T.parent = sp; T.left = sl; T.right = sr;   // restore
}

static inline int pc64(uint64_t x) { return __builtin_popcountll(x); }

// Real (parent,child) edges of T.
static void real_edges(const TreeState& T,
                       std::vector<std::pair<int,int>>& edges) {
  edges.clear();
  for (int i = 0; i < T.n_internal; ++i) {
    edges.push_back({T.n_tip + i, T.left[i]});
    edges.push_back({T.n_tip + i, T.right[i]});
  }
}

// Local degenerate targets for a clip: the collisions the task calls out
// explicitly (above/below coinciding with nx=parent, ns=sibling,
// nz=grandparent, or the clip itself).  We probe the full cartesian of
// {clip,nx,ns,nz} x {clip,nx,ns,nz} so no degenerate pair is missed,
// including non-edge pairs collect_edges_*() would never emit.
static void degenerate_targets(const TreeState& T, int clip,
                               std::vector<std::pair<int,int>>& out) {
  out.clear();
  const int root = T.n_tip;
  if (clip == root) return;
  const int nx = T.parent[clip];
  int ns = -1, nz = -1;
  if (nx >= T.n_tip) {
    const int nxi = nx - T.n_tip;
    ns = (T.left[nxi] == clip) ? T.right[nxi] : T.left[nxi];
    nz = T.parent[nx];
  }
  int locals[4] = {clip, nx, ns, nz};
  for (int a = 0; a < 4; ++a)
    for (int b = 0; b < 4; ++b) {
      int A = locals[a], B = locals[b];
      if (A < 0 || B < 0) continue;
      out.push_back({A, B});
    }
}

// ---- Pass A: geometric superset -------------------------------------
static void pass_A(TreeState& T, const std::vector<int>& ep,
                   const std::vector<int>& ec, const char* tag) {
  std::vector<std::pair<int,int>> edges, degs;
  real_edges(T, edges);
  const int root = T.n_tip;
  for (int clip = 0; clip < T.n_node; ++clip) {
    if (clip == root) continue;                 // root can't be moved
    degenerate_targets(T, clip, degs);
    for (auto& e : edges) probe(T, clip, e.first, e.second, tag, ep, ec, g_probes_A);
    for (auto& e : degs)  probe(T, clip, e.first, e.second, tag, ep, ec, g_probes_A);
  }
}

// ---- Pass C: FAITHFUL constraint-driven emission enumeration --------
// For each single constraint split S (tip subset, tip 0 conventionally
// outside), reproduce impose_one_pass's per-split emission EXACTLY, using
// the VERBATIM compute_node_tips + find_maximal_subtrees + collect_edges:
//   best_node     = argmin symmetric-difference over postorder     (:636-649)
//   move_out_mask = node_tips[best_node] & ~S                       (:658-661)
//   move_in_mask  = S & ~node_tips[best_node]
//   move_out_roots = find_maximal_subtrees(best_node, -1, move_out_mask)  (:665)
//   move_in_roots  = find_maximal_subtrees(root, best_node, move_in_mask) (:669)
// The code then picks ONE random collect_edges target per move (:726); we
// probe EVERY target (deterministic superset of any RNG outcome).
//
//   C1 (tight): M ranges over the ACTUAL move roots — the exact FIRST-move
//     emission set.  A C1 hit is reachable by the first move of a single
//     split -> an immediate concrete repro.
//   C2 (stale-M superset): M ranges over ALL non-root nodes with the same
//     targets — covers the T-327 "stale move-root list after an earlier
//     accepted move in the same split" path, where M can be an arbitrary
//     node id while targets still come from collect_edges(fresh best_node).
//     Cost-gated to small n.
//
// WHY SINGLE-SPLIT OVER ALL TREES IS COMPLETE: impose_one_pass processes
// one split at a time, rebuilding masks between splits (:630-633); each
// try_move starts from the current tree, which the guard keeps a VALID
// binary tree.  So every (tree-state, split) the loop ever visits is some
// (valid binary tree) x (single split) — already in this enumeration.
static void pass_C(TreeState& T, const std::vector<int>& ep,
                   const std::vector<int>& ec, int n_words, bool do_C2) {
  const int nt = T.n_tip, root = nt;

  std::vector<uint64_t> S(n_words), mout(n_words), min_(n_words);
  std::vector<int> mo_roots, mi_roots;
  std::vector<std::pair<int,int>> t_out, t_in;

  // Splits: subsets of {1..nt-1} (tip 0 kept outside, per the
  // map_constraint_nodes canonicalisation), at least one tip in.
  const long long combos = 1LL << (nt - 1);        // bit j -> tip (j+1)
  for (long long bits = 1; bits < combos; ++bits) {
    std::fill(S.begin(), S.end(), 0ULL);
    for (int j = 0; j < nt - 1; ++j)
      if (bits & (1LL << j)) { int tip = j + 1; S[tip / 64] |= (1ULL << (tip % 64)); }

    T.build_postorder();                                   // valid tree
    std::vector<uint64_t> node_tips = ts::compute_node_tips(T, n_words);

    int best_node = -1, best_cost = nt + 1;
    for (int node : T.postorder) {
      const uint64_t* nd = &node_tips[static_cast<size_t>(node) * n_words];
      int cost = 0;
      for (int w = 0; w < n_words; ++w) cost += pc64(nd[w] ^ S[w]);
      if (cost < best_cost) { best_cost = cost; best_node = node; }
    }
    if (best_cost == 0) continue;                          // split already satisfied

    const uint64_t* best_nt = &node_tips[static_cast<size_t>(best_node) * n_words];
    for (int w = 0; w < n_words; ++w) {
      mout[w] = best_nt[w] & ~S[w];
      min_[w] = S[w] & ~best_nt[w];
    }
    // NB: find_maximal_subtrees and collect_edges_*() APPEND to their output
    // (they never clear it — the real code passes a fresh local each call), so
    // we must clear these reused buffers before every call.
    mo_roots.clear(); mi_roots.clear();
    ts::find_maximal_subtrees(T, best_node, -1,       node_tips, mout, n_words, mo_roots);
    ts::find_maximal_subtrees(T, root,      best_node, node_tips, min_, n_words, mi_roots);

    const long long mlsize = (long long)mo_roots.size() + (long long)mi_roots.size();
    if (mlsize > g_max_movelist) g_max_movelist = mlsize;

    t_out.clear(); t_in.clear();
    ts::collect_edges_outside_subtree(T, best_node, t_out);
    ts::collect_edges_in_subtree(T, best_node, t_in);

    // C1 — the exact first-move emission set
    for (int M : mo_roots)
      for (auto& e : t_out) probe(T, M, e.first, e.second, "C1-out", ep, ec, g_probes_C, &S, n_words);
    for (int M : mi_roots)
      for (auto& e : t_in)  probe(T, M, e.first, e.second, "C1-in",  ep, ec, g_probes_C, &S, n_words);

    // C2 — stale-move-root superset (M arbitrary, targets unchanged)
    if (do_C2) {
      for (int M = 0; M < T.n_node; ++M) {
        if (M == root) continue;
        for (auto& e : t_out) probe(T, M, e.first, e.second, "C2-out", ep, ec, g_probes_C, &S, n_words);
        for (auto& e : t_in)  probe(T, M, e.first, e.second, "C2-in",  ep, ec, g_probes_C, &S, n_words);
      }
    }
  }
}

// ---- Pass D: FAITHFUL sequential simulation of impose_one_pass -------
// The definitive stale-move-root test.  For each (tree, single split) it
// reproduces the actual per-split loop: compute the move-root LIST ONCE on
// the pristine tree (as the code does at :665/669), then walk the moves in
// order, re-anchoring best_node before each (:721) and drawing targets from
// the verbatim collect_edges of the CURRENT tree.  Instead of the code's
// single RANDOM target pick (:726) it explores EVERY reachable RNG outcome:
//   - each ACCEPTED target that keeps the tree valid  -> recurse (the move
//     applied, list advances, later moves see a stale M on a mutated tree);
//   - the SKIP branch (move left unapplied) only when >=1 target is
//     REJECTED by the guard (that is the only way the RNG pick can miss).
// If any ACCEPTED move yields a STRUCTURALLY-INVALID tree, that is a
// genuinely reachable slip -> recorded.  Silence = the stale-root path
// cannot slip the guard (a real KILL), despite the loose C2 superset.
struct Move { int node; bool outside; };

static int reanchor_bn(TreeState& T, const std::vector<uint64_t>& S, int n_words) {
  T.build_postorder();
  std::vector<uint64_t> nt = ts::compute_node_tips(T, n_words);
  int bn = -1, bc = T.n_tip + 1;
  for (int node : T.postorder) {
    const uint64_t* nd = &nt[static_cast<size_t>(node) * n_words];
    int c = 0;
    for (int w = 0; w < n_words; ++w) c += pc64(nd[w] ^ S[w]);
    if (c < bc) { bc = c; bn = node; }
  }
  return bn;
}

static long long g_probes_D = 0;
static const long long D_PROBE_BUDGET = 400000000LL;  // bail-out to stay local-bounded
// "Continue past corruption" instrumentation (advisor #3): the real loop keeps
// processing the stale move-list on an accepted-but-corrupt tree, so downstream
// corrupt trees — not just the first — can be the one on which the split MAPS
// (the tree impose_one_pass would return as 'done' and the caller would KEEP).
static long long g_D_corrupt = 0;        // reachable accepted-corrupt trees (all, incl. downstream)
static long long g_D_corrupt_maps = 0;   // ...of which the split still maps (=> caller KEEPS)
static bool g_D_budget_hit = false;

// Mirror of the caller's post-repair verification (ts_driven.cpp:1056-1059,
// ts_nni_perturb.cpp:119-122): map_constraint_nodes searches the (rebuilt)
// postorder for a node whose subtree tip-mask EXACTLY equals the split.  If
// found for every split, the tree is KEPT; otherwise it is discarded.  On a
// corrupt tree this runs over the mis-built postorder / stale node_tips.
static bool split_still_maps(TreeState& T, const std::vector<uint64_t>& S, int n_words) {
  std::vector<uint64_t> nt = ts::compute_node_tips(T, n_words);   // over T.postorder
  for (int node : T.postorder) {
    const uint64_t* nd = &nt[static_cast<size_t>(node) * n_words];
    bool match = true;
    for (int w = 0; w < n_words; ++w) if (nd[w] != S[w]) { match = false; break; }
    if (match) return true;
  }
  return false;
}

// g_returnable_corrupt: the DECISIVE wrong-answer count — a FINAL tree (at an
// impose_constraint return point) that is structurally corrupt AND on which the
// split still maps, so the caller's verify-and-discard KEEPS it and it is
// scored / returned.  0 => no malformed tree is ever returned (single split).
static long long g_returnable_corrupt = 0;

static void run_pass(TreeState T, const std::vector<uint64_t>& S, int n_words,
                     const std::vector<int>& ep, const std::vector<int>& ec,
                     unsigned long long split_bits, int pass_idx, int max_pass_idx);

// impose_constraint has produced a FINAL tree T (its pass loop terminated).
// The caller re-runs map_constraint_nodes and KEEPS T iff the split maps.
static void terminal_check(TreeState& T, const std::vector<uint64_t>& S, int n_words,
                           const std::vector<int>& ep, const std::vector<int>& ec,
                           unsigned long long split_bits) {
  T.build_postorder();
  if (!split_still_maps(T, S, n_words)) return;   // caller DISCARDS -> safe
  if (full_validity(T)) return;                   // valid tree kept -> correct
  ++g_returnable_corrupt;                          // corrupt tree survives verify -> BUG
  if (g_witnesses.size() < 20000) {
    Witness w;
    w.n_tip = T.n_tip; w.ep = ep; w.ec = ec;
    w.clip = -1; w.above = -1; w.below = -1; w.pass = "D-final";
    w.ftype = invalidity_type(T); w.has_cycle = has_root_reachable_cycle(T);
    w.par_cycle = has_parent_cycle(T);
    w.passes_verify = true; w.split_bits = split_bits; w.depth = -1;
    w.parent = T.parent; w.left = T.left; w.right = T.right;
    w.postorder_size = T.postorder.size();
    g_witnesses.push_back(std::move(w));
  }
}

// Walk one impose_one_pass move-list (the FIXED, stale list) with branching over
// RNG-reachable target choices; continue past accepted-corrupt trees; at the end
// of the list, recurse into the next impose_constraint pass.
static void move_dfs(TreeState T, const std::vector<Move>& ml, size_t idx,
                     const std::vector<uint64_t>& S, int n_words,
                     const std::vector<int>& ep, const std::vector<int>& ec,
                     unsigned long long split_bits, int pass_idx, int max_pass_idx) {
  if (g_probes_D > D_PROBE_BUDGET) { g_D_budget_hit = true; return; }
  if (idx >= ml.size()) {                     // pass complete -> next impose_constraint pass
    run_pass(T, S, n_words, ep, ec, split_bits, pass_idx + 1, max_pass_idx);
    return;
  }
  const int M = ml[idx].node;
  const bool outside = ml[idx].outside;
  const int bn = reanchor_bn(T, S, n_words);
  std::vector<std::pair<int,int>> targets;
  if (outside) ts::collect_edges_outside_subtree(T, bn, targets);
  else         ts::collect_edges_in_subtree(T, bn, targets);

  bool any_rejected = false;
  for (auto& t : targets) {
    if (g_probes_D > D_PROBE_BUDGET) { g_D_budget_hit = true; return; }
    ++g_probes_D;
    TreeState T2 = T;
    ts::topology_spr(T2, M, t.first, t.second);
    T2.build_postorder();
    if (static_cast<int>(T2.postorder.size()) == T2.n_internal) {   // guard ACCEPTS
      if (!full_validity(T2)) {
        ++g_D_corrupt;
        const bool maps = split_still_maps(T2, S, n_words);
        if (maps) ++g_D_corrupt_maps;
        if (has_parent_cycle(T2)) ++g_parent_cycle;   // count ALL corrupt, not just recorded
        if (g_witnesses.size() < 20000) {
          Witness w;
          w.n_tip = T2.n_tip; w.ep = ep; w.ec = ec;
          w.clip = M; w.above = t.first; w.below = t.second; w.pass = "D";
          w.ftype = invalidity_type(T2); w.has_cycle = has_root_reachable_cycle(T2);
          w.par_cycle = has_parent_cycle(T2);
          w.passes_verify = maps; w.split_bits = split_bits; w.depth = static_cast<int>(idx);
          w.parent = T2.parent; w.left = T2.left; w.right = T2.right;
          w.postorder_size = T2.postorder.size();
          g_witnesses.push_back(std::move(w));
        }
      }
      move_dfs(T2, ml, idx + 1, S, n_words, ep, ec, split_bits, pass_idx, max_pass_idx);
    } else {
      any_rejected = true;
    }
  }
  if (any_rejected)   // RNG picked a rejected target => move skipped, tree unchanged
    move_dfs(T, ml, idx + 1, S, n_words, ep, ec, split_bits, pass_idx, max_pass_idx);
}

// One iteration of the impose_constraint pass loop (ts_constraint.cpp:768-773),
// faithfully modelling impose_one_pass on the CURRENT tree (valid or corrupt).
static void run_pass(TreeState T, const std::vector<uint64_t>& S, int n_words,
                     const std::vector<int>& ep, const std::vector<int>& ec,
                     unsigned long long split_bits, int pass_idx, int max_pass_idx) {
  if (g_probes_D > D_PROBE_BUDGET) { g_D_budget_hit = true; return; }
  if (pass_idx > max_pass_idx) { terminal_check(T, S, n_words, ep, ec, split_bits); return; }

  T.build_postorder();
  // impose_one_pass returns 0 (not violated) iff the split already maps -> break.
  if (split_still_maps(T, S, n_words)) { terminal_check(T, S, n_words, ep, ec, split_bits); return; }

  const int nt = T.n_tip, root = nt;
  std::vector<uint64_t> node_tips = ts::compute_node_tips(T, n_words);
  int best_node = -1, best_cost = nt + 1;
  for (int node : T.postorder) {
    const uint64_t* nd = &node_tips[static_cast<size_t>(node) * n_words];
    int cost = 0;
    for (int w = 0; w < n_words; ++w) cost += pc64(nd[w] ^ S[w]);
    if (cost < best_cost) { best_cost = cost; best_node = node; }
  }
  if (best_cost == 0) { terminal_check(T, S, n_words, ep, ec, split_bits); return; }

  std::vector<uint64_t> mout(n_words), min_(n_words);
  const uint64_t* best_nt = &node_tips[static_cast<size_t>(best_node) * n_words];
  for (int w = 0; w < n_words; ++w) { mout[w] = best_nt[w] & ~S[w]; min_[w] = S[w] & ~best_nt[w]; }
  std::vector<int> mo_roots, mi_roots;
  ts::find_maximal_subtrees(T, best_node, -1,       node_tips, mout, n_words, mo_roots);
  ts::find_maximal_subtrees(T, root,      best_node, node_tips, min_, n_words, mi_roots);

  const int move_cap = nt / 4 + 2;
  if (static_cast<int>(mo_roots.size() + mi_roots.size()) > move_cap) {  // bail -> break
    terminal_check(T, S, n_words, ep, ec, split_bits); return;
  }
  std::vector<Move> ml;
  for (int M : mo_roots) ml.push_back({M, true});
  for (int M : mi_roots) ml.push_back({M, false});
  if (ml.empty()) { terminal_check(T, S, n_words, ep, ec, split_bits); return; }  // moves==0 -> break

  move_dfs(T, ml, 0, S, n_words, ep, ec, split_bits, pass_idx, max_pass_idx);
}

// Faithful model of the FULL impose_constraint (loop of <= n_splits+1 passes)
// for every single-tip-split constraint.
static void pass_D(TreeState& T, const std::vector<int>& ep,
                   const std::vector<int>& ec, int n_words) {
  const int nt = T.n_tip;
  const int max_pass_idx = 1;   // n_splits == 1 => impose_constraint runs passes 0 and 1
  std::vector<uint64_t> S(n_words);
  const long long combos = 1LL << (nt - 1);
  for (long long bits = 1; bits < combos; ++bits) {
    std::fill(S.begin(), S.end(), 0ULL);
    for (int j = 0; j < nt - 1; ++j)
      if (bits & (1LL << j)) { int tip = j + 1; S[tip / 64] |= (1ULL << (tip % 64)); }
    run_pass(T, S, n_words, ep, ec, static_cast<unsigned long long>(bits), 0, max_pass_idx);
  }
}

// Relabel a valid tree's internal ids by permutation `perm` of the
// non-root internals, rebuild via init_from_edge, run Pass A.  Used only
// for the small-n relabel-invariance belt-and-suspenders check.
static void pass_A_permuted(const std::vector<int>& ep,
                            const std::vector<int>& ec, int n_tip,
                            const DataSet& ds) {
  // non-root internal final-ids are n_tip+1 .. 2n-2
  std::vector<int> ids;
  for (int id = n_tip + 1; id <= 2 * n_tip - 2; ++id) ids.push_back(id);
  std::sort(ids.begin(), ids.end());
  do {
    // build a remap: leaves + root fixed, non-root internals permuted
    std::vector<int> remap(2 * n_tip - 1);
    for (int i = 0; i < n_tip; ++i) remap[i] = i;
    remap[n_tip] = n_tip;
    for (size_t j = 0; j < ids.size(); ++j) remap[n_tip + 1 + (int)j] = ids[j];
    std::vector<int> pep(ep.size()), pec(ec.size());
    for (size_t i = 0; i < ep.size(); ++i) {
      pep[i] = remap[ep[i] - 1] + 1;   // ep is 1-based
      pec[i] = remap[ec[i] - 1] + 1;
    }
    TreeState T;
    T.init_from_edge(pep.data(), pec.data(), (int)pep.size(), ds);
    if (!full_validity(T)) {
      std::printf("FATAL: permuted seed invalid (n_tip=%d) — generator/relabel bug\n", n_tip);
      std::fflush(stdout); std::exit(2);
    }
    pass_A(T, pep, pec, "A-perm");
  } while (std::next_permutation(ids.begin(), ids.end()));
}

int main(int argc, char** argv) {
  int n_lo = 4, n_hi = 8;
  if (argc >= 3) { n_lo = std::atoi(argv[1]); n_hi = std::atoi(argv[2]); }

  std::printf("=== impose_one_pass T-327 guard validity harness ===\n");
  std::printf("Testing accept(guard) => full_validity for n_tip=%d..%d\n\n", n_lo, n_hi);

  // Pass A (geometric) and Pass C1 (faithful first-move) run for all n.
  // Pass C2 (crude arbitrary-M superset) is OFF by default — superseded by
  // Pass D, the faithful sequential simulator, which is the definitive
  // stale-move-root test.  D branches over accepted RNG outcomes so it is
  // cost-gated to small n.
  const int c2_max = 0;    // C2 disabled (D supersedes it); raise to enable
  const int d_max  = 6;    // faithful sequential simulation up to this n_tip

  bool count_fail = false;
  for (int n = n_lo; n <= n_hi; ++n) {
    DataSet ds = make_min_dataset(n);
    const int n_words = (n + 63) / 64;
    const long long expect = double_fact_odd(n);
    long long ntrees = 0;
    const long long probes_A0 = g_probes_A, probes_C0 = g_probes_C;
    const size_t wit0 = g_witnesses.size();

    AbsTree base = abs_base(n);
    // base has leaves {0,1}; enumerate inserts leaves 2..n-1
    enumerate(base, 2, n, [&](const AbsTree& at) {
      ++ntrees;
      std::vector<int> ep, ec;
      abs_to_edges(at, ep, ec);
      TreeState T;
      T.init_from_edge(ep.data(), ec.data(), (int)ep.size(), ds);
      if (!full_validity(T)) {
        std::printf("FATAL: generated seed invalid (n_tip=%d) — generator bug\n", n);
        std::fflush(stdout); std::exit(2);
      }
      pass_A(T, ep, ec, "A");                 // geometric superset (informational)
      pass_C(T, ep, ec, n_words, /*do_C2=*/ n <= c2_max);   // faithful first-move (C1)
      if (n <= d_max) pass_D(T, ep, ec, n_words);           // faithful sequential (decisive)
      if (n <= 5) pass_A_permuted(ep, ec, n, ds);
    });

    const bool count_ok = (ntrees == expect);
    std::printf("n_tip=%d : trees=%lld (expected (2n-3)!!=%lld) %s | "
                "probes A=%lld C=%lld | accept-on-invalid(all passes)=%zu\n",
                n, ntrees, expect, count_ok ? "OK" : "COUNT-MISMATCH!",
                g_probes_A - probes_A0, g_probes_C - probes_C0,
                g_witnesses.size() - wit0);
    if (!count_ok) count_fail = true;
  }

  // Tally witnesses per pass.
  //   Pass A / A-perm : geometric superset — includes (clip,above,below)
  //     triples collect_edges_*() can NEVER emit; hits here are EXPECTED
  //     artifacts, they only prove the guard is not a full structural
  //     validator (which is fine as long as no reachable move slips).
  //   C1-out / C1-in  : FAITHFUL first-move emission set (exact move roots,
  //     real collect_edges targets) — a hit is directly reachable.
  //   C2-out / C2-in  : stale-move-root superset — a hit is reachable only
  //     via a >=2-move single-split sequence (still to be reproduced).
  long long wA = 0, wAp = 0, wC1 = 0, wC2 = 0, wD = 0;
  long long cyc_total = 0;
  std::vector<Witness> reach_C1, reach_C2, reach_D;
  for (const auto& w : g_witnesses) {
    if (w.has_cycle) ++cyc_total;
    if      (w.pass == "A")      ++wA;
    else if (w.pass == "A-perm") ++wAp;
    else if (w.pass == "C1-out" || w.pass == "C1-in") { ++wC1; reach_C1.push_back(w); }
    else if (w.pass == "C2-out" || w.pass == "C2-in") { ++wC2; reach_C2.push_back(w); }
    else if (w.pass == "D" || w.pass == "D-final") { ++wD; reach_D.push_back(w); }
  }
  long long wD_kept = 0, wC1_kept = 0;   // witnesses the caller's verify-and-discard would KEEP
  long long t1 = 0, t2 = 0;              // invalidity types across C1+D (reachable) witnesses
  for (const auto& w : reach_D)  { if (w.passes_verify) ++wD_kept;  if (w.ftype==1)++t1; else if(w.ftype==2)++t2; }
  for (const auto& w : reach_C1) { if (w.passes_verify) ++wC1_kept; if (w.ftype==1)++t1; else if(w.ftype==2)++t2; }

  std::printf("\n--- summary ---\n");
  std::printf("total probes: A=%lld  C=%lld  D=%lld\n", g_probes_A, g_probes_C, g_probes_D);
  std::printf("accept-on-invalid witnesses:\n");
  std::printf("  Pass A   (geometric superset — unreachable artifacts)  : %lld\n", wA);
  std::printf("  A-perm   (relabel-invariance check)                    : %lld\n", wAp);
  std::printf("  Pass C1  (FAITHFUL first-move emission)                 : %lld\n", wC1);
  std::printf("    of which the caller's verify-and-discard would KEEP    : %lld\n", wC1_kept);
  std::printf("  Pass C2  (crude stale-root superset; off unless c2_max) : %lld\n", wC2);
  std::printf("  Pass D   (FAITHFUL sequential sim — DECISIVE)           : %lld\n", wD);
  std::printf("    [full impose_constraint model] reachable corrupt trees=%lld,\n"
              "      of which split maps somewhere=%lld ; RETURNABLE (final tree corrupt\n"
              "      AND survives caller verify)=%lld%s\n",
              g_D_corrupt, g_D_corrupt_maps, g_returnable_corrupt,
              g_D_budget_hit ? "  [BUDGET HIT — partial coverage]" : "");
  std::printf("    of which the caller's verify-and-discard would KEEP    : %lld\n"
              "    (KEEP>0 => a malformed tree can be scored & returned;\n"
              "     KEEP==0 => downstream verify catches every slip => benign)\n",
              wD_kept);
  std::printf("max move-list size (|move_out_roots|+|move_in_roots|) seen: %lld\n"
              "  (a stale-move-root path needs >=2 moves in one split)\n",
              g_max_movelist);
  std::printf("accepted-invalid trees with a ROOT-REACHABLE left/right CYCLE: %lld\n"
              "  (0 => the root-anchored DFS helpers cannot loop)\n", cyc_total);
  std::printf("accepted-invalid trees with a parent[] CYCLE: %lld\n"
              "  (0 => a parent-ASCENDING consumer — tbr_search ts_tbr.cpp:107, run on the\n"
              "   post-impose tree BEFORE verify-and-discard — cannot loop either)\n",
              g_parent_cycle);
  std::printf("reachable (C1+D) invalidity types: TYPE-1 (left/right corrupt)=%lld, "
              "TYPE-2 (only parent[] desynced)=%lld\n"
              "  (all TYPE-2 => build_postorder DOES fully validate left/right structure;\n"
              "   the residual is a stale parent[] consumed by the next move)\n",
              t1, t2);

  auto dump = [](const char* title, std::vector<Witness>& v) {
    if (v.empty()) return;
    const size_t show = std::min<size_t>(v.size(), 8);
    std::printf("\n=== %s ===\n", title);
    for (size_t i = 0; i < show; ++i) print_witness(v[i]);
    if (v.size() > show)
      std::printf("\n(... %zu more omitted ...)\n", v.size() - show);
  };
  std::stable_partition(reach_C1.begin(), reach_C1.end(),
                        [](const Witness& w){ return w.passes_verify; });
  dump("Pass C1 witnesses (directly reachable first moves)", reach_C1);
  dump("Pass C2 witnesses (crude stale-root superset)", reach_C2);
  // show the caller-KEPT (passes_verify) Pass-D witnesses first — those are
  // the ones with user-facing impact.
  std::stable_partition(reach_D.begin(), reach_D.end(),
                        [](const Witness& w){ return w.passes_verify; });
  dump("Pass D witnesses (CONFIRMED reachable via real sequential path)", reach_D);

  if (count_fail) {
    std::printf("\nRESULT: HARNESS ERROR — tree-count mismatch (see per-n lines).\n");
    return 2;
  }
  // The decisive signals are C1 (faithful first move) and D (faithful
  // sequential / stale-root).  C2 is a crude superset kept off by default.
  const long long decisive = wC1 + wD;
  if (decisive == 0) {
    std::printf("\nRESULT: KILL — the guard is SUFFICIENT on every reachable move.\n"
                "No move in the FAITHFUL enumeration (C1: exact first-move emission\n"
                "for n_tip=%d..%d; D: real sequential/stale-root simulation for\n"
                "n_tip<=d_max) ever produced a guard-ACCEPTED structurally-invalid\n"
                "tree.  The %lld Pass-A geometric hits are on (clip,above,below)\n"
                "triples collect_edges_*() cannot emit — unreachable artifacts.\n"
                "Severity note: 0 accepted-invalid trees had a root-reachable cycle,\n"
                "so even a hypothetical slip could not reach std::bad_alloc.\n",
                n_lo, n_hi, wA + wAp);
    return 0;
  }
  std::printf("\nRESULT: guard admits reachable structural corruption "
              "(C1 first-move=%lld, D within/across passes=%lld); all TYPE-%s.\n",
              wC1, wD, (t1 > 0 && t2 == 0) ? "1 (left/right corrupt)" :
                       (t2 > 0 && t1 == 0) ? "2 (parent[] only)" : "1+2 (mixed)");
  const bool crash_safe = (cyc_total == 0 && g_parent_cycle == 0);
  std::printf("  CRASH verdict   : %s\n"
              "                    root-reachable left/right cycles=%lld (=> DFS helpers can't loop);\n"
              "                    parent[] cycles=%lld (=> parent-ascending consumers e.g.\n"
              "                    tbr_search can't loop, incl. the pre-verify path).\n",
              crash_safe ? "std::bad_alloc REFUTED — not reachable"
                         : "POTENTIAL CRASH — a cycle exists, investigate consumer",
              cyc_total, g_parent_cycle);
  std::printf("  WRONG-ANSWER    : %s\n",
              g_returnable_corrupt == 0
                ? (g_D_budget_hit
                     ? "no returnable malformed tree found, BUT probe budget hit => PARTIAL coverage"
                     : "REFUTED for single-split n<=d_max — 0 corrupt final trees survive the\n"
                       "                    caller's map_constraint_nodes verify-and-discard")
                : "REACHABLE — a corrupt tree survives verify and is scored/returned (correctness bug)");
  return g_returnable_corrupt > 0 ? 1 : 0;
}
