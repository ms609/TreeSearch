#include "ts_fitch.h"
#include "ts_hsj.h"
#include "ts_sankoff.h"
#include <vector>
#include <cassert>
#include <R.h>
#ifdef TS_AUDIT_PROBE
#include <chrono>
#include <cstdio>
#endif

namespace ts {

// --- Node-level helpers ---

int fitch_downpass_node(
    const uint64_t* left_state,
    const uint64_t* right_state,
    uint64_t* node_state,
    int n_states,
    uint64_t active_mask)
{
  if (active_mask == 0) return 0;

  // Pass 1: compute any_intersect = OR( left[s] & right[s] )
  uint64_t any_intersect = simd::any_hit_reduce(left_state, right_state,
                                                 n_states);

  uint64_t needs_union = ~any_intersect & active_mask;
  int steps = popcount64(needs_union);

  // Pass 2: compute output states with broadcast masks
  simd::fitch_combine(left_state, right_state, node_state, n_states,
                      any_intersect, needs_union);

  return steps;
}

// Compute final states at node from ancestor's final and node's prelim.
// Returns true if any final state word changed.
static bool uppass_node(TreeState& tree, const DataSet& ds, int node) {
  int anc = tree.parent[node];
  size_t node_base = static_cast<size_t>(node) * tree.total_words;
  size_t anc_base = static_cast<size_t>(anc) * tree.total_words;
  bool changed = false;

  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    if (blk.active_mask == 0) continue;
    int offset = ds.block_word_offset[b];

    uint64_t any_intersect = simd::any_hit_reduce(
        &tree.final_[anc_base + offset],
        &tree.prelim[node_base + offset],
        blk.n_states);
    uint64_t has_isect = any_intersect;
    uint64_t no_isect = ~any_intersect & blk.active_mask;

    for (int s = 0; s < blk.n_states; ++s) {
      uint64_t isect = tree.final_[anc_base + offset + s]
                     & tree.prelim[node_base + offset + s];
      uint64_t new_val = (isect & has_isect)
                       | (tree.prelim[node_base + offset + s] & no_isect);
      if (new_val != tree.final_[node_base + offset + s]) changed = true;
      tree.final_[node_base + offset + s] = new_val;
    }
  }

  return changed;
}

// --- Full passes ---

int fitch_downpass(TreeState& tree, const DataSet& ds) {
  int total_steps = 0;

  for (int node : tree.postorder) {
    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];

    for (int b = 0; b < ds.n_blocks; ++b) {
      const CharBlock& blk = ds.blocks[b];
      if (blk.active_mask == 0) continue;
      int offset = ds.block_word_offset[b];

      const uint64_t* left_state =
          &tree.prelim[static_cast<size_t>(lc) * tree.total_words + offset];
      const uint64_t* right_state =
          &tree.prelim[static_cast<size_t>(rc) * tree.total_words + offset];
      uint64_t* node_state =
          &tree.prelim[static_cast<size_t>(node) * tree.total_words + offset];

      uint64_t any_intersect = simd::any_hit_reduce(
          left_state, right_state, blk.n_states);

      uint64_t needs_union = ~any_intersect & blk.active_mask;
      int nu = popcount64(needs_union);
      if (blk.upweight_mask) nu += popcount64(needs_union & blk.upweight_mask);
      total_steps += blk.weight * nu;

      // Store local cost
      tree.local_cost[static_cast<size_t>(node) * tree.n_blocks + b] =
          needs_union;

      // Compute output states with broadcast masks
      simd::fitch_combine(left_state, right_state, node_state,
                          blk.n_states, any_intersect, needs_union);
    }
  }

  return total_steps;
}

void fitch_uppass(TreeState& tree, const DataSet& ds) {
  // Root: final = prelim
  int root = tree.n_tip;
  size_t root_base = static_cast<size_t>(root) * tree.total_words;
  for (int w = 0; w < tree.total_words; ++w) {
    tree.final_[root_base + w] = tree.prelim[root_base + w];
  }

  // Reverse postorder (root to leaves) — skip root itself
  for (int i = static_cast<int>(tree.postorder.size()) - 1; i >= 0; --i) {
    int node = tree.postorder[i];
    if (node == root) continue;
    uppass_node(tree, ds, node);
  }
}

int fitch_score(TreeState& tree, const DataSet& ds) {
  int score = fitch_downpass(tree, ds);
  fitch_uppass(tree, ds);
  return score;
}

// --- Incremental passes for SPR clipping ---

int fitch_incremental_downpass(TreeState& tree, const DataSet& ds,
                               int start_node,
                               std::vector<int>* cs_delta) {
  int length_delta = 0;
  int node = start_node;
  if (cs_delta) std::fill(cs_delta->begin(), cs_delta->end(), 0);

  while (true) {
    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];

    // Save old state before overwriting
    tree.save_node_state(node);

    bool changed = false;

    for (int b = 0; b < ds.n_blocks; ++b) {
      const CharBlock& blk = ds.blocks[b];
      if (blk.active_mask == 0) continue;
      int offset = ds.block_word_offset[b];

      const uint64_t* left_state =
          &tree.prelim[static_cast<size_t>(lc) * tree.total_words + offset];
      const uint64_t* right_state =
          &tree.prelim[static_cast<size_t>(rc) * tree.total_words + offset];
      uint64_t* node_state =
          &tree.prelim[static_cast<size_t>(node) * tree.total_words + offset];

      // Subtract old local cost
      size_t cost_idx = static_cast<size_t>(node) * tree.n_blocks + b;
      uint64_t old_cost = tree.local_cost[cost_idx];
      int old_nu = popcount64(old_cost);
      if (blk.upweight_mask) old_nu += popcount64(old_cost & blk.upweight_mask);
      length_delta -= blk.weight * old_nu;

      // Compute new prelim
      uint64_t any_intersect = simd::any_hit_reduce(
          left_state, right_state, blk.n_states);

      uint64_t needs_union = ~any_intersect & blk.active_mask;
      int new_nu = popcount64(needs_union);
      if (blk.upweight_mask) new_nu += popcount64(needs_union & blk.upweight_mask);
      length_delta += blk.weight * new_nu;

      // Store new local cost
      tree.local_cost[cost_idx] = needs_union;

      // IW dirty-region: accumulate per-pattern char_steps change for this
      // node (old bits removed, new bits added). Matches extract_char_steps'
      // raw-local_cost counting (already active-masked).
      if (cs_delta) {
        uint64_t oc = old_cost;
        while (oc) { int c = ctz64(oc); (*cs_delta)[blk.pattern_index[c]]--; oc &= oc - 1; }
        uint64_t nc = needs_union;
        while (nc) { int c = ctz64(nc); (*cs_delta)[blk.pattern_index[c]]++; nc &= nc - 1; }
      }

      for (int s = 0; s < blk.n_states; ++s) {
        uint64_t isect = left_state[s] & right_state[s];
        uint64_t uni = left_state[s] | right_state[s];
        uint64_t new_val = (isect & any_intersect) | (uni & needs_union);
        if (new_val != node_state[s]) changed = true;
        node_state[s] = new_val;
      }
    }

    // Stop if at root or prelim didn't change
    if (!changed || node == tree.n_tip) break;
    int p = tree.parent[node];
    if (p == node) break;
    node = p;
  }

  return length_delta;
}

void fitch_incremental_uppass(TreeState& tree, const DataSet& ds,
                              int clip_ancestor) {
  // Recompute final_ for nodes affected by the clip.
  // Start at root, propagate downward through changed nodes.
  // Only nodes whose ancestor's final changed need recomputation.

  // Root: final = prelim (may have changed during incremental downpass)
  int root = tree.n_tip;
  size_t root_base = static_cast<size_t>(root) * tree.total_words;
  bool root_changed = false;
  for (int w = 0; w < tree.total_words; ++w) {
    if (tree.final_[root_base + w] != tree.prelim[root_base + w]) {
      root_changed = true;
    }
    tree.final_[root_base + w] = tree.prelim[root_base + w];
  }

  if (!root_changed) {
    // Check if any node on the path had its prelim changed but root didn't.
    // The incremental downpass may have stopped before reaching root.
    // In that case, final_ only needs updating for the sister of the clip
    // point (whose new ancestor may have different final states).
  }

  // Use reverse postorder, but only visit nodes whose ancestor's final
  // may have changed. We track this with a "dirty" flag per node.
  // Reusable per-thread scratch (S-PROF round 3 / Tier 1): this function runs
  // once per clip in the TBR hot loop, so a fresh vector<bool> here was a
  // per-clip heap allocation. thread_local keeps it per-thread-safe (each
  // search thread owns its TreeState); char avoids vector<bool> proxy-bit
  // access in the reverse scan below. assign() reuses capacity after warmup.
  std::vector<char> dirty;
  dirty.assign(tree.n_node, 0);

  // Mark root as dirty (we just updated it; its children need checking)
  dirty[root] = root_changed;

  // Also mark the clip ancestor — its children definitely need checking
  // because the topology changed around it.
  if (clip_ancestor >= tree.n_tip) {
    dirty[clip_ancestor] = true;
  }

  // Reverse postorder traversal
  for (int i = static_cast<int>(tree.postorder.size()) - 1; i >= 0; --i) {
    int node = tree.postorder[i];
    if (node == root) continue;

    int anc = tree.parent[node];
    if (!dirty[anc]) continue;

    // Save state before modifying (if not already saved during downpass)
    // The downpass only saved nodes on the rootward path; the uppass may
    // touch additional nodes (siblings and their descendants).
    tree.save_node_state(node);

    bool changed = uppass_node(tree, ds, node);

    if (changed && node >= tree.n_tip) {
      dirty[node] = true;
    }
  }
}

// --- Dirty-set rescore (T-300) ---
//
// After an SPR move, the set of nodes whose prelim is stale is the union of
// the rootward paths from the two clip endpoints.  Walking each chain
// independently (as the prior implementation did) double-reads shared
// ancestors and is brittle in the IW path.  The dirty-set approach visits
// each affected node exactly once in postorder, reading current children's
// prelims — which are guaranteed correct because postorder processes
// children before parents.

int fitch_dirty_downpass(TreeState& tree, const DataSet& ds,
                         int start_a, int start_b) {
  std::vector<char> dirty(tree.n_node, 0);

  // Mark the rootward path from `node` up to (and including) the root.
  // Mirrors the parent-self termination used elsewhere in this file.
  auto mark_path = [&](int node) {
    while (node >= tree.n_tip && !dirty[node]) {
      dirty[node] = 1;
      int p = tree.parent[node];
      if (p == node) break;  // root
      node = p;
    }
  };
  mark_path(start_a);
  mark_path(start_b);

  int length_delta = 0;

  for (int node : tree.postorder) {
    if (node < tree.n_tip) continue;
    if (!dirty[node]) continue;

    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];

    tree.save_node_state(node);

    for (int b = 0; b < ds.n_blocks; ++b) {
      const CharBlock& blk = ds.blocks[b];
      if (blk.active_mask == 0) continue;
      int offset = ds.block_word_offset[b];

      const uint64_t* left_state =
          &tree.prelim[static_cast<size_t>(lc) * tree.total_words + offset];
      const uint64_t* right_state =
          &tree.prelim[static_cast<size_t>(rc) * tree.total_words + offset];
      uint64_t* node_state =
          &tree.prelim[static_cast<size_t>(node) * tree.total_words + offset];

      size_t cost_idx = static_cast<size_t>(node) * tree.n_blocks + b;
      uint64_t old_cost = tree.local_cost[cost_idx];
      int old_nu = popcount64(old_cost);
      if (blk.upweight_mask) old_nu += popcount64(old_cost & blk.upweight_mask);
      length_delta -= blk.weight * old_nu;

      uint64_t any_intersect = simd::any_hit_reduce(
          left_state, right_state, blk.n_states);
      uint64_t needs_union = ~any_intersect & blk.active_mask;
      int new_nu = popcount64(needs_union);
      if (blk.upweight_mask) new_nu += popcount64(needs_union & blk.upweight_mask);
      length_delta += blk.weight * new_nu;

      tree.local_cost[cost_idx] = needs_union;

      simd::fitch_combine(left_state, right_state, node_state,
                          blk.n_states, any_intersect, needs_union);
    }
  }

  return length_delta;
}

void fitch_dirty_uppass(TreeState& tree, const DataSet& ds,
                        int start_a, int start_b) {
  // Step 1: root final_ = prelim (root prelim may have changed in downpass).
  int root = tree.n_tip;
  size_t root_base = static_cast<size_t>(root) * tree.total_words;
  for (int w = 0; w < tree.total_words; ++w) {
    tree.final_[root_base + w] = tree.prelim[root_base + w];
  }

  // Step 2: mark dirty_up = same rootward paths as the downpass.  Every node
  // on these paths had its prelim updated, so its final_ may shift and its
  // children must be re-checked.
  std::vector<char> dirty_up(tree.n_node, 0);
  auto mark_path = [&](int node) {
    while (node >= tree.n_tip && !dirty_up[node]) {
      dirty_up[node] = 1;
      int p = tree.parent[node];
      if (p == node) break;
      node = p;
    }
  };
  mark_path(start_a);
  mark_path(start_b);

  // Step 3: reverse postorder — visit any node whose parent is dirty_up.
  // If that node's final_ changes, propagate the flag to it.
  for (int i = static_cast<int>(tree.postorder.size()) - 1; i >= 0; --i) {
    int node = tree.postorder[i];
    if (node == root) continue;

    int anc = tree.parent[node];
    if (!dirty_up[anc]) continue;

    tree.save_node_state(node);
    bool changed = uppass_node(tree, ds, node);
    if (changed && node >= tree.n_tip) {
      dirty_up[node] = 1;
    }
  }
}

// --- Indirect tree length calculation (Goloboff 1996) ---

int fitch_indirect_length(const uint64_t* clip_prelim,
                          const TreeState& tree,
                          const DataSet& ds,
                          int node_a, int node_d) {
  // Approximate the virtual root Y at edge (A, D) as the union of the endpoint
  // final states:
  //   Y = final(A) | final(D)
  // Extra steps = count of characters where clip_prelim & Y == 0.
  //
  // NOTE: this is an APPROXIMATION, not exact.  The union of the two endpoints'
  // final sets is a superset of the true directional Fitch edge set, so it makes
  // more states appear "available" on the edge than any most-parsimonious
  // reconstruction allows — hence it UNDER-counts the true insertion cost (never
  // over-counts).  The exact cost uses the directional edge set
  //   edge_set[D] = combine(prelim[D], up[D])   (per-character intersect-else-union)
  // via compute_insertion_edge_sets() + fitch_indirect_length_cached(); see
  // ts_fitch.h.  This cheaper union variant is retained for callers (temper) that
  // rank candidates approximately and then re-score the chosen move exactly.

  int extra_steps = 0;

  size_t a_base = static_cast<size_t>(node_a) * tree.total_words;
  size_t d_base = static_cast<size_t>(node_d) * tree.total_words;

  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    if (blk.active_mask == 0) continue;
    int offset = ds.block_word_offset[b];

    uint64_t any_hit = simd::any_hit_reduce3(
        &clip_prelim[offset],
        &tree.final_[a_base + offset],
        &tree.final_[d_base + offset],
        blk.n_states);

    uint64_t needs_step = ~any_hit & blk.active_mask;
    int ns = popcount64(needs_step);
    if (blk.upweight_mask) ns += popcount64(needs_step & blk.upweight_mask);
    extra_steps += blk.weight * ns;
  }

  return extra_steps;
}


// Early-termination variant of fitch_indirect_length.
int fitch_indirect_length_bounded(const uint64_t* clip_prelim,
                                  const TreeState& tree,
                                  const DataSet& ds,
                                  int node_a, int node_d,
                                  int cutoff) {
  int extra_steps = 0;

  size_t a_base = static_cast<size_t>(node_a) * tree.total_words;
  size_t d_base = static_cast<size_t>(node_d) * tree.total_words;

  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    if (blk.active_mask == 0) continue;
    int offset = ds.block_word_offset[b];

    uint64_t any_hit = simd::any_hit_reduce3(
        &clip_prelim[offset],
        &tree.final_[a_base + offset],
        &tree.final_[d_base + offset],
        blk.n_states);

    uint64_t needs_step = ~any_hit & blk.active_mask;
    int ns = popcount64(needs_step);
    if (blk.upweight_mask) ns += popcount64(needs_step & blk.upweight_mask);
    extra_steps += blk.weight * ns;
    if (extra_steps >= cutoff) return extra_steps;
  }

  return extra_steps;
}

// Precomputed-vroot variant with early termination.
int fitch_indirect_length_cached(const uint64_t* clip_prelim,
                                 const uint64_t* vroot,
                                 const DataSet& ds,
                                 int cutoff) {
  int extra_steps = 0;

  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    if (blk.active_mask == 0) continue;
    int offset = ds.block_word_offset[b];

    uint64_t any_hit = simd::any_hit_reduce(
        &clip_prelim[offset], &vroot[offset], blk.n_states);

    uint64_t needs_step = ~any_hit & blk.active_mask;
    int ns = popcount64(needs_step);
    if (blk.upweight_mask) ns += popcount64(needs_step & blk.upweight_mask);
    extra_steps += blk.weight * ns;
    if (extra_steps >= cutoff) return extra_steps;
  }

  return extra_steps;
}

// Exact per-node insertion edge sets via directional Fitch messages.
// See the header for the formula.  O(n * chars): one preorder up-pass plus one
// combine per node.
#ifdef TS_AUDIT_PROBE
// Audit #56: time the no-bail precompute to measure its share of SECTOR wall
// (run via ts_rss_search, where every call is a sector call). precompute_share
// x 0.40 (measured block reduction) x 0.30 (sectorial mission share) = the
// realizable mission saving from column reduction. RAII so all exits are timed.
static long long g_precompute_ns = 0;
static long long g_precompute_calls = 0;
struct PrecomputeTimer {
  std::chrono::steady_clock::time_point t0;
  PrecomputeTimer() : t0(std::chrono::steady_clock::now()) {}
  ~PrecomputeTimer() {
    g_precompute_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now() - t0).count();
    if ((++g_precompute_calls % 50000LL) == 0)
      std::fprintf(stderr, "PRECOMPUTE_NS calls=%lld ns=%lld\n",
                   g_precompute_calls, g_precompute_ns);
  }
};
#endif
void compute_insertion_edge_sets(const TreeState& tree, const DataSet& ds,
                                 std::vector<uint64_t>& edge_set,
                                 std::vector<uint64_t>& up,
                                 std::vector<int>& pre) {
#ifdef TS_AUDIT_PROBE
  PrecomputeTimer _pt;
#endif
  const int n_tip = tree.n_tip;
  const int tw    = tree.total_words;
  const int nb    = ds.n_blocks;
  const int root  = n_tip;

  // Non-zeroing size-ensure on caller-owned scratch.  `up` and `edge_set` grow
  // monotonically across calls, so after the first call no zero-fill happens
  // (resize value-inits only NEW elements).  Every slot a downstream reader
  // touches is edge_set[D] for a non-root in-tree node D, and the two combine
  // loops below overwrite exactly those slots before any read; the stale
  // contents of grown-but-unwritten slots (the root slot, and slots for
  // clipped-out nodes that are not edges of the current tree) are never
  // observed.  This removes the per-call assign() zero-fill and the per-call
  // up/pre heap allocations that VTune flagged as ~27% of EW Fitch CPU.
  const size_t N = static_cast<size_t>(tree.n_node) * tw;
  if (edge_set.size() < N) edge_set.resize(N);
  if (up.size() < N) up.resize(N);

  // Preorder over current in-tree nodes (parents before children).
  pre.clear();
  {
    std::vector<int> st;
    st.push_back(root);
    while (!st.empty()) {
      int nd = st.back(); st.pop_back();
      pre.push_back(nd);
      if (nd >= n_tip) {
        int ni = nd - n_tip;
        st.push_back(tree.left[ni]);
        st.push_back(tree.right[ni]);
      }
    }
  }

  // Fitch combine (per character intersect-else-union) of a & b into dst.
  auto combine = [&](uint64_t* dst, const uint64_t* a, const uint64_t* b) {
    for (int bi = 0; bi < nb; ++bi) {
      const CharBlock& blk = ds.blocks[bi];
      int off = ds.block_word_offset[bi];
      uint64_t any_isect = 0;
      for (int s = 0; s < blk.n_states; ++s) any_isect |= a[off + s] & b[off + s];
      uint64_t needs_union = ~any_isect & blk.active_mask;
      for (int s = 0; s < blk.n_states; ++s)
        dst[off + s] = ((a[off + s] & b[off + s]) & any_isect)
                     | ((a[off + s] | b[off + s]) & needs_union);
    }
  };

  // Directional up-pass: up[D] = combine(up[parent], prelim[sibling]);
  // root is a degree-2 vertex so up[child] = prelim[other child].
  for (int D : pre) {
    if (D == root) continue;
    int A   = tree.parent[D];
    int ai  = A - n_tip;
    int Sib = (tree.left[ai] == D) ? tree.right[ai] : tree.left[ai];
    uint64_t* uD       = &up[static_cast<size_t>(D) * tw];
    const uint64_t* pS = &tree.prelim[static_cast<size_t>(Sib) * tw];
    if (A == root) {
      for (int w = 0; w < tw; ++w) uD[w] = pS[w];
    } else {
      combine(uD, &up[static_cast<size_t>(A) * tw], pS);
    }
  }

  // Edge set above each non-root node: E[D] = combine(prelim[D], up[D]).
#ifndef NDEBUG
  // Debug-only write-before-read guard: the non-zeroing size-ensure relies on
  // every reader index (edge_set[D] for a non-root in-tree node D) being
  // overwritten this call.  Record which non-root nodes are written and assert
  // completeness against the in-tree node set, so a stale slot can never be
  // read.  Release builds (NDEBUG) compile none of this.
  std::vector<char> written(static_cast<size_t>(tree.n_node), 0);
#endif
  for (int D : pre) {
    if (D == root) continue;
    combine(&edge_set[static_cast<size_t>(D) * tw],
            &tree.prelim[static_cast<size_t>(D) * tw],
            &up[static_cast<size_t>(D) * tw]);
#ifndef NDEBUG
    written[static_cast<size_t>(D)] = 1;
#endif
  }
#ifndef NDEBUG
  // Every in-tree node except the root must have its edge_set slot written.
  for (int D : pre) {
    if (D == root) continue;
    assert(written[static_cast<size_t>(D)] &&
           "compute_insertion_edge_sets: in-tree node left unwritten");
  }
#endif
}


// --- Flat EW specializations ---
//
// These eliminate per-block overhead from the indirect scoring hot path:
// - No CharBlock struct dereference (288 bytes each; FlatBlock is 16 bytes)
// - No upweight_mask check (always 0 during normal search)
// - No weight multiply (always 1 when all_weight_one)
// - No active_mask==0 check (empty blocks never exist after build_dataset)

int fitch_indirect_bounded_flat(const uint64_t* clip_prelim,
                                const TreeState& tree,
                                const DataSet& ds,
                                int node_a, int node_d,
                                int cutoff) {
  int extra_steps = 0;
  const FlatBlock* fb = ds.flat_blocks.data();
  size_t a_base = static_cast<size_t>(node_a) * tree.total_words;
  size_t d_base = static_cast<size_t>(node_d) * tree.total_words;

  for (int b = 0; b < ds.n_blocks; ++b) {
    uint64_t any_hit = simd::any_hit_reduce3(
        &clip_prelim[fb[b].offset],
        &tree.final_[a_base + fb[b].offset],
        &tree.final_[d_base + fb[b].offset],
        fb[b].n_states);
    extra_steps += popcount64(~any_hit & fb[b].active_mask);
    if (extra_steps >= cutoff) return extra_steps;
  }
  return extra_steps;
}

int fitch_indirect_cached_flat(const uint64_t* clip_prelim,
                               const uint64_t* vroot,
                               const DataSet& ds,
                               int cutoff) {
  int extra_steps = 0;
  const FlatBlock* fb = ds.flat_blocks.data();

  for (int b = 0; b < ds.n_blocks; ++b) {
    uint64_t any_hit = simd::any_hit_reduce(
        &clip_prelim[fb[b].offset], &vroot[fb[b].offset],
        fb[b].n_states);
    extra_steps += popcount64(~any_hit & fb[b].active_mask);
    if (extra_steps >= cutoff) return extra_steps;
  }
  return extra_steps;
}

// NA-aware flat bounded indirect (SPR candidates).
// Handles mixed standard + inapplicable blocks using FlatBlock metadata.
int fitch_na_indirect_bounded_flat(const uint64_t* clip_prelim,
                                   const uint64_t* clip_actives,
                                   const TreeState& tree,
                                   const DataSet& ds,
                                   int node_a, int node_d,
                                   int cutoff) {
  int extra_steps = 0;
  const FlatBlock* fb = ds.flat_blocks.data();
  size_t a_base = static_cast<size_t>(node_a) * tree.total_words;
  size_t d_base = static_cast<size_t>(node_d) * tree.total_words;

  for (int b = 0; b < ds.n_blocks; ++b) {
    uint64_t needs_step;
    if (!fb[b].has_inapplicable) {
      uint64_t any_hit = simd::any_hit_reduce3(
          &clip_prelim[fb[b].offset],
          &tree.final_[a_base + fb[b].offset],
          &tree.final_[d_base + fb[b].offset],
          fb[b].n_states);
      needs_step = ~any_hit & fb[b].active_mask;
    } else {
      uint64_t any_hit = simd::any_hit_reduce3_from1(
          &clip_prelim[fb[b].offset],
          &tree.final_[a_base + fb[b].offset],
          &tree.final_[d_base + fb[b].offset],
          fb[b].n_states);
      uint64_t clip_has_active =
          simd::or_reduce(&clip_actives[fb[b].offset], fb[b].n_states, 1);
      uint64_t below_has_active =
          simd::or_reduce(&tree.subtree_actives[d_base + fb[b].offset],
                          fb[b].n_states, 1);
      needs_step = ~any_hit & clip_has_active & below_has_active
                 & fb[b].active_mask;
    }
    extra_steps += popcount64(needs_step);
    if (extra_steps >= cutoff) return extra_steps;
  }
  return extra_steps;
}

// NA-aware flat cached indirect (TBR rerooting candidates).
int fitch_na_indirect_cached_flat(const uint64_t* clip_prelim,
                                  const uint64_t* clip_actives,
                                  const uint64_t* vroot,
                                  const uint64_t* below_actives,
                                  const DataSet& ds,
                                  int cutoff) {
  int extra_steps = 0;
  const FlatBlock* fb = ds.flat_blocks.data();

  for (int b = 0; b < ds.n_blocks; ++b) {
    uint64_t needs_step;
    if (!fb[b].has_inapplicable) {
      uint64_t any_hit = simd::any_hit_reduce(
          &clip_prelim[fb[b].offset], &vroot[fb[b].offset],
          fb[b].n_states);
      needs_step = ~any_hit & fb[b].active_mask;
    } else {
      uint64_t any_hit = simd::any_hit_reduce_from1(
          &clip_prelim[fb[b].offset], &vroot[fb[b].offset],
          fb[b].n_states);
      uint64_t clip_has_active =
          simd::or_reduce(&clip_actives[fb[b].offset], fb[b].n_states, 1);
      needs_step = ~any_hit & clip_has_active & below_actives[b]
                 & fb[b].active_mask;
    }
    extra_steps += popcount64(needs_step);
    if (extra_steps >= cutoff) return extra_steps;
  }
  return extra_steps;
}

// --- 4-wide TBR rerooting batch functions (T-245) ---
//
// Both functions process all blocks for all 4 candidates in lockstep.
// Within each block the 4 simd::any_hit_reduce() calls are data-independent,
// so the out-of-order CPU can issue 4 separate load streams concurrently,
// hiding L2 latency for the vroot_cache array.
//
// Early-exit: once ALL 4 accumulators exceed cutoff we stop iterating blocks.
// Using '&' (bitwise AND) rather than '&&' in the combined cutoff check avoids
// branch-prediction overhead on the hot path — all 4 comparisons are always
// evaluated and combined into a single bitmask test.

#ifdef TS_AUDIT_PROBE
// Audit #57: x4 reroot-batch wasted-block counter. The x4 scans all four members
// to the DEEPEST-bailing member's depth (breaks only when ALL four exceed cutoff).
// "wasted" = per-member blocks scanned AFTER that member individually crossed the
// cutoff; frac = wasted / total-scanned. Measures the ceiling for a force-scalar
// reroot (ILP-confounded, so a large frac still needs a wall A/B to settle sign).
static long long g_x4_waste = 0;
static long long g_x4_total = 0;
static unsigned long long g_x4_calls = 0;
#endif
void fitch_indirect_cached_flat_x4(
    const uint64_t* clip_prelim,
    const uint64_t* vroot0, const uint64_t* vroot1,
    const uint64_t* vroot2, const uint64_t* vroot3,
    const DataSet& ds, int cutoff, int out[4]) {
  const FlatBlock* fb = ds.flat_blocks.data();
  int es0 = 0, es1 = 0, es2 = 0, es3 = 0;
#ifdef TS_AUDIT_PROBE
  int bail[4] = {-1, -1, -1, -1};
  int nb = 0;
#endif

  for (int b = 0; b < ds.n_blocks; ++b) {
    const int  off  = fb[b].offset;
    const int  nst  = fb[b].n_states;
    const uint64_t mask = fb[b].active_mask;

    // 4 independent loads from distinct vroot_cache rows.
    uint64_t a0 = simd::any_hit_reduce(&clip_prelim[off], &vroot0[off], nst);
    uint64_t a1 = simd::any_hit_reduce(&clip_prelim[off], &vroot1[off], nst);
    uint64_t a2 = simd::any_hit_reduce(&clip_prelim[off], &vroot2[off], nst);
    uint64_t a3 = simd::any_hit_reduce(&clip_prelim[off], &vroot3[off], nst);

    es0 += popcount64(~a0 & mask);
    es1 += popcount64(~a1 & mask);
    es2 += popcount64(~a2 & mask);
    es3 += popcount64(~a3 & mask);
#ifdef TS_AUDIT_PROBE
    nb = b + 1;
    if (bail[0] < 0 && es0 >= cutoff) bail[0] = nb;
    if (bail[1] < 0 && es1 >= cutoff) bail[1] = nb;
    if (bail[2] < 0 && es2 >= cutoff) bail[2] = nb;
    if (bail[3] < 0 && es3 >= cutoff) bail[3] = nb;
#endif

    // Bitwise & avoids short-circuit, keeping branch count low.
    if ((es0 >= cutoff) & (es1 >= cutoff) & (es2 >= cutoff) & (es3 >= cutoff))
      break;
  }
#ifdef TS_AUDIT_PROBE
  {
    int B = (nb > 0) ? nb : ds.n_blocks;
    for (int i = 0; i < 4; ++i) {
      int bi = (bail[i] < 0) ? B : bail[i];
      g_x4_waste += (B - bi);
    }
    g_x4_total += 4LL * B;
    if ((++g_x4_calls % 5000000ULL) == 0) {
      std::fprintf(stderr, "X4_WASTE calls=%llu wasted=%lld total=%lld frac=%.4f\n",
                   (unsigned long long)g_x4_calls, (long long)g_x4_waste,
                   (long long)g_x4_total,
                   g_x4_total ? (double)g_x4_waste / (double)g_x4_total : 0.0);
    }
  }
#endif

  out[0] = es0; out[1] = es1; out[2] = es2; out[3] = es3;
}

void fitch_na_indirect_cached_flat_x4(
    const uint64_t* clip_prelim,
    const uint64_t* clip_actives,
    const uint64_t* vroot0, const uint64_t* vroot1,
    const uint64_t* vroot2, const uint64_t* vroot3,
    const uint64_t* ba0, const uint64_t* ba1,
    const uint64_t* ba2, const uint64_t* ba3,
    const DataSet& ds, int cutoff, int out[4]) {
  const FlatBlock* fb = ds.flat_blocks.data();
  int es0 = 0, es1 = 0, es2 = 0, es3 = 0;

  for (int b = 0; b < ds.n_blocks; ++b) {
    const int  off  = fb[b].offset;
    const int  nst  = fb[b].n_states;
    const uint64_t mask = fb[b].active_mask;

    uint64_t ns0, ns1, ns2, ns3;

    if (!fb[b].has_inapplicable) {
      // Standard block: plain any_hit_reduce (skip inapplicable state 0).
      uint64_t a0 = simd::any_hit_reduce(&clip_prelim[off], &vroot0[off], nst);
      uint64_t a1 = simd::any_hit_reduce(&clip_prelim[off], &vroot1[off], nst);
      uint64_t a2 = simd::any_hit_reduce(&clip_prelim[off], &vroot2[off], nst);
      uint64_t a3 = simd::any_hit_reduce(&clip_prelim[off], &vroot3[off], nst);
      ns0 = ~a0 & mask;
      ns1 = ~a1 & mask;
      ns2 = ~a2 & mask;
      ns3 = ~a3 & mask;
    } else {
      // NA block: from1 skips the inapplicable state; AND with active masks.
      // clip_has_active is shared across all 4 candidates for this block.
      uint64_t clip_ha = simd::or_reduce(&clip_actives[off], nst, 1);
      uint64_t a0 = simd::any_hit_reduce_from1(&clip_prelim[off], &vroot0[off], nst);
      uint64_t a1 = simd::any_hit_reduce_from1(&clip_prelim[off], &vroot1[off], nst);
      uint64_t a2 = simd::any_hit_reduce_from1(&clip_prelim[off], &vroot2[off], nst);
      uint64_t a3 = simd::any_hit_reduce_from1(&clip_prelim[off], &vroot3[off], nst);
      ns0 = ~a0 & clip_ha & ba0[b] & mask;
      ns1 = ~a1 & clip_ha & ba1[b] & mask;
      ns2 = ~a2 & clip_ha & ba2[b] & mask;
      ns3 = ~a3 & clip_ha & ba3[b] & mask;
    }

    es0 += popcount64(ns0);
    es1 += popcount64(ns1);
    es2 += popcount64(ns2);
    es3 += popcount64(ns3);

    if ((es0 >= cutoff) & (es1 >= cutoff) & (es2 >= cutoff) & (es3 >= cutoff))
      break;
  }

  out[0] = es0; out[1] = es1; out[2] = es2; out[3] = es3;
}

// 4-wide IW batch (pure-IW TBR rerooting). Mirrors fitch_indirect_cached_flat_x4
// but accumulates weighted iw_delta via a per-candidate ctz gather instead of a
// popcount. 4 data-independent any_hit_reduce streams hide vroot_cache L2 latency
// (the same ILP win T-245 gave EW); each es{k} keeps the scalar add order of
// indirect_iw_length_cached, so per-candidate results are bit-identical.
void indirect_iw_cached_flat_x4(
    const uint64_t* clip_prelim,
    const uint64_t* vroot0, const uint64_t* vroot1,
    const uint64_t* vroot2, const uint64_t* vroot3,
    const DataSet& ds, double base_iw,
    const std::vector<double>& iw_delta,
    double cutoff, double out[4]) {
  double es0 = base_iw, es1 = base_iw, es2 = base_iw, es3 = base_iw;

  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    if (blk.active_mask == 0) continue;
    const int off = ds.block_word_offset[b];
    const int nst = blk.n_states;

    // 4 independent loads from distinct vroot_cache rows.
    uint64_t a0 = simd::any_hit_reduce(&clip_prelim[off], &vroot0[off], nst);
    uint64_t a1 = simd::any_hit_reduce(&clip_prelim[off], &vroot1[off], nst);
    uint64_t a2 = simd::any_hit_reduce(&clip_prelim[off], &vroot2[off], nst);
    uint64_t a3 = simd::any_hit_reduce(&clip_prelim[off], &vroot3[off], nst);

    uint64_t ns0 = ~a0 & blk.active_mask;
    uint64_t ns1 = ~a1 & blk.active_mask;
    uint64_t ns2 = ~a2 & blk.active_mask;
    uint64_t ns3 = ~a3 & blk.active_mask;

    while (ns0) { int c = ctz64(ns0); es0 += iw_delta[blk.pattern_index[c]]; ns0 &= ns0 - 1; }
    while (ns1) { int c = ctz64(ns1); es1 += iw_delta[blk.pattern_index[c]]; ns1 &= ns1 - 1; }
    while (ns2) { int c = ctz64(ns2); es2 += iw_delta[blk.pattern_index[c]]; ns2 &= ns2 - 1; }
    while (ns3) { int c = ctz64(ns3); es3 += iw_delta[blk.pattern_index[c]]; ns3 &= ns3 - 1; }

    if ((es0 >= cutoff) & (es1 >= cutoff) & (es2 >= cutoff) & (es3 >= cutoff))
      break;
  }

  out[0] = es0; out[1] = es1; out[2] = es2; out[3] = es3;
}

// 4-wide NA-IW batch (inapplicable characters + implied-weights TBR rerooting).
// Fuses the NA active-mask candidate selection of fitch_na_indirect_cached_flat_x4
// (standard vs has_inapplicable blocks; from1 reduce; shared clip_has_active;
// per-candidate below_actives AND) with the per-candidate iw_delta ctz-gather of
// indirect_iw_cached_flat_x4. Each es{k} keeps the exact scalar add order of
// indirect_na_iw_length_cached, so per-candidate results are bit-identical; the
// shared (all-4 >= cutoff) bail only changes early-exit on cutoff-LOSING
// candidates (their partial sum may then differ from the scalar early-return) --
// any sub-cutoff candidate the caller actually reads is fully accumulated, so the
// search outcome (argmin) is byte-identical. 4 data-independent any_hit_reduce
// streams over distinct vroot_cache rows hide the L2 load latency (the T-245 ILP
// win, now on the IW-NA path).
void indirect_na_iw_cached_flat_x4(
    const uint64_t* clip_prelim,
    const uint64_t* clip_actives,
    const uint64_t* vroot0, const uint64_t* vroot1,
    const uint64_t* vroot2, const uint64_t* vroot3,
    const uint64_t* ba0, const uint64_t* ba1,
    const uint64_t* ba2, const uint64_t* ba3,
    const DataSet& ds, double base_iw,
    const std::vector<double>& iw_delta,
    double cutoff, double out[4]) {
  double es0 = base_iw, es1 = base_iw, es2 = base_iw, es3 = base_iw;

  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    if (blk.active_mask == 0) continue;   // mirror scalar indirect_na_iw_length_cached
    const int off = ds.block_word_offset[b];
    const int k = blk.n_states;

    uint64_t ns0, ns1, ns2, ns3;

    if (!blk.has_inapplicable) {
      // Standard block: plain any_hit_reduce (states 0..k-1), invert & mask.
      uint64_t a0 = simd::any_hit_reduce(&clip_prelim[off], &vroot0[off], k);
      uint64_t a1 = simd::any_hit_reduce(&clip_prelim[off], &vroot1[off], k);
      uint64_t a2 = simd::any_hit_reduce(&clip_prelim[off], &vroot2[off], k);
      uint64_t a3 = simd::any_hit_reduce(&clip_prelim[off], &vroot3[off], k);
      ns0 = ~a0 & blk.active_mask;
      ns1 = ~a1 & blk.active_mask;
      ns2 = ~a2 & blk.active_mask;
      ns3 = ~a3 & blk.active_mask;
    } else {
      // NA block: from1 reduce skips the inapplicable state; clip_has_active is
      // shared across all 4 candidates; AND with each candidate's below_actives.
      uint64_t clip_ha = simd::or_reduce(&clip_actives[off], k, 1);
      uint64_t a0 = simd::any_hit_reduce_from1(&clip_prelim[off], &vroot0[off], k);
      uint64_t a1 = simd::any_hit_reduce_from1(&clip_prelim[off], &vroot1[off], k);
      uint64_t a2 = simd::any_hit_reduce_from1(&clip_prelim[off], &vroot2[off], k);
      uint64_t a3 = simd::any_hit_reduce_from1(&clip_prelim[off], &vroot3[off], k);
      ns0 = ~a0 & clip_ha & ba0[b] & blk.active_mask;
      ns1 = ~a1 & clip_ha & ba1[b] & blk.active_mask;
      ns2 = ~a2 & clip_ha & ba2[b] & blk.active_mask;
      ns3 = ~a3 & clip_ha & ba3[b] & blk.active_mask;
    }

    while (ns0) { int c = ctz64(ns0); es0 += iw_delta[blk.pattern_index[c]]; ns0 &= ns0 - 1; }
    while (ns1) { int c = ctz64(ns1); es1 += iw_delta[blk.pattern_index[c]]; ns1 &= ns1 - 1; }
    while (ns2) { int c = ctz64(ns2); es2 += iw_delta[blk.pattern_index[c]]; ns2 &= ns2 - 1; }
    while (ns3) { int c = ctz64(ns3); es3 += iw_delta[blk.pattern_index[c]]; ns3 &= ns3 - 1; }

    if ((es0 >= cutoff) & (es1 >= cutoff) & (es2 >= cutoff) & (es3 >= cutoff))
      break;
  }

  out[0] = es0; out[1] = es1; out[2] = es2; out[3] = es3;
}

// --- Per-character step extraction ---

void extract_char_steps(const TreeState& tree, const DataSet& ds,
                        std::vector<int>& char_steps) {
  std::fill(char_steps.begin(), char_steps.end(), 0);
  // Standard blocks: count from local_cost
  for (int node : tree.postorder) {
    for (int b = 0; b < ds.n_blocks; ++b) {
      const CharBlock& blk = ds.blocks[b];
      if (blk.has_inapplicable || blk.active_mask == 0) continue;
      uint64_t mask =
          tree.local_cost[static_cast<size_t>(node) * tree.n_blocks + b];
      while (mask) {
        int c = ctz64(mask);
        char_steps[blk.pattern_index[c]] += 1;
        mask &= mask - 1;
      }
    }
  }

  // NA blocks: use same formula as Pass 3 of three-pass algorithm.
  // Requires subtree_actives (computed by fitch_na_score).
  for (int node : tree.postorder) {
    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];
    size_t nb = static_cast<size_t>(node) * tree.total_words;
    size_t lb = static_cast<size_t>(lc) * tree.total_words;
    size_t rb = static_cast<size_t>(rc) * tree.total_words;

    for (int b = 0; b < ds.n_blocks; ++b) {
      if (!ds.blocks[b].has_inapplicable || ds.blocks[b].active_mask == 0)
        continue;
      const CharBlock& blk = ds.blocks[b];
      int off = ds.block_word_offset[b];
      int k = blk.n_states;

      const uint64_t* L2 = &tree.down2[lb + off];
      const uint64_t* R2 = &tree.down2[rb + off];
      const uint64_t* D2 = &tree.down2[nb + off];

      // Node's applicable second-downpass states
      uint64_t ss_app = 0;
      for (int s = 1; s < k; ++s) ss_app |= D2[s];

      // D2 intersection (including NA state)
      uint64_t any_d2_isect = 0;
      for (int s = 0; s < k; ++s) any_d2_isect |= (L2[s] & R2[s]);

      // subtree_actives: applicable tips anywhere in children's subtrees
      const uint64_t* la = &tree.subtree_actives[lb + off];
      const uint64_t* ra = &tree.subtree_actives[rb + off];
      uint64_t l_act = 0, r_act = 0;
      for (int s = 1; s < k; ++s) { l_act |= la[s]; r_act |= ra[s]; }

      // Same formula as fitch_na_score Pass 3
      uint64_t needs_step = l_act & r_act
                          & ~(ss_app & any_d2_isect) & blk.active_mask;

      while (needs_step) {
        int c = ctz64(needs_step);
        char_steps[blk.pattern_index[c]] += 1;
        needs_step &= needs_step - 1;
      }
    }
  }
}

// --- IW computation ---

double compute_iw(const DataSet& ds, const std::vector<int>& char_steps) {
  double score = 0.0;
  for (int p = 0; p < ds.n_patterns; ++p) {
    int extra = char_steps[p] - ds.min_steps[p];
    if (extra > 0) {
      score += ds.pattern_freq[p] * ds.phi[p] *
               (static_cast<double>(extra) / (ds.eff_k[p] + extra));
    }
  }
  return score;
}

void precompute_iw_delta(const DataSet& ds,
                         const std::vector<int>& divided_steps,
                         std::vector<double>& iw_delta) {
  for (int p = 0; p < ds.n_patterns; ++p) {
    int e = divided_steps[p] - ds.min_steps[p];
    if (e < 0) {
      // Below minimum in divided tree: adding one step still leaves at or
      // below minimum, so reconnection incurs no IW cost increase.
      iw_delta[p] = 0.0;
      continue;
    }
    double k = ds.eff_k[p];
    // Guard e == 0: old_cost = 0 for any k >= 0, avoiding 0/0 NaN when k == 0
    double old_cost = (e == 0) ? 0.0
                                : static_cast<double>(e) / (k + e);
    double new_cost = static_cast<double>(e + 1) / (k + e + 1);
    iw_delta[p] = ds.pattern_freq[p] * ds.phi[p] * (new_cost - old_cost);
  }
}

double indirect_iw_length(
    const uint64_t* clip_prelim,
    const TreeState& tree, const DataSet& ds,
    int node_a, int node_d,
    double base_iw,
    const std::vector<double>& iw_delta) {

  double candidate_iw = base_iw;

  size_t a_base = static_cast<size_t>(node_a) * tree.total_words;
  size_t d_base = static_cast<size_t>(node_d) * tree.total_words;

  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    if (blk.active_mask == 0) continue;
    int offset = ds.block_word_offset[b];

    uint64_t any_hit = simd::any_hit_reduce3(
        &clip_prelim[offset],
        &tree.final_[a_base + offset],
        &tree.final_[d_base + offset],
        blk.n_states);

    uint64_t needs_step = ~any_hit & blk.active_mask;
    while (needs_step) {
      int c = ctz64(needs_step);
      candidate_iw += iw_delta[blk.pattern_index[c]];
      needs_step &= needs_step - 1;
    }
  }

  return candidate_iw;
}

// Early-termination IW variant.
double indirect_iw_length_bounded(
    const uint64_t* clip_prelim,
    const TreeState& tree, const DataSet& ds,
    int node_a, int node_d,
    double base_iw,
    const std::vector<double>& iw_delta,
    double cutoff) {

  double candidate_iw = base_iw;

  size_t a_base = static_cast<size_t>(node_a) * tree.total_words;
  size_t d_base = static_cast<size_t>(node_d) * tree.total_words;

  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    if (blk.active_mask == 0) continue;
    int offset = ds.block_word_offset[b];

    uint64_t any_hit = simd::any_hit_reduce3(
        &clip_prelim[offset],
        &tree.final_[a_base + offset],
        &tree.final_[d_base + offset],
        blk.n_states);

    uint64_t needs_step = ~any_hit & blk.active_mask;
    while (needs_step) {
      int c = ctz64(needs_step);
      candidate_iw += iw_delta[blk.pattern_index[c]];
      needs_step &= needs_step - 1;
    }
    if (candidate_iw >= cutoff) return candidate_iw;
  }

  return candidate_iw;
}

// Precomputed-vroot IW variant with early termination.
double indirect_iw_length_cached(
    const uint64_t* clip_prelim,
    const uint64_t* vroot,
    const DataSet& ds,
    double base_iw,
    const std::vector<double>& iw_delta,
    double cutoff) {

  double candidate_iw = base_iw;

  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    if (blk.active_mask == 0) continue;
    int offset = ds.block_word_offset[b];

    uint64_t any_hit = simd::any_hit_reduce(
        &clip_prelim[offset], &vroot[offset], blk.n_states);

    uint64_t needs_step = ~any_hit & blk.active_mask;
    while (needs_step) {
      int c = ctz64(needs_step);
      candidate_iw += iw_delta[blk.pattern_index[c]];
      needs_step &= needs_step - 1;
    }
    if (candidate_iw >= cutoff) return candidate_iw;
  }

  return candidate_iw;
}

// --- Profile parsimony scoring ---

double compute_profile(const DataSet& ds, const std::vector<int>& char_steps) {
  double score = 0.0;
  for (int p = 0; p < ds.n_patterns; ++p) {
    // Add back precomputed_steps to get the original total step count,
    // since info_amounts is indexed by the full (unsimplified) step count.
    int s = char_steps[p];
    if (!ds.precomputed_steps.empty()) s += ds.precomputed_steps[p];
    // info_amounts is column-major: [(s-1) + info_max_steps * p]
    // s is total step count (1-based); row 0 = 1 total step
    int idx = s - 1;
    if (idx >= 0 && idx < ds.info_max_steps) {
      score += ds.pattern_freq[p] * ds.info_amounts[idx + ds.info_max_steps * p];
    }
    // s == 0: invariant character → 0 cost
    // s > info_max_steps: beyond table → treat as max cost
    else if (s > ds.info_max_steps && ds.info_max_steps > 0) {
      score += ds.pattern_freq[p] *
               ds.info_amounts[(ds.info_max_steps - 1) + ds.info_max_steps * p];
    }
  }
  return score;
}

void precompute_profile_delta(const DataSet& ds,
                               const std::vector<int>& divided_steps,
                               std::vector<double>& delta) {
  for (int p = 0; p < ds.n_patterns; ++p) {
    int s = divided_steps[p];  // reduced steps in divided tree
    // Add back precomputed_steps to get the original total step count,
    // since info_amounts is indexed by the full (unsimplified) step count.
    if (!ds.precomputed_steps.empty()) s += ds.precomputed_steps[p];
    int idx_old = s - 1;       // 0-based row for current step count
    int idx_new = s;           // 0-based row for step count + 1

    // old_cost: 0 if invariant (s<=0), capped at max if beyond table.
    // Note: must mirror the capping in compute_profile() to avoid overestimating
    // delta when divided_steps already exceeds info_max_steps (S-RED focus 10).
    double old_cost;
    if (idx_old < 0) {
      old_cost = 0.0;
    } else if (idx_old < ds.info_max_steps) {
      old_cost = ds.info_amounts[idx_old + ds.info_max_steps * p];
    } else {
      old_cost = ds.info_amounts[(ds.info_max_steps - 1) + ds.info_max_steps * p];
    }

    double new_cost;
    if (idx_new >= 0 && idx_new < ds.info_max_steps) {
      new_cost = ds.info_amounts[idx_new + ds.info_max_steps * p];
    } else if (idx_new >= ds.info_max_steps && ds.info_max_steps > 0) {
      // Beyond table: cap at maximum cost (delta = 0 from this point)
      new_cost = ds.info_amounts[(ds.info_max_steps - 1) + ds.info_max_steps * p];
    } else {
      new_cost = 0.0;
    }

    delta[p] = ds.pattern_freq[p] * (new_cost - old_cost);
  }
}

// --- Weighted scoring dispatch (IW or profile) ---

double compute_weighted_score(const DataSet& ds,
                               const std::vector<int>& char_steps) {
  if (ds.scoring_mode == ScoringMode::PROFILE)
    return compute_profile(ds, char_steps);
  return compute_iw(ds, char_steps);
}

void precompute_weighted_delta(const DataSet& ds,
                                const std::vector<int>& divided_steps,
                                std::vector<double>& delta) {
  if (ds.scoring_mode == ScoringMode::PROFILE)
    precompute_profile_delta(ds, divided_steps, delta);
  else
    precompute_iw_delta(ds, divided_steps, delta);
}

// --- Unified scoring ---

double fitch_score_ew(TreeState& tree, const DataSet& ds) {
  // Check if any block has inapplicable characters
  bool has_na = false;
  for (int b = 0; b < ds.n_blocks; ++b) {
    if (ds.blocks[b].has_inapplicable) { has_na = true; break; }
  }

  if (ds.scoring_mode == ScoringMode::EW ||
      ds.scoring_mode == ScoringMode::HSJ ||
      ds.scoring_mode == ScoringMode::XFORM) {
    // Equal weights — add back precomputed topology-independent steps
    // (HSJ/XFORM non-hierarchy chars also use EW scoring)
    if (has_na) {
      return static_cast<double>(fitch_na_score(tree, ds)) + ds.ew_offset;
    } else {
      return static_cast<double>(fitch_score(tree, ds)) + ds.ew_offset;
    }
  }

  // Weighted scoring (IW or profile): run Fitch, extract per-pattern steps,
  // transform. Reusable scratch: this fires per accept / convergence check /
  // sector+Wagner rescore AND per exact_verify candidate (IW/profile only — EW
  // returns above). thread_local keeps capacity across calls.
  static thread_local std::vector<int> char_steps;
  char_steps.resize(ds.n_patterns);

  // FUSED extraction (NA path): fitch_na_score fills char_steps per-pattern
  // during its own Pass1/Pass3, avoiding a redundant full-tree extract_char_steps
  // re-walk (~14.6% of the per-candidate NA-IW rescore; see native-na-tbr-frontier).
  // Byte-identical. TS_NA_NOFUSE (read once) restores the separate-walk path for A/B.
  static const bool na_nofuse = std::getenv("TS_NA_NOFUSE") != nullptr;
  if (has_na && !na_nofuse) {
    fitch_na_score(tree, ds, &char_steps);
    return compute_weighted_score(ds, char_steps);
  }

  if (has_na) {
    fitch_na_score(tree, ds);
  } else {
    fitch_score(tree, ds);
  }
  extract_char_steps(tree, ds, char_steps);
  return compute_weighted_score(ds, char_steps);
}

double score_tree(TreeState& tree, const DataSet& ds) {
  if (ds.scoring_mode == ScoringMode::HSJ) {
    // HSJ: Fitch on non-hierarchy chars + HSJ DP on hierarchy blocks.
    // hsj_score() calls fitch_score_ew() internally, avoiding recursion.
    return hsj_score(tree, ds);
  }
  if (ds.scoring_mode == ScoringMode::XFORM) {
    // Xform: Fitch on non-hierarchy chars + Sankoff on recoded hierarchy chars.
    double fitch_part = fitch_score_ew(tree, ds);
    if (ds.sankoff_n_chars > 0) {
      // Build SankoffData on the fly from DataSet fields.
      // Reuse the Sankoff engine's multi-char scoring.
      SankoffData sd;
      sd.n_tips = tree.n_tip;
      sd.n_chars = ds.sankoff_n_chars;
      sd.max_states = ds.sankoff_max_states;
      sd.chars.resize(sd.n_chars);
      for (int ch = 0; ch < sd.n_chars; ++ch) {
        const int ns = ds.sankoff_n_states[ch];   // per-character state count
        const int ms = ds.sankoff_max_states;     // source block row-stride
        sd.chars[ch].n_states = ns;
        sd.chars[ch].forced_root_state = ds.sankoff_forced_root[ch];
        // sankoff_score_char reads cost_matrix at the per-character stride,
        // cost_matrix[s * ns + t]; the source block in ds.sankoff_cost_matrices
        // is laid out at the max_states stride (unpack_xform writes
        // dst[r * max_states + c]). Compact the top-left ns x ns block down to
        // the ns stride here, exactly as ts_sankoff_test does. Copying the whole
        // max_states^2 block verbatim left a stride mismatch whenever
        // ns < max_states (two recoded Sankoff chars of differing state counts):
        // rows s > 0 were then read from the zero-padded gap, silently treating
        // transition (loss) costs as 0 and undercounting the tree's score.
        sd.chars[ch].cost_matrix.resize(static_cast<size_t>(ns) * ns);
        const double* src = ds.sankoff_cost_matrices.data() +
            static_cast<size_t>(ch) * ms * ms;
        for (int r = 0; r < ns; ++r) {
          for (int c = 0; c < ns; ++c) {
            sd.chars[ch].cost_matrix[static_cast<size_t>(r) * ns + c] =
                src[static_cast<size_t>(r) * ms + c];
          }
        }
      }
      sd.tip_costs = ds.sankoff_tip_costs;

      fitch_part += sankoff_score(
          tree.left.data(), tree.right.data(),
          tree.postorder.data(), tree.n_internal,
          tree.n_tip, sd);
    }
    return fitch_part;
  }
  return fitch_score_ew(tree, ds);
}


// =========================================================================
// Inapplicable (NA) three-pass scoring (Brazeau et al. 2019)
// =========================================================================

#include "ts_fitch_na.h"

// =========================================================================
// Incremental NA-aware scoring for SPR/TBR candidate evaluation
// =========================================================================

#include "ts_fitch_na_incr.h"

// =========================================================================
// NA-aware dirty-set incremental rescore (T-300 NA variant)
// =========================================================================

#include "ts_fitch_na_dirty.h"

} // namespace ts
