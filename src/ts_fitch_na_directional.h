// This file is #included at the end of ts_fitch.cpp, after ts_fitch_na_incr.h.
//
// EXACT directional (Regime-C) three-pass scoring for NA TBR candidates: an
// O(1)/character per-candidate combine that replaces the O(n) full_rescore in the
// NA TBR path.  Validated bit-for-bit in R (dev/prototype/rt_tbr.R) and in a
// standalone C++ port (dev/prototype/na_directional.cpp, ~51k candidates, 0 mis).
//
// A per-CHARACTER message carries the m2-indexed 2^k response table
// m2 -> (steps, down2, sact).  This table is load-bearing for exactness: a
// fragment's internal steps are a FUNCTION of the boundary message pushed from
// above (a single-state summary is the APPROXIMATE fitch_na_indirect_length).
// Layout is per-character flat tables (NOT block-bitwise) -- the perf-proven
// architecture from bench_directional.cpp.
//
// M2 status: the AUDIT path (na_dir_audit_score) is wired into exact_verify_sweep
// behind TS_EV_AUDIT to assert directional == full_rescore on every candidate,
// under base AND perturbed (ratchet) weighting.  It is NOT yet the scorer.

#ifndef TS_NA_DIRECTIONAL_H
#define TS_NA_DIRECTIONAL_H

namespace nadir {

// ---- node rules (verbatim port of the M1-validated na_directional.cpp) -------
static inline int down1_rule(int Lm, int Rm, int appmask) {
  int LR = Lm & Rm;
  int I_app = LR & appmask, L_app = Lm & appmask, R_app = Rm & appmask;
  int I0 = LR & 1;
  bool both_app = (L_app != 0) && (R_app != 0);
  bool ia = (I_app != 0), I0b = (I0 != 0);
  bool case_keep = ia || (I0b && !ia && !both_app);
  bool case_strip = (!I0b) && (!ia) && both_app;
  int N_app = case_keep ? (LR & appmask) : ((Lm | Rm) & appmask);
  int NAu = (Lm | Rm) & 1;
  int N0 = case_keep ? I0 : (!case_strip ? NAu : 0);
  return N_app | N0;
}
static inline int up1_rule(int A, int Np, int child_union, int appmask) {
  bool npre_has_NA = (Np & 1) != 0, npre_has_app = (Np & appmask) != 0;
  bool anc_app = (A & appmask) != 0, anc_is_NA = ((A & 1) != 0) && !anc_app;
  bool children_app = (child_union & appmask) != 0;
  bool case_pass = !npre_has_NA;
  bool case_strip = npre_has_NA && npre_has_app && !anc_is_NA;
  bool case_children = npre_has_NA && !npre_has_app && !anc_is_NA && children_app;
  bool case_force = !case_pass && !case_strip && !case_children;
  int F0 = case_pass ? (Np & 1) : 0;
  if (case_force) F0 = 1;
  int Fapp = (case_pass || case_strip) ? (Np & appmask)
           : (case_children ? (child_union & appmask) : 0);
  return F0 | Fapp;
}
struct D2 { int step, down2, sact; };
static inline D2 down2_rule(int Fm, int L2, int R2, int sL, int sR, int appmask) {
  bool ss_app = (Fm & appmask) != 0;
  int isect = L2 & R2;
  bool any_isect = (isect != 0), I_app = (isect & appmask) != 0;
  bool l_act = (sL & appmask) != 0, r_act = (sR & appmask) != 0;
  int step = (l_act && r_act && !(ss_app && any_isect)) ? 1 : 0;
  bool na_only = any_isect && !I_app;
  int D2app = ss_app ? (any_isect ? (isect & appmask) : ((L2 | R2) & appmask)) : 0;
  int down2 = D2app | ((!ss_app || na_only) ? 1 : 0);
  int sact = (sL | sR) & appmask;
  return { step, down2, sact };
}

struct Msg {
  bool is_tip;
  int prelim_root, child_union_root;
  std::vector<int> steps, down2, sact;   // 2^k entries; empty for tips
};
struct PR { int steps, down2, sact; };

static inline PR push_to_child(int parent_final, const Msg& m, int appmask) {
  if (m.is_tip) {
    int Tm = m.prelim_root, A = parent_final;
    bool any_isect = (Tm & A) != 0, anc_app = (A & appmask) != 0;
    bool strip_na = any_isect && anc_app;
    int final_t = (strip_na ? 0 : (Tm & 1)) | (Tm & appmask);
    int sact_t = any_isect ? (Tm & A & appmask) : (Tm & appmask);
    return { 0, final_t, sact_t };
  }
  int m2c = up1_rule(parent_final, m.prelim_root, m.child_union_root, appmask);
  return { m.steps[m2c], m.down2[m2c], m.sact[m2c] };
}
static inline Msg tip_msg(int Tm) {
  Msg m; m.is_tip = true; m.prelim_root = Tm; m.child_union_root = 0; return m;
}
static Msg combine_msgs(const Msg& mL, const Msg& mR, int k, int appmask) {
  int full = (1 << k) - 1;
  Msg out; out.is_tip = false;
  out.prelim_root = down1_rule(mL.prelim_root, mR.prelim_root, appmask);
  out.child_union_root = mL.prelim_root | mR.prelim_root;
  out.steps.resize(full + 1); out.down2.resize(full + 1); out.sact.resize(full + 1);
  for (int m2 = 0; m2 <= full; ++m2) {
    PR rL = push_to_child(m2, mL, appmask), rR = push_to_child(m2, mR, appmask);
    D2 d = down2_rule(m2, rL.down2, rR.down2, rL.sact, rR.sact, appmask);
    out.steps[m2] = rL.steps + rR.steps + d.step;
    out.down2[m2] = d.down2; out.sact[m2] = d.sact;
  }
  return out;
}
static int combine_at_root_msgs(const Msg& mA, const Msg& mB, int appmask) {
  int prelim_r = down1_rule(mA.prelim_root, mB.prelim_root, appmask);
  PR rA = push_to_child(prelim_r, mA, appmask), rB = push_to_child(prelim_r, mB, appmask);
  D2 d = down2_rule(prelim_r, rA.down2, rB.down2, rA.sact, rB.sact, appmask);
  return rA.steps + rB.steps + d.step;
}

// ---- a "piece" of the tree (F or H), in shared original node ids ------------
struct Piece {
  int root;
  std::vector<int> post;                 // internal nodes, children before parents
  std::vector<int> left, right, parent;  // size n_node; valid only for piece nodes
  std::vector<char> in_piece;            // membership, size n_node
};

// Per-character extracted tip mask for character c of block b at tip `tip`, plus
// the combine alphabet size.  NA blocks: word w -> bit w (bit0 = inapplicable).
// Standard blocks: applicable state words shifted into bits 1.. (bit0 stays 0, so
// the three-pass reduces to Fitch).
static inline void char_alphabet(const CharBlock& blk, int& k_combine, int& appmask) {
  k_combine = blk.has_inapplicable ? blk.n_states : blk.n_states + 1;
  appmask = (1 << k_combine) - 1 - 1;
}
// Reads the CANONICAL tip data (ds.tip_states, already POST genuine-inapp strip
// from ts_data.cpp), so it is valid regardless of transient tree.prelim state.
static inline int extract_mask(const DataSet& ds, int b, int off, int c, int tip) {
  const CharBlock& blk = ds.blocks[b];
  int k = blk.n_states;
  size_t tb = static_cast<size_t>(tip) * ds.total_words + off;
  int mask = 0;
  uint64_t bit = 1ULL << c;
  if (blk.has_inapplicable) {
    for (int w = 0; w < k; ++w)
      if (ds.tip_states[tb + w] & bit) mask |= (1 << w);
  } else {
    for (int w = 0; w < k; ++w)
      if (ds.tip_states[tb + w] & bit) mask |= (1 << (w + 1));
  }
  return mask;
}

// fold_down + fold_up for a piece, for ONE character (tip masks supplied by id).
static void fold_piece(const Piece& pc, const std::vector<int>& tipmask_byid,
                       int n_tip, int k, int appmask,
                       std::vector<Msg>& Rd, std::vector<Msg>& Ru) {
  int n_node = static_cast<int>(pc.in_piece.size());
  Rd.assign(n_node, Msg());
  Ru.assign(n_node, Msg());
  for (int v = 0; v < n_node; ++v)
    if (pc.in_piece[v] && v < n_tip) Rd[v] = tip_msg(tipmask_byid[v]);
  for (int v : pc.post) Rd[v] = combine_msgs(Rd[pc.left[v]], Rd[pc.right[v]], k, appmask);
  // fold_up: preorder (parents first) = reverse postorder of internals
  for (auto it = pc.post.rbegin(); it != pc.post.rend(); ++it) {
    int p = *it; int a = pc.left[p], b = pc.right[p];
    if (p == pc.root) { Ru[a] = Rd[b]; Ru[b] = Rd[a]; }
    else { Ru[a] = combine_msgs(Rd[b], Ru[p], k, appmask);
           Ru[b] = combine_msgs(Rd[a], Ru[p], k, appmask); }
  }
}

// Build F (subtree at clip_node) and H (divided host) as pieces, on the CLEAN tree.
// Assumes clip_node's parent nx is not the display root (caller's guard).
// NB the kernel stores tree.left/right by INTERNAL index (node - n_tip) but
// tree.parent by NODE id; the Piece arrays are all NODE-indexed for uniformity.
static void build_pieces(const TreeState& tree, int clip_node, Piece& F, Piece& H) {
  int n_node = tree.n_node, n_tip = tree.n_tip;
  int nx = tree.parent[clip_node];
  int nz = tree.parent[nx];
  int ns = (tree.left[nx - n_tip] == clip_node) ? tree.right[nx - n_tip]
                                                : tree.left[nx - n_tip];

  auto node_indexed = [&](Piece& P) {
    P.left.assign(n_node, -1); P.right.assign(n_node, -1);
    P.parent = tree.parent;                         // already node-indexed
    for (int v = n_tip; v < n_node; ++v) {
      P.left[v]  = tree.left[v - n_tip];
      P.right[v] = tree.right[v - n_tip];
    }
    P.in_piece.assign(n_node, 0);
  };
  auto dfs_post = [&](Piece& P) {                   // postorder of internal nodes
    P.post.clear();
    std::vector<std::pair<int,int>> st; st.push_back({P.root, 1});
    while (!st.empty()) {
      auto pr = st.back(); st.pop_back();
      int v = pr.first, stage = pr.second;
      P.in_piece[v] = 1;
      if (v < n_tip) continue;
      if (stage == 1) { st.push_back({v, 2});
        st.push_back({P.left[v], 1}); st.push_back({P.right[v], 1}); }
      else P.post.push_back(v);
    }
  };

  // ---- F: subtree rooted at clip_node ----
  node_indexed(F); F.root = clip_node; dfs_post(F);

  // ---- H: whole tree minus F, with ns spliced onto nz, nx removed ----
  node_indexed(H); H.root = n_tip;
  H.parent[ns] = nz;                                // splice: ns under nz
  if (H.left[nz] == nx) H.left[nz] = ns; else H.right[nz] = ns;   // nz>=n_tip always
  dfs_post(H);
}

// Per-clip fold cache for all blocks/characters.
struct ClipFolds {
  int clip_node;
  Piece F, H;
  // index: [global_char][node] -> Msg.  global_char enumerates (b,c) in order.
  std::vector<std::vector<Msg>> RdF, RuF, RdH, RuH;
  std::vector<int> ch_block, ch_index, ch_k, ch_appmask;   // per global_char
};

static void build_clip_folds(const TreeState& tree, const DataSet& ds,
                             int clip_node, ClipFolds& cf) {
  cf.clip_node = clip_node;
  build_pieces(tree, clip_node, cf.F, cf.H);
  cf.RdF.clear(); cf.RuF.clear(); cf.RdH.clear(); cf.RuH.clear();
  cf.ch_block.clear(); cf.ch_index.clear(); cf.ch_k.clear(); cf.ch_appmask.clear();
  std::vector<int> tmF(tree.n_node, 0), tmH(tree.n_node, 0);
  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    if (blk.active_mask == 0) continue;
    int off = ds.block_word_offset[b];
    int k, appmask; char_alphabet(blk, k, appmask);
    for (int c = 0; c < blk.n_chars; ++c) {
      if (!((blk.active_mask >> c) & 1ULL)) continue;   // inactive char: contributes 0
      // tip masks by original id
      for (int tip = 0; tip < tree.n_tip; ++tip) {
        int m = extract_mask(ds, b, off, c, tip);
        if (cf.F.in_piece[tip]) tmF[tip] = m;
        else tmH[tip] = m;
      }
      std::vector<Msg> RdFc, RuFc, RdHc, RuHc;
      fold_piece(cf.F, tmF, tree.n_tip, k, appmask, RdFc, RuFc);
      fold_piece(cf.H, tmH, tree.n_tip, k, appmask, RdHc, RuHc);
      cf.RdF.push_back(std::move(RdFc)); cf.RuF.push_back(std::move(RuFc));
      cf.RdH.push_back(std::move(RdHc)); cf.RuH.push_back(std::move(RuHc));
      cf.ch_block.push_back(b); cf.ch_index.push_back(c);
      cf.ch_k.push_back(k); cf.ch_appmask.push_back(appmask);
    }
  }
}

// Regime-correct directional score of a single TBR candidate.  Emits per-pattern
// step counts (the regime-independent core), then aggregates exactly as
// score_tree: EW -> weight*(1+upweight) sum + ew_offset; IW/profile ->
// compute_weighted_score (the min_steps subtraction + concave transform).
//   rp,rc = reroot edge (-1,-1 = SPR); below = host regraft edge's lower node.
static double candidate_score(const DataSet& ds,
                              const ClipFolds& cf, int rp, int rc, int below) {
  static thread_local std::vector<int> char_steps;
  char_steps.assign(ds.n_patterns, 0);
  double ew_sum = 0.0;
  int nchar = static_cast<int>(cf.ch_block.size());
  for (int g = 0; g < nchar; ++g) {
    int b = cf.ch_block[g], c = cf.ch_index[g];
    int k = cf.ch_k[g], appmask = cf.ch_appmask[g];
    const CharBlock& blk = ds.blocks[b];
    const std::vector<Msg>& RdF = cf.RdF[g];
    const std::vector<Msg>& RuF = cf.RuF[g];
    const std::vector<Msg>& RdH = cf.RdH[g];
    const std::vector<Msg>& RuH = cf.RuH[g];

    // ON-DEMAND single-entry combine (NO Msg allocation, NO 2^k loop): compute
    // only the table entries the score actually reads.  Validated identical to the
    // full-table path (M1 mode 1).  This is the path the microbench timed.
    bool spr = (rp < 0);
    const Msg* fA; const Msg* fB; int frag_prelim, frag_cu;
    if (spr) {                                           // SPR: F at its root
      fA = &RdF[cf.clip_node]; fB = nullptr;
      frag_prelim = fA->prelim_root; frag_cu = fA->child_union_root;
    } else {                                             // TBR: F rooted on edge (rp,rc)
      fA = &RdF[rc]; fB = &RuF[rc];
      frag_prelim = down1_rule(fA->prelim_root, fB->prelim_root, appmask);
      frag_cu = fA->prelim_root | fB->prelim_root;
    }
    const Msg& Hd = RdH[below];
    const Msg& Hu = RuH[below];
    int inner_prelim = down1_rule(frag_prelim, Hd.prelim_root, appmask);
    int inner_cu = frag_prelim | Hd.prelim_root;
    int prelim_r = down1_rule(inner_prelim, Hu.prelim_root, appmask);
    int mi = up1_rule(prelim_r, inner_prelim, inner_cu, appmask);
    PR fe;                                               // fragment response at mi
    if (spr) {
      fe = push_to_child(mi, *fA, appmask);
    } else {
      int mf = up1_rule(mi, frag_prelim, frag_cu, appmask);
      PR pA = push_to_child(mf, *fA, appmask), pB = push_to_child(mf, *fB, appmask);
      D2 df = down2_rule(mf, pA.down2, pB.down2, pA.sact, pB.sact, appmask);
      fe = { pA.steps + pB.steps + df.step, df.down2, df.sact };
    }
    PR he = push_to_child(mi, Hd, appmask);
    D2 di = down2_rule(mi, fe.down2, he.down2, fe.sact, he.sact, appmask);
    PR inner_r = { fe.steps + he.steps + di.step, di.down2, di.sact };
    PR rB = push_to_child(prelim_r, Hu, appmask);
    D2 dr = down2_rule(prelim_r, inner_r.down2, rB.down2, inner_r.sact, rB.sact, appmask);
    int steps = inner_r.steps + rB.steps + dr.step;

    char_steps[blk.pattern_index[c]] += steps;          // per-pattern (for IW/profile)
    bool up = blk.upweight_mask && ((blk.upweight_mask >> c) & 1ULL);
    ew_sum += static_cast<double>(blk.weight) * steps * (up ? 2 : 1);
  }
  if (ds.scoring_mode == ScoringMode::EW)
    return ew_sum + ds.ew_offset;                       // ts_fitch.cpp:1063
  return compute_weighted_score(ds, char_steps);        // IW / profile (1078)
}

// Whole-tree directional score (diagnostic): fold the un-clipped tree and
// combine at the root edge.  Must equal full_rescore in EW regime.
static double whole_score(const TreeState& tree, const DataSet& ds) {
  int n_node = tree.n_node, n_tip = tree.n_tip;
  Piece W;
  W.left.assign(n_node, -1); W.right.assign(n_node, -1); W.parent = tree.parent;
  for (int v = n_tip; v < n_node; ++v) { W.left[v] = tree.left[v - n_tip]; W.right[v] = tree.right[v - n_tip]; }
  W.in_piece.assign(n_node, 0); W.root = n_tip;
  { std::vector<std::pair<int,int>> st; st.push_back({W.root, 1});
    while (!st.empty()) { auto pr = st.back(); st.pop_back(); int v = pr.first, stg = pr.second;
      W.in_piece[v] = 1; if (v < n_tip) continue;
      if (stg == 1) { st.push_back({v, 2}); st.push_back({W.left[v], 1}); st.push_back({W.right[v], 1}); }
      else W.post.push_back(v); } }
  int rc = W.left[n_tip];
  std::vector<int> char_steps(ds.n_patterns, 0);
  double ew_sum = 0.0;
  std::vector<int> tm(n_node, 0);
  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    if (blk.active_mask == 0) continue;
    int off = ds.block_word_offset[b];
    int k, appmask; char_alphabet(blk, k, appmask);
    for (int c = 0; c < blk.n_chars; ++c) {
      if (!((blk.active_mask >> c) & 1ULL)) continue;
      for (int tip = 0; tip < n_tip; ++tip) tm[tip] = extract_mask(ds, b, off, c, tip);
      std::vector<Msg> Rd, Ru; fold_piece(W, tm, n_tip, k, appmask, Rd, Ru);
      int steps = combine_at_root_msgs(Rd[rc], Ru[rc], appmask);
      char_steps[blk.pattern_index[c]] += steps;
      bool up = blk.upweight_mask && ((blk.upweight_mask >> c) & 1ULL);
      ew_sum += static_cast<double>(blk.weight) * steps * (up ? 2 : 1);
    }
  }
  if (ds.scoring_mode == ScoringMode::EW)
    return ew_sum + ds.ew_offset;
  return compute_weighted_score(ds, char_steps);
}

}  // namespace nadir

#endif  // TS_NA_DIRECTIONAL_H
