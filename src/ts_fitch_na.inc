// This file is #included at the end of ts_fitch.cpp
// Inapplicable (NA) three-pass scoring (Brazeau et al. 2019)

int fitch_na_score(TreeState& tree, const DataSet& ds) {
  bool has_any_na = false;
  for (int b = 0; b < ds.n_blocks; ++b) {
    if (ds.blocks[b].has_inapplicable) { has_any_na = true; break; }
  }
  if (!has_any_na) return fitch_score(tree, ds);

  int total_steps = 0;

  // ==== Pass 1: First downpass ====
  for (int node : tree.postorder) {
    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];
    size_t nb = static_cast<size_t>(node) * tree.total_words;
    size_t lb = static_cast<size_t>(lc) * tree.total_words;
    size_t rb = static_cast<size_t>(rc) * tree.total_words;

    for (int b = 0; b < ds.n_blocks; ++b) {
      const CharBlock& blk = ds.blocks[b];
      int off = ds.block_word_offset[b];
      int k = blk.n_states;
      const uint64_t* L = &tree.prelim[lb + off];
      const uint64_t* R = &tree.prelim[rb + off];
      uint64_t* N = &tree.prelim[nb + off];

      if (!blk.has_inapplicable) {
        // Standard Fitch
        uint64_t any_isect = 0;
        for (int s = 0; s < k; ++s) any_isect |= (L[s] & R[s]);
        uint64_t needs_union = ~any_isect & blk.active_mask;
        total_steps += blk.weight * popcount64(needs_union);
        tree.local_cost[static_cast<size_t>(node) * tree.n_blocks + b] = needs_union;
        for (int s = 0; s < k; ++s) {
          uint64_t isect = L[s] & R[s];
          N[s] = (isect & any_isect) | ((L[s] | R[s]) & needs_union);
        }
      } else {
        // NA-aware first downpass (Brazeau et al. 2019)
        // Matches morphy's mpl_NA_fitch_first_downpass exactly.
        uint64_t I_app = 0, L_app = 0, R_app = 0;
        for (int s = 1; s < k; ++s) {
          I_app |= (L[s] & R[s]);
          L_app |= L[s];
          R_app |= R[s];
        }
        uint64_t I0 = L[0] & R[0];
        uint64_t both_app = L_app & R_app;
        // case_keep: use intersection (applicable states intersect,
        //   OR only NA intersects and not both children applicable)
        uint64_t case_keep = I_app | (I0 & ~I_app & ~both_app);
        // case_strip: strip NA from union (no intersection at all,
        //   both children have applicable states)
        uint64_t case_strip = ~I0 & ~I_app & both_app;
        for (int s = 1; s < k; ++s) {
          uint64_t isect = L[s] & R[s];
          N[s] = (isect & case_keep) | ((L[s] | R[s]) & ~case_keep);
        }
        N[0] = (I0 & case_keep) | ((L[0] | R[0]) & ~case_keep & ~case_strip);

        // Subtree actives
        const uint64_t* la = &tree.subtree_actives[lb + off];
        const uint64_t* ra = &tree.subtree_actives[rb + off];
        uint64_t* na = &tree.subtree_actives[nb + off];
        na[0] = 0;
        for (int s = 1; s < k; ++s) na[s] = la[s] | ra[s];
      }
    }
  }

  // ==== Pass 2: First uppass (applicability propagation) ====
  int root = tree.n_tip;
  size_t root_b = static_cast<size_t>(root) * tree.total_words;
  for (int w = 0; w < tree.total_words; ++w) {
    tree.final_[root_b + w] = tree.prelim[root_b + w];
  }

  // Internal nodes in reverse postorder
  for (int i = static_cast<int>(tree.postorder.size()) - 1; i >= 0; --i) {
    int node = tree.postorder[i];
    if (node == root) continue;
    int anc = tree.parent[node];
    size_t nb = static_cast<size_t>(node) * tree.total_words;
    size_t ab = static_cast<size_t>(anc) * tree.total_words;

    for (int b = 0; b < ds.n_blocks; ++b) {
      const CharBlock& blk = ds.blocks[b];
      int off = ds.block_word_offset[b];
      int k = blk.n_states;

      if (!blk.has_inapplicable) {
        // Standard uppass
        uint64_t any_isect = 0;
        for (int s = 0; s < k; ++s) {
          any_isect |= (tree.final_[ab + off + s] & tree.prelim[nb + off + s]);
        }
        uint64_t no_isect = ~any_isect & blk.active_mask;
        for (int s = 0; s < k; ++s) {
          uint64_t isect = tree.final_[ab + off + s] & tree.prelim[nb + off + s];
          tree.final_[nb + off + s] = (isect & any_isect)
                                    | (tree.prelim[nb + off + s] & no_isect);
        }
      } else {
        // NA-aware first uppass
        const uint64_t* Np = &tree.prelim[nb + off];
        const uint64_t* A = &tree.final_[ab + off];
        uint64_t* F = &tree.final_[nb + off];

        uint64_t npre_has_NA = Np[0];
        uint64_t npre_has_app = 0;
        for (int s = 1; s < k; ++s) npre_has_app |= Np[s];
        uint64_t anc_app = 0;
        for (int s = 1; s < k; ++s) anc_app |= A[s];
        uint64_t anc_is_NA = A[0] & ~anc_app;

        uint64_t children_app = 0;
        if (node >= tree.n_tip) {
          int ni2 = node - tree.n_tip;
          size_t cl = static_cast<size_t>(tree.left[ni2]) * tree.total_words + off;
          size_t cr = static_cast<size_t>(tree.right[ni2]) * tree.total_words + off;
          for (int s = 1; s < k; ++s) {
            children_app |= (tree.prelim[cl + s] | tree.prelim[cr + s]);
          }
        }

        uint64_t case_pass = ~npre_has_NA & blk.active_mask;
        uint64_t case_strip = npre_has_NA & npre_has_app & ~anc_is_NA;
        uint64_t case_children = npre_has_NA & ~npre_has_app & ~anc_is_NA & children_app;
        uint64_t case_force = blk.active_mask & ~case_pass & ~case_strip & ~case_children;

        F[0] = (Np[0] & case_pass) | case_force;
        for (int s = 1; s < k; ++s) {
          uint64_t child_union = 0;
          if (node >= tree.n_tip) {
            int ni2 = node - tree.n_tip;
            child_union = tree.prelim[static_cast<size_t>(tree.left[ni2]) * tree.total_words + off + s]
                        | tree.prelim[static_cast<size_t>(tree.right[ni2]) * tree.total_words + off + s];
          }
          F[s] = (Np[s] & (case_pass | case_strip)) | (child_union & case_children);
        }

        // For tips in postorder: also set down2 = final
        if (node < tree.n_tip) {
          uint64_t* D2 = &tree.down2[nb + off];
          for (int s = 0; s < k; ++s) D2[s] = F[s];
        }
      }
    }
  }

  // Process tips for uppass (they are NOT in postorder)
  for (int tip = 0; tip < tree.n_tip; ++tip) {
    int anc = tree.parent[tip];
    size_t tb = static_cast<size_t>(tip) * tree.total_words;
    size_t ab = static_cast<size_t>(anc) * tree.total_words;

    for (int b = 0; b < ds.n_blocks; ++b) {
      const CharBlock& blk = ds.blocks[b];
      int off = ds.block_word_offset[b];
      int k = blk.n_states;

      if (!blk.has_inapplicable) {
        uint64_t any_isect = 0;
        for (int s = 0; s < k; ++s) {
          any_isect |= (tree.final_[ab + off + s] & tree.prelim[tb + off + s]);
        }
        uint64_t no_isect = ~any_isect & blk.active_mask;
        for (int s = 0; s < k; ++s) {
          uint64_t isect = tree.final_[ab + off + s] & tree.prelim[tb + off + s];
          tree.final_[tb + off + s] = (isect & any_isect)
                                    | (tree.prelim[tb + off + s] & no_isect);
        }
      } else {
        // NA-aware tip update (matches morphy's mpl_fitch_NA_tip_update)
        const uint64_t* T = &tree.prelim[tb + off];
        const uint64_t* A = &tree.final_[ab + off];
        uint64_t* F = &tree.final_[tb + off];

        // Check if tip intersects ancestor (any state word)
        uint64_t any_isect = 0;
        for (int s = 0; s < k; ++s) any_isect |= (T[s] & A[s]);

        // Ancestor has applicable states?
        uint64_t anc_app = 0;
        for (int s = 1; s < k; ++s) anc_app |= A[s];

        // Strip NA from tip only when intersection exists AND
        // ancestor is applicable
        uint64_t strip_na = any_isect & anc_app;

        F[0] = T[0] & ~strip_na;
        for (int s = 1; s < k; ++s) F[s] = T[s];

        // down2 = final for tips
        uint64_t* D2 = &tree.down2[tb + off];
        for (int s = 0; s < k; ++s) D2[s] = F[s];

        // Update tip subtree_actives (morphy's tip update logic):
        // intersection exists → applicable part of intersection
        // no intersection → tip's applicable states
        uint64_t* sa = &tree.subtree_actives[tb + off];
        for (int s = 1; s < k; ++s) {
          sa[s] = (T[s] & A[s] & any_isect) | (T[s] & ~any_isect);
        }
        sa[0] = 0;
      }
    }
  }

  // ==== Pass 3: Second downpass (corrected scoring) ====
  // Matches morphy's mpl_NA_fitch_second_downpass exactly.
  //
  // For each internal node, three step-counting conditions:
  //   (a) Applicable node, no D2 intersection at all, both children
  //       have applicable D2 states (standard Fitch step).
  //   (b) Applicable node, no D2 intersection, but both SUBTREES
  //       have applicable tips (region-separation step).
  //   (c) Inapplicable node, but both subtrees have applicable tips
  //       (region-separation step through inapplicable node).
  //
  // Simplified: step when l_act & r_act & ~(ss_app & any_d2_isect).
  for (int node : tree.postorder) {
    int ni = node - tree.n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];
    size_t nb = static_cast<size_t>(node) * tree.total_words;
    size_t lb = static_cast<size_t>(lc) * tree.total_words;
    size_t rb = static_cast<size_t>(rc) * tree.total_words;

    for (int b = 0; b < ds.n_blocks; ++b) {
      if (!ds.blocks[b].has_inapplicable) continue;
      const CharBlock& blk = ds.blocks[b];
      int off = ds.block_word_offset[b];
      int k = blk.n_states;

      const uint64_t* F = &tree.final_[nb + off];
      const uint64_t* L2 = &tree.down2[lb + off];
      const uint64_t* R2 = &tree.down2[rb + off];
      uint64_t* D2 = &tree.down2[nb + off];

      // Node applicability from first uppass
      uint64_t ss_app = 0;
      for (int s = 1; s < k; ++s) ss_app |= F[s];

      // Any intersection in children's D2 (including NA bit)
      uint64_t any_isect = 0;
      for (int s = 0; s < k; ++s) any_isect |= (L2[s] & R2[s]);

      // Applicable intersection
      uint64_t I_app = 0;
      for (int s = 1; s < k; ++s) I_app |= (L2[s] & R2[s]);

      // Children's subtree_actives: applicable tips anywhere below
      const uint64_t* la = &tree.subtree_actives[lb + off];
      const uint64_t* ra = &tree.subtree_actives[rb + off];
      uint64_t l_act = 0, r_act = 0;
      for (int s = 1; s < k; ++s) { l_act |= la[s]; r_act |= ra[s]; }

      // Step: both subtrees have applicable tips, AND NOT (node is
      //       applicable with some D2 intersection).
      uint64_t needs_step = l_act & r_act
                          & ~(ss_app & any_isect) & blk.active_mask;
      total_steps += blk.weight * popcount64(needs_step);

      // Compute down2 for this node (matches morphy exactly):
      // Applicable node with any intersection:
      //   - applicable isect → D2 = applicable part of intersection
      //   - NA-only isect → D2 = {NA} (node "loses" applicable status)
      // Applicable node with no intersection:
      //   - D2 = applicable-only union
      // Inapplicable node: D2 = first uppass state (= {NA})
      uint64_t na_only_isect = any_isect & ~I_app;
      for (int s = 1; s < k; ++s) {
        uint64_t isect = L2[s] & R2[s];
        uint64_t uni = L2[s] | R2[s];
        D2[s] = ss_app & ((isect & any_isect) | (uni & ~any_isect));
      }
      D2[0] = (~ss_app | na_only_isect) & blk.active_mask;

      // Update subtree_actives: union of children's applicable actives
      uint64_t* na_out = &tree.subtree_actives[nb + off];
      na_out[0] = 0;
      for (int s = 1; s < k; ++s) na_out[s] = (la[s] | ra[s]);
    }
  }

  return total_steps;
}
