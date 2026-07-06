#include "ts_simplify.h"
#include <algorithm>

namespace ts {

// Count how many tips have state s as their ONLY state (unambiguous).
// Also count how many tips have state s in ANY token (including ambiguous).
static void count_state_occurrences(
    const std::vector<uint32_t>& tips, int n_tips, int n_states,
    int inapp_state,
    std::vector<int>& unambig_count,    // [n_states]
    std::vector<int>& any_count,        // [n_states]
    std::vector<int>& unambig_tip_idx)  // [n_states] last tip with unambig
{
  std::fill(unambig_count.begin(), unambig_count.end(), 0);
  std::fill(any_count.begin(), any_count.end(), 0);
  std::fill(unambig_tip_idx.begin(), unambig_tip_idx.end(), -1);

  // Dataset-wide count of applicable (non-inapplicable) states. A tip whose
  // token contains EVERY applicable state is full-"?" missing data: in Fitch it
  // intersects any sibling set, so it never forces a step and never constrains
  // the reconstruction (provably inert). Such tips are ignored when counting
  // GENUINE occurrences (any_count), so missing data does not defeat the
  // singleton / redundant-state guards below (a state in one genuine tip plus
  // any number of "?" tips is still a lone-tip option, soundly removable).
  int n_applicable = 0;
  for (int s = 0; s < n_states; ++s) if (s != inapp_state) ++n_applicable;

  for (int tip = 0; tip < n_tips; ++tip) {
    uint32_t tok = tips[tip];
    // Count applicable states in this token
    int n_set = 0;
    int single_state = -1;
    for (int s = 0; s < n_states; ++s) {
      if (s == inapp_state) continue;
      if (tok & (1u << s)) {
        ++n_set;
        single_state = s;
      }
    }
    if (n_set == 1 && single_state >= 0) {
      unambig_count[single_state]++;
      unambig_tip_idx[single_state] = tip;
    }
    // Skip inert full-"?" tips for any_count (n_set >= 2 guards the single-state
    // case when n_applicable == 1, which is unambiguous, not "?").
    if (n_set >= 2 && n_set == n_applicable) continue;
    for (int s = 0; s < n_states; ++s) {
      if (s == inapp_state) continue;
      if (tok & (1u << s)) any_count[s]++;
    }
  }
}

// Check if a pattern is parsimony-uninformative: at most one state appears
// in 2+ tips (unambiguously). Only valid when all tips are unambiguous;
// caller must use verify_uninformative() when ambiguous tokens are present.
static bool is_uninformative_classical(const std::vector<int>& unambig_count,
                                        int n_states, int inapp_state) {
  int informative_states = 0;
  for (int s = 0; s < n_states; ++s) {
    if (s == inapp_state) continue;
    if (unambig_count[s] >= 2) ++informative_states;
  }
  return informative_states <= 1;
}

// Check if any tip has GENUINE ambiguity (more than one applicable state, but
// not ALL of them). A full-"?" tip (every applicable state) is inert missing
// data — it never affects the Fitch score — so it does not count as ambiguous;
// this lets the unambiguous fast path and the classical criterion see through
// missing data.
static bool has_ambiguous_tips(const std::vector<uint32_t>& tips,
                                int n_tips, int n_states, int inapp_state) {
  int n_applicable = 0;
  for (int s = 0; s < n_states; ++s) if (s != inapp_state) ++n_applicable;
  for (int tip = 0; tip < n_tips; ++tip) {
    int n_set = 0;
    for (int s = 0; s < n_states; ++s) {
      if (s == inapp_state) continue;
      if (tips[tip] & (1u << s)) ++n_set;
    }
    if (n_set > 1 && n_set < n_applicable) return true;
  }
  return false;
}

// Compute fixed step count for an uninformative pattern.
// For all-unambiguous characters: distinct states with count >= 1, minus 1.
// When ambiguous tokens are present, the caller should use
// verify_uninformative() which computes the exact fixed cost.
static int compute_fixed_steps(const std::vector<int>& unambig_count,
                                int n_states, int inapp_state) {
  int distinct = 0;
  for (int s = 0; s < n_states; ++s) {
    if (s == inapp_state) continue;
    if (unambig_count[s] >= 1) ++distinct;
  }
  return (distinct > 0) ? distinct - 1 : 0;
}

// Get the "all applicable states" mask (excluding inapp_state).
static uint32_t all_applicable_mask(const std::vector<uint32_t>& tips,
                                      int n_tips, int n_states,
                                      int inapp_state) {
  uint32_t all = 0;
  for (int tip = 0; tip < n_tips; ++tip) {
    all |= tips[tip];
  }
  if (inapp_state >= 0) all &= ~(1u << inapp_state);
  return all;
}

SimplificationResult simplify_patterns(
    const std::vector<uint32_t>& token_states,
    const int* tip_data_r, int n_tips, int n_patterns,
    const int* weight_r,
    int n_states, int inapp_state)
{
  SimplificationResult result;
  result.n_patterns_removed = 0;
  result.n_states_reduced = 0;
  result.total_offset_steps = 0;
  result.patterns.resize(n_patterns);

  std::vector<int> unambig_count(n_states);
  std::vector<int> any_count(n_states);
  std::vector<int> unambig_tip_idx(n_states);

  for (int p = 0; p < n_patterns; ++p) {
    auto& sp = result.patterns[p];
    sp.original_index = p;
    sp.precomputed_steps = 0;
    sp.informative = true;
    sp.has_genuine_inapp = false;

    // Skip zero-weight patterns (they'll be removed by build_dataset anyway)
    if (weight_r[p] == 0) {
      sp.informative = false;
      sp.n_states_remaining = 0;
      sp.tip_tokens.resize(n_tips, 0);
      result.n_patterns_removed++;
      continue;
    }

    // Build per-tip token bitmasks from the original data
    sp.tip_tokens.resize(n_tips);
    bool has_genuine_inapp = false;
    // n_states may equal MAX_STATES (32); `1u << 32` is undefined behaviour
    // (UBSAN), so build the all-ones mask directly in that case.
    uint32_t all_states_mask =
        n_states >= 32 ? ~0u : (1u << n_states) - 1;
    for (int tip = 0; tip < n_tips; ++tip) {
      int token = tip_data_r[tip + n_tips * p] - 1;  // 1-based to 0-based
      sp.tip_tokens[tip] = token_states[token];
      // Genuine inapplicable: token has inapp bit but is NOT full-?
      if (inapp_state >= 0 &&
          (token_states[token] & (1u << inapp_state)) &&
          token_states[token] != all_states_mask) {
        has_genuine_inapp = true;
      }
    }

    sp.has_genuine_inapp = has_genuine_inapp;

    // Phase 1: characters with genuine inapplicable ("-") tokens skip
    // transforms. Transforms 2/3 are not score-preserving for the NA
    // three-pass algorithm because they modify applicable state bits in
    // tokens that also carry the inapp bit.
    //
    // Characters where the inapp bit only appears in "?" (full missing
    // data) are scored with standard Fitch, so transforms ARE safe.
    // These fall through to the transform pipeline below.
    if (has_genuine_inapp) {
      // Count states for metadata only (transforms skipped)
      uint32_t all = all_applicable_mask(sp.tip_tokens, n_tips, n_states,
                                          inapp_state);
      int nc = 0;
      for (int s = 0; s < n_states; ++s) {
        if (s == inapp_state) continue;
        if (all & (1u << s)) ++nc;
      }
      sp.n_states_remaining = nc;
      continue;
    }

    // Fully-unambiguous fast path. With no ambiguous tokens the classical
    // informativeness criterion is EXACT (there is no topology-dependent
    // ambiguity to resolve), so classify uninformative characters directly and
    // skip the transforms entirely. This also prevents Transform 2 from turning
    // an all-autapomorphy character (e.g. states 0..7, each in one tip) into a
    // wildcard-bearing token that the conservative ambiguous check can no
    // longer recognise as uninformative — the score would stay correct, but the
    // character would be needlessly rescored per tree. Informative unambiguous
    // characters fall through so genuine singleton autapomorphies are still
    // removed by Transform 2.
    if (!has_ambiguous_tips(sp.tip_tokens, n_tips, n_states, inapp_state)) {
      count_state_occurrences(sp.tip_tokens, n_tips, n_states, inapp_state,
                              unambig_count, any_count, unambig_tip_idx);
      if (is_uninformative_classical(unambig_count, n_states, inapp_state)) {
        sp.precomputed_steps += compute_fixed_steps(unambig_count, n_states,
                                                    inapp_state);
        sp.informative = false;
        sp.n_states_remaining = 0;
        result.n_patterns_removed++;
        result.total_offset_steps += sp.precomputed_steps;
        continue;
      }
    }

    // --- Transforms 2 & 3: iterate until stable ---
    // Transform 2 (singleton removal) can create situations where
    // Transform 3 (ambiguity removal) applies, and vice versa.
    bool changed_outer = true;
    while (changed_outer) {
      changed_outer = false;

      // Transform 3: Remove redundant ambiguity states
      // A state is redundant if it never appears unambiguously AND every
      // token containing it also contains at least one other state that
      // appears unambiguously somewhere.
      bool changed_t3 = true;
      while (changed_t3) {
        changed_t3 = false;
        count_state_occurrences(sp.tip_tokens, n_tips, n_states,
                                inapp_state, unambig_count, any_count,
                                unambig_tip_idx);
        for (int s = 0; s < n_states; ++s) {
          if (s == inapp_state) continue;
          if (any_count[s] == 0) continue;
          if (unambig_count[s] > 0) continue;
          // Sound ONLY when state s occurs in exactly one tip. A state present
          // in >=2 tips can form a cost-free clade on SOME tree (the optimal
          // resolution of an ambiguous token is topology-dependent), so
          // stripping it there changes the score. Simplification must be
          // topology-independent, so we may only drop a state that appears in a
          // single tip (a lone leaf option, provably free to remove).
          if (any_count[s] != 1) continue;

          bool safe_to_remove = true;
          for (int tip = 0; tip < n_tips; ++tip) {
            if (!(sp.tip_tokens[tip] & (1u << s))) continue;
            uint32_t other = sp.tip_tokens[tip] & ~(1u << s);
            if (inapp_state >= 0) other &= ~(1u << inapp_state);
            bool has_other_unambig = false;
            for (int s2 = 0; s2 < n_states; ++s2) {
              if (s2 == inapp_state || s2 == s) continue;
              if ((other & (1u << s2)) && unambig_count[s2] > 0) {
                has_other_unambig = true;
                break;
              }
            }
            if (!has_other_unambig) {
              safe_to_remove = false;
              break;
            }
          }

          if (safe_to_remove) {
            uint32_t mask = ~(1u << s);
            for (int tip = 0; tip < n_tips; ++tip) {
              sp.tip_tokens[tip] &= mask;
            }
            changed_t3 = true;
            changed_outer = true;
            result.n_states_reduced++;
          }
        }
      }

      // Transform 2: Remove singleton states
      // A state is a removable singleton if it appears unambiguously in
      // exactly 1 tip, AND every ambiguity token containing it also has
      // at least one other state with unambig_count >= 2 (so the
      // optimizer never needs to resolve the ambiguity to this state).
      bool changed_t2 = true;
      while (changed_t2) {
        changed_t2 = false;
        count_state_occurrences(sp.tip_tokens, n_tips, n_states,
                                inapp_state, unambig_count, any_count,
                                unambig_tip_idx);

        for (int s = 0; s < n_states; ++s) {
          if (s == inapp_state) continue;
          if (unambig_count[s] != 1) continue;
          // Sound ONLY when state s occurs in exactly one tip total, i.e. it is
          // a genuine autapomorphy (fixed +1 step on every tree). A singleton
          // that ALSO appears in ambiguous tokens is NOT fixed: the optimal
          // reconstruction may resolve an ambiguous neighbour TO s to form a
          // cost-free clade (e.g. tip {s} beside {s,x} join at zero cost), so
          // removing s and charging +1 over-counts. The previous `has_dominant`
          // guard only proved another state EXISTED, not that the optimiser
          // would avoid s — an unsound test that over-counted multistate
          // ambiguity. Ambiguity resolution is topology-dependent;
          // simplification must be topology-independent.
          if (any_count[s] != 1) continue;

          int tip = unambig_tip_idx[s];

          uint32_t all = all_applicable_mask(sp.tip_tokens, n_tips,
                                              n_states, inapp_state);
          uint32_t wildcard = all & ~(1u << s);
          if (wildcard == 0) continue;

          sp.tip_tokens[tip] = wildcard;
          sp.precomputed_steps += 1;
          changed_t2 = true;
          changed_outer = true;
          result.n_states_reduced++;
          break;
        }
      }
    }

    // --- Transform 1: Check informativeness ---
    count_state_occurrences(sp.tip_tokens, n_tips, n_states,
                            inapp_state, unambig_count, any_count,
                            unambig_tip_idx);

    bool uninformative = false;
    int fixed_steps = 0;

    if (is_uninformative_classical(unambig_count, n_states, inapp_state)) {
      if (has_ambiguous_tips(sp.tip_tokens, n_tips, n_states, inapp_state)) {
        // The classical criterion is unreliable with ambiguous tokens, and a
        // fixed cost > 0 CANNOT be proven cheaply: the optimal resolution of an
        // ambiguous token is topology-dependent. (The old caterpillar-sampling
        // heuristic was unsound — 4 orderings can agree on a cost while a
        // balanced tree scores differently, over- OR under-counting.) The one
        // provably topology-independent uninformative case is a state shared by
        // EVERY tip: then all tips resolve to it for 0 steps on every tree.
        // Otherwise keep the character; the per-tree downpass scores it exactly.
        uint32_t common = ~0u;
        for (int tip = 0; tip < n_tips; ++tip) {
          uint32_t tok = sp.tip_tokens[tip];
          if (inapp_state >= 0) tok &= ~(1u << inapp_state);
          common &= tok;
        }
        if (common != 0) {
          uninformative = true;
          fixed_steps = 0;
        }
        // else: uninformative stays false -> scored per-tree (sound).
      } else {
        // All tips unambiguous — classical criterion is correct.
        uninformative = true;
        fixed_steps = compute_fixed_steps(unambig_count, n_states,
                                           inapp_state);
      }
    }

    if (uninformative) {
      // Add the fixed steps of the reduced character to the accumulated
      // singleton steps from Transform 2
      sp.precomputed_steps += fixed_steps;
      sp.informative = false;
      sp.n_states_remaining = 0;
      result.n_patterns_removed++;
      result.total_offset_steps += sp.precomputed_steps;
    } else {
      // Count remaining states
      uint32_t all = all_applicable_mask(sp.tip_tokens, n_tips, n_states,
                                          inapp_state);
      int nc = 0;
      for (int s = 0; s < n_states; ++s) {
        if (s == inapp_state) continue;
        if (all & (1u << s)) ++nc;
      }
      sp.n_states_remaining = nc;
      result.total_offset_steps += sp.precomputed_steps;
    }
  }

  return result;
}

} // namespace ts
