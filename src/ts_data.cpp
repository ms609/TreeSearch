#include "ts_data.h"
#include "ts_simplify.h"
#include <R.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <map>

namespace ts {

DataSet build_dataset(
    const double* contrast_r, int n_tokens, int n_states,
    const int* tip_data_r, int n_tips, int n_patterns,
    const int* weight_r,
    const char** levels_r,
    const int* min_steps_r,
    double concavity,
    const double* info_amounts_r,
    int info_max_steps,
    bool xpiwe,
    double xpiwe_r,
    double xpiwe_max_f,
    const int* obs_count_r)
{
  DataSet ds;
  ds.n_tips = n_tips;
  ds.n_patterns = n_patterns;

  // Guard: state bitmasks use uint32_t, so n_states must fit in 32 bits.
  if (n_states > MAX_STATES) {
    Rf_error("TreeSearch C++ engine: n_states (%d) exceeds MAX_STATES (%d)",
             n_states, MAX_STATES);
  }

  // Identify the inapplicable state (column index where level == "-")
  int inapp_state = -1;
  for (int s = 0; s < n_states; ++s) {
    if (std::strcmp(levels_r[s], "-") == 0) {
      inapp_state = s;
      break;
    }
  }
  ds.inapp_state = inapp_state;

  // Precompute: for each token, the set of possible states as a bitmask.
  std::vector<uint32_t> token_states(n_tokens, 0);
  for (int t = 0; t < n_tokens; ++t) {
    for (int s = 0; s < n_states; ++s) {
      if (contrast_r[t + n_tokens * s] > 0.5) {
        token_states[t] |= (1u << s);
      }
    }
  }

  // --- Character simplification ---
  SimplificationResult simpl = simplify_patterns(
      token_states, tip_data_r, n_tips, n_patterns,
      weight_r, n_states, inapp_state);

  // Store per-pattern precomputed_steps and compute ew_offset
  ds.precomputed_steps.resize(n_patterns, 0);
  ds.ew_offset = 0;
  for (int p = 0; p < n_patterns; ++p) {
    ds.precomputed_steps[p] = simpl.patterns[p].precomputed_steps;
    ds.ew_offset += simpl.patterns[p].precomputed_steps * weight_r[p];
  }

  // Character cost-ordering mode (tertiary sort key after has_inapp, weight).
  // Concentrates high-homoplasy characters into the earliest blocks so the
  // bounded indirect-length scorers reach `cutoff` in fewer blocks on rejected
  // candidates (the converged-search majority) — a static, score-invariant
  // speedup on both the SPR-scan and reroot-x4 hot paths.  The accumulators sum
  // *homoplasy* (the min-steps floor is already folded into ew_offset), so the
  // measure estimates expected per-move homoplasy, NOT the step floor:
  //   MINORITY = informative tips off the plurality state (integer tally) —
  //              tracks expected step-count directly, folds in #states + balance
  //   ENTROPY  = Shannon entropy of applicable-state frequencies (nats) —
  //              information-theoretic sibling of MINORITY
  //   MIN_STEPS = #applicable states − 1 — the offset-out floor (weak; ties all
  //              2-state chars regardless of balance); kept for the bake-off
  //   NONE     = preserve input order within each (has_inapp, weight) group
  enum class CharOrder { NONE, MIN_STEPS, MINORITY, ENTROPY };
  // Default MINORITY: order characters within each (has_inapp, weight) group by
  // descending homoplasy so the bounded indirect-length scorers reach `cutoff`
  // in fewer blocks on rejected candidates (~3% fewer blocks on the EW x4 path,
  // measured 2026-07-13).  No valid move is ever dropped — a candidate below
  // cutoff is fully scanned and accepted regardless of order; only rejected
  // candidates bail early, and only their discarded partial changes.  The
  // search is stochastic, so a faster bail realigns the RNG stream to a
  // different-but-equally-valid path (reach neutral-to-better across seeds, not
  // a regression).  `minority` is an integer tally → reproducible across
  // platforms, unlike `entropy` whose log() risks ULP-level sort-key flips.
  CharOrder char_order = CharOrder::MINORITY;
  if (const char* e = std::getenv("TS_CHAR_ORDER")) {
    if (std::strcmp(e, "none") == 0)      char_order = CharOrder::NONE;
    else if (std::strcmp(e, "min_steps") == 0) char_order = CharOrder::MIN_STEPS;
    else if (std::strcmp(e, "minority") == 0)  char_order = CharOrder::MINORITY;
    else if (std::strcmp(e, "entropy") == 0)   char_order = CharOrder::ENTROPY;
  }

  // Classify each pattern: has_inapp + number of applicable states + weight
  // Use simplified tokens to determine has_inapp and n_applicable.
  struct PatternInfo {
    int pattern_idx;
    bool has_inapp;
    int n_applicable;
    int weight;
    double cost;   // homoplasy-ordering key (higher = scanned earlier); see above
    uint32_t alpha_mask;  // TS_PACK_LOCAL: applicable states the char can actually
                          // take (OR of NON-"?" tokens). "?" == this mask (Fitch-
                          // exact: a state absent from every tip can't be in any MPR),
                          // so a char using k<<n_states states needs only k planes.
  };

  // Full applicable-state mask (all applicable positions). A token equal to this
  // (all applicable bits) is a "?" / full-missing tip: it constrains nothing and is
  // excluded from a char's alpha_mask.
  uint32_t full_app_mask = 0;
  for (int s = 0; s < n_states; ++s)
    if (s != inapp_state) full_app_mask |= (1u << s);

  // TS_PACK_LOCAL (default OFF): per-block LOCAL state alphabets. Each block carries
  // only the union of its chars' actual state alphabets (not the global n_states),
  // cutting words gathered/reduced per candidate on low-state morphological data.
  // EXACT (Fitch length invariant under the "?"->char-alphabet substitution and a
  // bijective per-block state relabel). Scorer/combine are index-agnostic (loop
  // blk.n_states planes) so they are UNCHANGED; only the data build + plane indexing
  // change. Reconstruction/output map planes back to global labels via
  // CharBlock::plane_state[] (built below), so packed output is correct too.
  // DEFAULT-ON (2026-07-16): score-exact per tree, net-faster (1.1x small -> 1.5x
  // char-rich), reach-neutral across the corpus. Set TS_PACK_LOCAL=0 to force off.
  const char* pack_local_env = std::getenv("TS_PACK_LOCAL");
  const bool pack_local =
      (pack_local_env == nullptr) ||
      !(pack_local_env[0] == '0' && pack_local_env[1] == '\0');
  // TS_PACK_SORT: also cluster chars by alphabet (unlocks the union headroom but
  // reorders away from cost-ordering). Split from pack_local so we can measure
  // packing with cost-ordering PRESERVED (word-reduction isolated from the
  // reorder's bail cost).
  const bool pack_sort = pack_local && std::getenv("TS_PACK_SORT") != nullptr;

  std::vector<PatternInfo> patterns(n_patterns);
  for (int p = 0; p < n_patterns; ++p) {
    const auto& sp = simpl.patterns[p];
    patterns[p].pattern_idx = p;
    patterns[p].weight = weight_r[p];
    // Use the simplification-phase flag: only genuine "-" tokens trigger
    // inapplicable scoring (BGS three-pass). Characters where the inapp
    // bit only appeared in "?" (full missing data) use standard Fitch.
    patterns[p].has_inapp = sp.has_genuine_inapp;

    patterns[p].cost = 0.0;

    // Skip uninformative patterns (they're fully accounted for by ew_offset)
    if (!sp.informative) {
      patterns[p].weight = 0;  // will be removed below
      patterns[p].n_applicable = 0;
      patterns[p].alpha_mask = 0;
      continue;
    }

    // Single pass over tips: union of states (for n_applicable) plus a tally of
    // unambiguous applicable-state occurrences (for the homoplasy cost).  Only
    // single-state tips contribute to the tally — ambiguous/missing tips do not
    // force steps, so they add nothing to expected homoplasy.
    int state_count[MAX_STATES] = {0};
    const uint32_t inapp_bit = (inapp_state >= 0) ? (1u << inapp_state) : 0u;
    uint32_t all_states = 0;
    uint32_t alpha = 0;   // OR of applicable bits over NON-"?" (constraining) tips
    for (int tip = 0; tip < n_tips; ++tip) {
      uint32_t tok = sp.tip_tokens[tip];
      all_states |= tok;
      uint32_t app = tok & ~inapp_bit;
      if (app && app != full_app_mask) alpha |= app;  // "?" (== full_app_mask) constrains nothing
      if (app && (app & (app - 1)) == 0) {  // exactly one applicable bit set
        ++state_count[__builtin_ctz(app)];
      }
    }
    patterns[p].alpha_mask = alpha ? alpha : (all_states & ~inapp_bit);  // all-? fallback
    int n_app = 0;
    for (int s = 0; s < n_states; ++s) {
      if (s != inapp_state && (all_states & (1u << s))) ++n_app;
    }
    patterns[p].n_applicable = n_app;

    // Homoplasy-ordering cost (higher = scanned earlier).  See CharOrder above.
    double cost = 0.0;
    switch (char_order) {
      case CharOrder::NONE:
        break;
      case CharOrder::MIN_STEPS:
        cost = (n_app > 0) ? (n_app - 1) : 0;  // offset-out floor (weak)
        break;
      case CharOrder::MINORITY: {
        int total = 0, maxc = 0;
        for (int s = 0; s < n_states; ++s) {
          if (s == inapp_state) continue;
          total += state_count[s];
          if (state_count[s] > maxc) maxc = state_count[s];
        }
        cost = total - maxc;  // informative tips off the plurality state
        break;
      }
      case CharOrder::ENTROPY: {
        int total = 0;
        for (int s = 0; s < n_states; ++s) {
          if (s != inapp_state) total += state_count[s];
        }
        if (total > 0) {
          const double invN = 1.0 / total;
          for (int s = 0; s < n_states; ++s) {
            if (s == inapp_state) continue;
            const int c = state_count[s];
            if (c > 0) {
              const double pr = c * invN;
              cost -= pr * std::log(pr);  // Shannon entropy (nats)
            }
          }
        }
        break;
      }
    }
    patterns[p].cost = cost;
  }

  // Per-pattern alphabet lookup (indexed by ORIGINAL pattern idx) for the
  // TS_PACK_LOCAL tip fill's "?"->char-alphabet substitution. Built before the
  // erase/sort while patterns[p].pattern_idx == p.
  std::vector<uint32_t> pat_alpha(n_patterns, 0);
  for (int p = 0; p < n_patterns; ++p) pat_alpha[p] = patterns[p].alpha_mask;

  // Remove zero-weight patterns before sorting — they contribute nothing
  // to scoring and would waste block space (especially after resampling).
  patterns.erase(
    std::remove_if(patterns.begin(), patterns.end(),
      [](const PatternInfo& p) { return p.weight == 0; }),
    patterns.end());

  // Sort by (has_inapp, weight desc, cost desc) so characters with the same
  // weight and inapplicability status end up in the same blocks.  Descending
  // weight puts expensive blocks first; the cost tertiary key then packs the
  // highest-homoplasy characters of each (has_inapp, weight) group into the
  // earliest blocks, so the bounded indirect-length scorers reach `cutoff` in
  // fewer blocks on rejected candidates.  In the common equal-weights case
  // (all weight 1) weight is inert and cost is the sole discriminator.
  // stable_sort keeps input order among equal-cost characters (determinism).
  std::stable_sort(patterns.begin(), patterns.end(),
    [pack_sort](const PatternInfo& a, const PatternInfo& b) {
      if (a.has_inapp != b.has_inapp) return a.has_inapp < b.has_inapp;
      if (a.weight != b.weight) return a.weight > b.weight;
      if (pack_sort && a.alpha_mask != b.alpha_mask) {
        // Cluster chars by their state alphabet so each block's union stays small.
        // Group by alphabet SIZE first (few-state chars together), then by the mask
        // value (identical alphabets adjacent). Cost is the within-group tiebreaker,
        // so cost-ordering is preserved inside each alphabet cluster.
        int pa = __builtin_popcount(a.alpha_mask), pb = __builtin_popcount(b.alpha_mask);
        if (pa != pb) return pa < pb;
        return a.alpha_mask < b.alpha_mask;
      }
      return a.cost > b.cost;
    });

  // Count total applicable states in the dataset.
  // All blocks must use the same number of applicable state words because
  // the state_remap assigns globally consecutive indices. A pattern using
  // state index k needs state word k, regardless of how many states that
  // individual pattern uses.
  int total_app_states = 0;
  for (int s = 0; s < n_states; ++s) {
    if (s != inapp_state) ++total_app_states;
  }
  int max_app_standard = total_app_states;
  int max_app_inapp = total_app_states;

  // Group into blocks: same has_inapp AND same weight, up to 64 per block.
  ds.n_blocks = 0;
  ds.blocks.clear();

  // TS_PACK_LOCAL: per-block union of applicable-state alphabets (global bit
  // positions present in the block). Plane p of the block = the p-th set bit.
  // Empty when pack_local is off. Consumed by the tip_states fill below.
  std::vector<uint32_t> block_app_mask;

  int i_pat = 0;
  int total_patterns_active = static_cast<int>(patterns.size());
  while (i_pat < total_patterns_active) {
    bool block_inapp = patterns[i_pat].has_inapp;
    int block_weight = patterns[i_pat].weight;
    int max_app = block_inapp ? max_app_inapp : max_app_standard;

    int block_size = 0;
    int start = i_pat;
    while (i_pat < total_patterns_active && block_size < MAX_CHARS_PER_BLOCK &&
           patterns[i_pat].has_inapp == block_inapp &&
           patterns[i_pat].weight == block_weight) {
      ++block_size;
      ++i_pat;
    }

    // Local alphabet = union of the block's chars' alpha_masks (applicable states).
    uint32_t blk_app = 0;
    for (int c = 0; c < block_size; ++c) blk_app |= patterns[start + c].alpha_mask;
    int local_app = __builtin_popcount(blk_app);

    CharBlock blk;
    blk.n_chars = block_size;
    blk.has_inapplicable = block_inapp;
    blk.weight = block_weight;
    blk.n_states = (pack_local ? local_app : max_app) + (block_inapp ? 1 : 0);
    block_app_mask.push_back(blk_app);
    blk.active_mask = (block_size == 64) ? ~0ULL
                      : ((1ULL << block_size) - 1);
    for (int c = 0; c < block_size; ++c) {
      blk.pattern_index[c] = patterns[start + c].pattern_idx;
    }
    for (int c = block_size; c < MAX_CHARS_PER_BLOCK; ++c) {
      blk.pattern_index[c] = -1;
    }
    ds.blocks.push_back(blk);
    ++ds.n_blocks;
  }

  // Compute word offsets
  ds.block_word_offset.resize(ds.n_blocks);
  ds.total_words = 0;
  for (int b = 0; b < ds.n_blocks; ++b) {
    ds.block_word_offset[b] = ds.total_words;
    ds.total_words += ds.blocks[b].n_states;
  }

  // Pad total_words to even count for SIMD safety (SSE2 loads 2 × uint64_t).
  // Padding words are zero-initialized and don't affect scoring.
  if (ds.total_words % 2 != 0) {
    ds.total_words += 1;
  }

  // Build cache-friendly flat block metadata for indirect scoring hot paths.
  ds.flat_blocks.resize(ds.n_blocks);
  ds.all_weight_one = true;
  for (int b = 0; b < ds.n_blocks; ++b) {
    ds.flat_blocks[b].offset = ds.block_word_offset[b];
    ds.flat_blocks[b].n_states = ds.blocks[b].n_states;
    ds.flat_blocks[b].active_mask = ds.blocks[b].active_mask;
    ds.flat_blocks[b].has_inapplicable = ds.blocks[b].has_inapplicable ? 1 : 0;
    std::memset(ds.flat_blocks[b]._pad, 0, sizeof(ds.flat_blocks[b]._pad));
    if (ds.blocks[b].weight != 1) ds.all_weight_one = false;
  }

  // Build state-to-word mapping (applicable states only, excluding inapp)
  std::vector<int> state_remap(n_states, -1);
  {
    int idx = 0;
    for (int s = 0; s < n_states; ++s) {
      if (s != inapp_state) {
        state_remap[s] = idx++;
      }
    }
  }

  // Build tip state data — use simplified tokens where available
  ds.tip_states.assign(
    static_cast<size_t>(n_tips) * ds.total_words, 0ULL);

  const uint32_t inapp_bit = (inapp_state >= 0) ? (1u << inapp_state) : 0u;
  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    int base = ds.block_word_offset[b];

    // Applicable-state -> local plane index for THIS block. Default: the global
    // state_remap (identity across blocks). TS_PACK_LOCAL: the rank of the state
    // among the block's local alphabet (block_app_mask[b]) — so a block spanning
    // few states uses few planes.
    int loc[MAX_STATES];
    for (int s = 0; s < n_states; ++s) loc[s] = (s == inapp_state) ? -1 : state_remap[s];
    if (pack_local) {
      uint32_t bm = block_app_mask[b]; int idx = 0;
      for (int s = 0; s < n_states; ++s) {
        if (s == inapp_state) { loc[s] = -1; continue; }
        loc[s] = (bm & (1u << s)) ? idx++ : -1;
      }
    }

    // Invert loc[] into the display map: applicable plane p shows the global
    // applicable label state_remap[s] of whichever state occupies it. Identity
    // in the unpacked layout (loc[s]==state_remap[s]); under packing it maps the
    // p-th local-alphabet plane back to its global label for reconstruction.
    for (int p = 0; p < MAX_STATES; ++p) ds.blocks[b].plane_state[p] = p;
    for (int s = 0; s < n_states; ++s) {
      if (loc[s] >= 0) ds.blocks[b].plane_state[loc[s]] = state_remap[s];
    }

    for (int c = 0; c < blk.n_chars; ++c) {
      int pat = blk.pattern_index[c];
      uint64_t bit = 1ULL << c;

      const auto& sp = simpl.patterns[pat];

      for (int tip = 0; tip < n_tips; ++tip) {
        // Use simplified tip tokens
        uint32_t tstates = sp.tip_tokens[tip];

        // TS_PACK_LOCAL: a "?" (all applicable bits set) maps to the CHAR's own
        // alphabet, not the global all-states mask (Fitch-exact: a state absent
        // from every real tip is never in an MPR). This is what lets a low-state
        // char occupy few planes even under heavy missing data.
        if (pack_local && (tstates & full_app_mask) == full_app_mask)
          tstates = (tstates & inapp_bit) | pat_alpha[pat];

        size_t tip_base = static_cast<size_t>(tip) * ds.total_words;

        if (blk.has_inapplicable) {
          // Encode {-,X} partial ambiguity faithfully: set the NA (word 0)
          // bit AND the applicable-state bit(s).  Under BGS the applicability
          // character (binary Fitch, 0=inapplicable/1=applicable) is
          // reconstructed parsimoniously — resolving each {-,X} tip in tree
          // context, applicable-preferred on ties (De Laet: maximise homology,
          // not minimise homoplasy) — and the three passes then count state
          // steps.
          //   pure inapplicable -> {NA};      {-,X} -> {NA, X};
          //   missing {-,all}   -> {NA, all applicable} (full wildcard);
          //   pure applicable   -> applicable bits only.
          if (inapp_state >= 0 && (tstates & (1u << inapp_state))) {
            ds.tip_states[tip_base + base] |= bit;
          }
          for (int s = 0; s < n_states; ++s) {
            if (s == inapp_state) continue;
            if ((tstates & (1u << s)) && loc[s] >= 0) {
              int w = loc[s] + 1;
              ds.tip_states[tip_base + base + w] |= bit;
            }
          }
        } else {
          for (int s = 0; s < n_states; ++s) {
            if (s == inapp_state) continue;
            if ((tstates & (1u << s)) && loc[s] >= 0) {
              int w = loc[s];
              ds.tip_states[tip_base + base + w] |= bit;
            }
          }
        }
      }
    }
  }

  // IW metadata — adjust min_steps by precomputed_steps offset
  ds.min_steps.resize(n_patterns, 0);
  ds.pattern_freq.resize(n_patterns, 0);
  ds.concavity = concavity;
  for (int p = 0; p < n_patterns; ++p) {
    ds.pattern_freq[p] = weight_r[p];
    if (min_steps_r) {
      ds.min_steps[p] = min_steps_r[p] - ds.precomputed_steps[p];
      if (ds.min_steps[p] < 0) ds.min_steps[p] = 0;
    }
  }

  // Profile parsimony: copy info_amounts table and set scoring mode.
  // info_amounts_r is column-major from R: [max_steps × n_patterns].
  if (info_amounts_r && info_max_steps > 0) {
    size_t table_size = static_cast<size_t>(info_max_steps) * n_patterns;
    ds.info_amounts.assign(info_amounts_r, info_amounts_r + table_size);
    ds.info_max_steps = info_max_steps;
    ds.scoring_mode = ScoringMode::PROFILE;
    // Set concavity to a finite value so isfinite() checks in search
    // modules activate the weighted (indirect IW) pipeline.
    ds.concavity = 1.0;
  } else if (std::isfinite(concavity)) {
    ds.scoring_mode = (xpiwe && obs_count_r) ? ScoringMode::XPIWE
                                             : ScoringMode::IW;
  } else {
    ds.scoring_mode = ScoringMode::EW;
  }

  // Populate per-pattern eff_k and phi (Goloboff 2014, §missing entries).
  ds.eff_k.resize(n_patterns);
  ds.phi.resize(n_patterns);
  if (ds.scoring_mode == ScoringMode::XPIWE) {
    for (int p = 0; p < n_patterns; ++p) {
      // Goloboff (2014) Extension 3, verified against TNT 1.6:
      // f = 1 + r * missing / obs  (NOT r * total / obs)
      int obs = obs_count_r[p];
      int missing = n_tips - obs;
      // An all-missing pattern (obs == 0) has no observed tips and so no
      // steps to weight; guard the division (obs in the denominator would give
      // f = Inf, capped to xpiwe_max_f -> a spurious finite penalty). Neutral
      // eff_k/phi leave its (zero) step count unweighted.
      if (obs == 0) {
        ds.eff_k[p] = concavity;
        ds.phi[p] = 1.0;
        continue;
      }
      double f = 1.0 + xpiwe_r * missing / static_cast<double>(obs);
      if (f < 1.0) f = 1.0;
      if (f > xpiwe_max_f) f = xpiwe_max_f;
      ds.eff_k[p] = concavity / f;
      // Φ = w(1, k_ref) / w(1, k_c) where w(1, k) = 1/(1+k)
      ds.phi[p] = (1.0 + ds.eff_k[p]) / (1.0 + concavity);
    }
  } else {
    std::fill(ds.eff_k.begin(), ds.eff_k.end(),
              std::isfinite(concavity) ? concavity : 0.0);
    std::fill(ds.phi.begin(), ds.phi.end(), 1.0);
  }

  return ds;
}

} // namespace ts
