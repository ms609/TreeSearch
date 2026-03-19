#include "ts_data.h"
#include "ts_simplify.h"
#include <R.h>
#include <algorithm>
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
    int info_max_steps)
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

  // Classify each pattern: has_inapp + number of applicable states + weight
  // Use simplified tokens to determine has_inapp and n_applicable.
  struct PatternInfo {
    int pattern_idx;
    bool has_inapp;
    int n_applicable;
    int weight;
  };

  std::vector<PatternInfo> patterns(n_patterns);
  for (int p = 0; p < n_patterns; ++p) {
    const auto& sp = simpl.patterns[p];
    patterns[p].pattern_idx = p;
    patterns[p].weight = weight_r[p];
    patterns[p].has_inapp = false;

    // Skip uninformative patterns (they're fully accounted for by ew_offset)
    if (!sp.informative) {
      patterns[p].weight = 0;  // will be removed below
      patterns[p].n_applicable = 0;
      continue;
    }

    uint32_t all_states = 0;
    for (int tip = 0; tip < n_tips; ++tip) {
      all_states |= sp.tip_tokens[tip];
      if (inapp_state >= 0 && (sp.tip_tokens[tip] & (1u << inapp_state))) {
        patterns[p].has_inapp = true;
      }
    }
    int n_app = 0;
    for (int s = 0; s < n_states; ++s) {
      if (s != inapp_state && (all_states & (1u << s))) ++n_app;
    }
    patterns[p].n_applicable = n_app;
  }

  // Remove zero-weight patterns before sorting — they contribute nothing
  // to scoring and would waste block space (especially after resampling).
  patterns.erase(
    std::remove_if(patterns.begin(), patterns.end(),
      [](const PatternInfo& p) { return p.weight == 0; }),
    patterns.end());

  // Sort by (has_inapp, weight desc) so characters with the same weight
  // and inapplicability status end up in the same blocks.
  // Descending weight puts expensive blocks first, improving early
  // termination in bounded indirect-length functions.
  std::stable_sort(patterns.begin(), patterns.end(),
    [](const PatternInfo& a, const PatternInfo& b) {
      if (a.has_inapp != b.has_inapp) return a.has_inapp < b.has_inapp;
      return a.weight > b.weight;
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

    CharBlock blk;
    blk.n_chars = block_size;
    blk.has_inapplicable = block_inapp;
    blk.weight = block_weight;
    blk.n_states = max_app + (block_inapp ? 1 : 0);
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

  // Mask with all applicable-state bits set (for missing-data detection)
  uint32_t all_applicable_mask = 0;
  for (int s = 0; s < n_states; ++s) {
    if (s != inapp_state) all_applicable_mask |= (1u << s);
  }

  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    int base = ds.block_word_offset[b];

    for (int c = 0; c < blk.n_chars; ++c) {
      int pat = blk.pattern_index[c];
      uint64_t bit = 1ULL << c;

      const auto& sp = simpl.patterns[pat];

      for (int tip = 0; tip < n_tips; ++tip) {
        // Use simplified tip tokens
        uint32_t tstates = sp.tip_tokens[tip];

        size_t tip_base = static_cast<size_t>(tip) * ds.total_words;

        if (blk.has_inapplicable) {
          bool has_inapp = inapp_state >= 0
                        && (tstates & (1u << inapp_state));
          uint32_t applicable_bits = tstates & all_applicable_mask;
          bool is_missing = has_inapp
                         && (applicable_bits == all_applicable_mask);
          if (has_inapp && !is_missing && applicable_bits != 0) {
            // Partial {-,X} ambiguity: strip applicable bits, treat as
            // pure inapplicable (MorphyLib convention). The three-pass
            // NA algorithm cannot correctly resolve partial {-,X}
            // ambiguity in tree context.
            ds.tip_states[tip_base + base] |= bit;
          } else if (has_inapp && applicable_bits == 0) {
            // Pure inapplicable: just set inapp bit
            ds.tip_states[tip_base + base] |= bit;
          } else {
            // Pure applicable OR full missing data: encode all bits
            if (has_inapp) {
              ds.tip_states[tip_base + base] |= bit;
            }
            for (int s = 0; s < n_states; ++s) {
              if (s == inapp_state) continue;
              if (tstates & (1u << s)) {
                int w = state_remap[s] + 1;
                ds.tip_states[tip_base + base + w] |= bit;
              }
            }
          }
        } else {
          for (int s = 0; s < n_states; ++s) {
            if (s == inapp_state) continue;
            if (tstates & (1u << s)) {
              int w = state_remap[s];
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
    ds.scoring_mode = ScoringMode::IW;
  } else {
    ds.scoring_mode = ScoringMode::EW;
  }

  return ds;
}

} // namespace ts
