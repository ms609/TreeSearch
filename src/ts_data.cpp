#include "ts_data.h"
#include <algorithm>
#include <cstring>
#include <map>

namespace ts {

DataSet build_dataset(
    const double* contrast_r, int n_tokens, int n_states,
    const int* tip_data_r, int n_tips, int n_patterns,
    const int* weight_r,
    const char** levels_r)
{
  DataSet ds;
  ds.n_tips = n_tips;
  ds.n_patterns = n_patterns;

  // Identify the inapplicable state (column index where level == "-")
  int inapp_state = -1;
  for (int s = 0; s < n_states; ++s) {
    if (std::strcmp(levels_r[s], "-") == 0) {
      inapp_state = s;
      break;
    }
  }

  // For each pattern, determine:
  //   - which states are possible for each token
  //   - whether any taxon can be inapplicable (→ has_inapplicable block)
  //   - the number of applicable states (for grouping into blocks)

  // Precompute: for each token, the set of possible states as a bitmask.
  // Bit s set ⟹ state s is possible.
  std::vector<uint32_t> token_states(n_tokens, 0);
  for (int t = 0; t < n_tokens; ++t) {
    for (int s = 0; s < n_states; ++s) {
      // contrast_r is column-major: [t + n_tokens * s]
      if (contrast_r[t + n_tokens * s] > 0.5) {
        token_states[t] |= (1u << s);
      }
    }
  }

  // Classify each pattern: has_inapp + number of states
  struct PatternInfo {
    int pattern_idx;      // original 0-based pattern index
    bool has_inapp;       // any taxon can be inapplicable
    int n_applicable;     // number of applicable states used
    int weight;           // expansion count
  };

  std::vector<PatternInfo> patterns(n_patterns);
  for (int p = 0; p < n_patterns; ++p) {
    patterns[p].pattern_idx = p;
    patterns[p].weight = weight_r[p];
    patterns[p].has_inapp = false;
    uint32_t all_states = 0;
    for (int tip = 0; tip < n_tips; ++tip) {
      // tip_data_r is column-major: [tip + n_tips * p], 1-based token index
      int token = tip_data_r[tip + n_tips * p] - 1;
      all_states |= token_states[token];
      if (inapp_state >= 0 && (token_states[token] & (1u << inapp_state))) {
        patterns[p].has_inapp = true;
      }
    }
    // Count applicable states (exclude inapp_state)
    int n_app = 0;
    for (int s = 0; s < n_states; ++s) {
      if (s != inapp_state && (all_states & (1u << s))) ++n_app;
    }
    patterns[p].n_applicable = n_app;
  }

  // Expand patterns by weight and sort: first by has_inapp, then by n_states.
  // This groups characters that share a block type.
  struct ExpandedChar {
    int pattern_idx;
    bool has_inapp;
    int n_applicable;
  };

  std::vector<ExpandedChar> expanded;
  expanded.reserve(n_patterns * 2);  // rough guess
  for (int p = 0; p < n_patterns; ++p) {
    for (int w = 0; w < patterns[p].weight; ++w) {
      expanded.push_back({
        patterns[p].pattern_idx,
        patterns[p].has_inapp,
        patterns[p].n_applicable
      });
    }
  }

  // Sort: non-inapplicable first (inapplicable characters get their own blocks).
  // Within each partition, all blocks use the same number of state words
  // (the max observed across all characters in that partition) so the global
  // state_remap stays valid.
  std::stable_sort(expanded.begin(), expanded.end(),
    [](const ExpandedChar& a, const ExpandedChar& b) {
      return a.has_inapp < b.has_inapp;  // false (standard) before true (NA)
    });

  // Find the max n_applicable for each partition type
  int max_app_standard = 0, max_app_inapp = 0;
  for (const auto& ec : expanded) {
    if (ec.has_inapp) {
      if (ec.n_applicable > max_app_inapp) max_app_inapp = ec.n_applicable;
    } else {
      if (ec.n_applicable > max_app_standard) max_app_standard = ec.n_applicable;
    }
  }

  // Group into blocks of up to 64, all sharing the same has_inapp.
  ds.n_blocks = 0;
  ds.blocks.clear();

  int i = 0;
  int total_chars = static_cast<int>(expanded.size());
  while (i < total_chars) {
    bool block_inapp = expanded[i].has_inapp;
    int max_app = block_inapp ? max_app_inapp : max_app_standard;

    int block_size = 0;
    int start = i;
    while (i < total_chars && block_size < MAX_CHARS_PER_BLOCK &&
           expanded[i].has_inapp == block_inapp) {
      ++block_size;
      ++i;
    }

    CharBlock blk;
    blk.n_chars = block_size;
    blk.has_inapplicable = block_inapp;
    // All blocks in a partition use the max state count
    blk.n_states = max_app + (block_inapp ? 1 : 0);
    blk.active_mask = (block_size == 64) ? ~0ULL
                      : ((1ULL << block_size) - 1);
    for (int c = 0; c < block_size; ++c) {
      blk.pattern_index[c] = expanded[start + c].pattern_idx;
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

  // Build tip state data.
  // For each block, we need to map the original state indices to the
  // block's state word indices.
  //
  // If has_inapplicable: word 0 = NA state, words 1..n_app = applicable states
  // If !has_inapplicable: words 0..n_app-1 = applicable states
  //
  // We need a mapping: original_state_index → block_state_word for each block.
  // All characters in a block share the same applicable states (by construction),
  // so we build a global mapping of applicable states.

  // First, identify which original states are "applicable" (i.e., not the inapp state)
  // and assign them indices.
  std::vector<int> state_remap(n_states, -1);
  {
    int idx = 0;
    for (int s = 0; s < n_states; ++s) {
      if (s != inapp_state) {
        state_remap[s] = idx++;
      }
    }
    // inapp_state always maps to -1 (handled separately)
  }

  ds.tip_states.assign(
    static_cast<size_t>(n_tips) * ds.total_words, 0ULL);

  for (int b = 0; b < ds.n_blocks; ++b) {
    const CharBlock& blk = ds.blocks[b];
    int base = ds.block_word_offset[b];

    for (int c = 0; c < blk.n_chars; ++c) {
      int pat = blk.pattern_index[c];
      uint64_t bit = 1ULL << c;

      for (int tip = 0; tip < n_tips; ++tip) {
        int token = tip_data_r[tip + n_tips * pat] - 1;
        uint32_t tstates = token_states[token];

        int word_offset = base;  // start of this block's words for this tip
        size_t tip_base = static_cast<size_t>(tip) * ds.total_words;

        if (blk.has_inapplicable) {
          // Word 0 = NA
          if (inapp_state >= 0 && (tstates & (1u << inapp_state))) {
            ds.tip_states[tip_base + word_offset] |= bit;
          }
          // Words 1.. = applicable states in state_remap order
          for (int s = 0; s < n_states; ++s) {
            if (s == inapp_state) continue;
            if (tstates & (1u << s)) {
              int w = state_remap[s] + 1; // +1 because word 0 is NA
              ds.tip_states[tip_base + word_offset + w] |= bit;
            }
          }
        } else {
          // No NA; words 0.. = applicable states in state_remap order
          for (int s = 0; s < n_states; ++s) {
            if (s == inapp_state) continue;
            if (tstates & (1u << s)) {
              int w = state_remap[s];
              ds.tip_states[tip_base + word_offset + w] |= bit;
            }
          }
        }
      }
    }
  }

  // IW metadata
  ds.min_steps.resize(n_patterns, 0);
  ds.pattern_freq.resize(n_patterns, 0);
  for (int p = 0; p < n_patterns; ++p) {
    ds.pattern_freq[p] = weight_r[p];
    // min_steps will be computed later (requires tree structure knowledge
    // or separate calculation). Initialize to 0 for now.
  }

  return ds;
}

} // namespace ts
