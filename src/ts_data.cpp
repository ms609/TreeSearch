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

  // Precompute: for each token, the set of possible states as a bitmask.
  std::vector<uint32_t> token_states(n_tokens, 0);
  for (int t = 0; t < n_tokens; ++t) {
    for (int s = 0; s < n_states; ++s) {
      if (contrast_r[t + n_tokens * s] > 0.5) {
        token_states[t] |= (1u << s);
      }
    }
  }

  // Classify each pattern: has_inapp + number of applicable states + weight
  struct PatternInfo {
    int pattern_idx;
    bool has_inapp;
    int n_applicable;
    int weight;
  };

  std::vector<PatternInfo> patterns(n_patterns);
  for (int p = 0; p < n_patterns; ++p) {
    patterns[p].pattern_idx = p;
    patterns[p].weight = weight_r[p];
    patterns[p].has_inapp = false;
    uint32_t all_states = 0;
    for (int tip = 0; tip < n_tips; ++tip) {
      int token = tip_data_r[tip + n_tips * p] - 1;
      all_states |= token_states[token];
      if (inapp_state >= 0 && (token_states[token] & (1u << inapp_state))) {
        patterns[p].has_inapp = true;
      }
    }
    int n_app = 0;
    for (int s = 0; s < n_states; ++s) {
      if (s != inapp_state && (all_states & (1u << s))) ++n_app;
    }
    patterns[p].n_applicable = n_app;
  }

  // Sort by (has_inapp, weight) so characters with the same weight
  // and inapplicability status end up in the same blocks.
  // No weight expansion: each pattern appears exactly once.
  std::stable_sort(patterns.begin(), patterns.end(),
    [](const PatternInfo& a, const PatternInfo& b) {
      if (a.has_inapp != b.has_inapp) return a.has_inapp < b.has_inapp;
      return a.weight < b.weight;
    });

  // Find the max n_applicable for each partition type
  int max_app_standard = 0, max_app_inapp = 0;
  for (const auto& pi : patterns) {
    if (pi.has_inapp) {
      if (pi.n_applicable > max_app_inapp) max_app_inapp = pi.n_applicable;
    } else {
      if (pi.n_applicable > max_app_standard) max_app_standard = pi.n_applicable;
    }
  }

  // Group into blocks: same has_inapp AND same weight, up to 64 per block.
  ds.n_blocks = 0;
  ds.blocks.clear();

  int i_pat = 0;
  int total_patterns = static_cast<int>(patterns.size());
  while (i_pat < total_patterns) {
    bool block_inapp = patterns[i_pat].has_inapp;
    int block_weight = patterns[i_pat].weight;
    int max_app = block_inapp ? max_app_inapp : max_app_standard;

    int block_size = 0;
    int start = i_pat;
    while (i_pat < total_patterns && block_size < MAX_CHARS_PER_BLOCK &&
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

  // Build tip state data
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

        size_t tip_base = static_cast<size_t>(tip) * ds.total_words;

        if (blk.has_inapplicable) {
          if (inapp_state >= 0 && (tstates & (1u << inapp_state))) {
            ds.tip_states[tip_base + base] |= bit;
          }
          for (int s = 0; s < n_states; ++s) {
            if (s == inapp_state) continue;
            if (tstates & (1u << s)) {
              int w = state_remap[s] + 1;
              ds.tip_states[tip_base + base + w] |= bit;
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

  // IW metadata
  ds.min_steps.resize(n_patterns, 0);
  ds.pattern_freq.resize(n_patterns, 0);
  for (int p = 0; p < n_patterns; ++p) {
    ds.pattern_freq[p] = weight_r[p];
  }

  return ds;
}

} // namespace ts
