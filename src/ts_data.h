#ifndef TS_DATA_H
#define TS_DATA_H

// Bit-packed character data for fast parsimony scoring.
//
// Characters are grouped into blocks of up to 64, where bit i of state-word j
// means "character i can be in state j". State 0 is the inapplicable (NA)
// state when has_inapplicable is true.
//
// Pattern weights are eliminated by expanding each pattern by its frequency,
// so all characters have weight 1 and popcount gives the score directly.

#include <cstdint>
#include <vector>

namespace ts {

// Portable popcount for uint64_t
inline int popcount64(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_popcountll(x);
#elif defined(_MSC_VER)
  return static_cast<int>(__popcnt64(x));
#else
  // Fallback: Hamming weight
  x = x - ((x >> 1) & 0x5555555555555555ULL);
  x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
  return static_cast<int>(
    (((x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL) * 0x0101010101010101ULL) >> 56);
#endif
}

static constexpr int MAX_CHARS_PER_BLOCK = 64;
static constexpr int MAX_STATES = 32;  // practical limit for morphological data

struct CharBlock {
  int n_chars;             // characters in this block (1..64)
  int n_states;            // number of states (including NA if has_inapplicable)
  bool has_inapplicable;   // state 0 is inapplicable → use NA-aware algorithm
  uint64_t active_mask;    // bits 0..n_chars-1 set, rest clear

  // For IW: map each character back to its original pattern index
  // (multiple characters may share the same pattern after weight expansion)
  int pattern_index[MAX_CHARS_PER_BLOCK];
};

struct DataSet {
  int n_tips;
  int n_blocks;
  int total_words;          // sum of n_states across all blocks

  std::vector<CharBlock> blocks;

  // Tip state data, flattened for cache locality:
  //   tip_states[tip * total_words + word_offset(block, state)]
  // where word_offset for block b, state s = block_word_offset[b] + s
  std::vector<uint64_t> tip_states;
  std::vector<int> block_word_offset;  // cumulative offset for each block

  // IW metadata (per original pattern)
  int n_patterns;                      // number of unique patterns
  std::vector<int> min_steps;          // minimum steps per pattern
  std::vector<int> pattern_freq;       // original weight (for reporting)
};

// Build a DataSet from R-side data.
//
// contrast_r: n_tokens x n_states matrix (doubles, 0/1) — the phyDat contrast
// tip_data_r: n_tips x n_patterns integer matrix — phyDat token indices (1-based)
// weight_r:   n_patterns integer vector — pattern frequencies
// levels_r:   character vector of state labels; "-" marks the inapplicable state
//
// Returns a fully populated DataSet with patterns expanded by weight.
DataSet build_dataset(
    const double* contrast_r, int n_tokens, int n_states,
    const int* tip_data_r, int n_tips, int n_patterns,
    const int* weight_r,
    const char** levels_r);

} // namespace ts

#endif // TS_DATA_H
