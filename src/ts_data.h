#ifndef TS_DATA_H
#define TS_DATA_H

// Bit-packed character data for fast parsimony scoring.
//
// Characters are grouped into blocks of up to 64, where bit i of state-word j
// means "character i can be in state j". State 0 is the inapplicable (NA)
// state when has_inapplicable is true.
//
// Characters with the same weight are grouped into the same blocks.
// Block score = weight * popcount(needs_union), avoiding redundant expansion.

#include <cmath>
#include <cstdint>
#include <vector>
#include <unordered_set>

namespace ts {

// Hardware POPCNT via inline asm (no -mpopcnt flag needed).
// Same approach as TreeDist::popcnt64 — emits the instruction directly,
// avoiding the software Hamming weight fallback that __builtin_popcountll
// compiles to without -mpopcnt.
inline int popcount64(uint64_t x) {
#if (defined(__GNUC__) || defined(__clang__)) && defined(__x86_64__)
  uint64_t result;
  __asm__ ("popcnt %1, %0" : "=r" (result) : "r" (x));
  return static_cast<int>(result);
#elif defined(_MSC_VER) && defined(_M_X64)
  return static_cast<int>(__popcnt64(x));
#else
  // Fallback: software Hamming weight (non-x86-64 platforms)
  x = x - ((x >> 1) & 0x5555555555555555ULL);
  x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
  return static_cast<int>(
    (((x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL) * 0x0101010101010101ULL) >> 56);
#endif
}

// Portable count-trailing-zeros for uint64_t (undefined for x == 0)
inline int ctz64(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_ctzll(x);
#elif defined(_MSC_VER)
  unsigned long idx;
  _BitScanForward64(&idx, x);
  return static_cast<int>(idx);
#else
  // Fallback: de Bruijn sequence
  static const int debruijn[64] = {
     0,  1,  2, 53,  3,  7, 54, 27,  4, 38, 41,  8, 34, 55, 48, 28,
    62,  5, 39, 46, 44, 42, 22,  9, 24, 35, 59, 56, 49, 18, 29, 11,
    63, 52,  6, 26, 37, 40, 33, 47, 61, 45, 43, 21, 23, 58, 17, 10,
    51, 25, 36, 32, 60, 20, 57, 16, 50, 31, 19, 15, 30, 14, 13, 12
  };
  return debruijn[((x & -x) * 0x022FDD63CC95386DULL) >> 58];
#endif
}

static constexpr int MAX_CHARS_PER_BLOCK = 64;
static constexpr int MAX_STATES = 32;  // practical limit for morphological data

enum class ScoringMode { EW, IW, XPIWE, PROFILE, HSJ, XFORM };

// A hierarchy block describes one controlling primary + its secondaries
// (Hopkins & St. John 2021). Used by HSJ scoring.
struct HierarchyBlock {
  int primary_char;               // original character index (0-based)
  std::vector<int> secondary_chars; // original character indices (0-based)
  int n_secondaries;              // = secondary_chars.size()
  int absent_state;               // state index meaning "absent" in primary
};

struct CharBlock {
  int n_chars;             // characters in this block (1..64)
  int n_states;            // number of states (including NA if has_inapplicable)
  int weight;              // block weight (all chars in block share same weight)
  bool has_inapplicable;   // state 0 is inapplicable → use NA-aware algorithm
  uint64_t active_mask;    // bits 0..n_chars-1 set, rest clear

  // Ratchet upweighting: bits set here count double during perturbed scoring.
  // Must be a subset of active_mask. Default 0 (no upweighting).
  uint64_t upweight_mask = 0;

  // For IW: map each character back to its original pattern index
  // (multiple characters may share the same pattern after weight expansion)
  int pattern_index[MAX_CHARS_PER_BLOCK];
};

// Cache-friendly metadata for indirect scoring hot paths.
// Packs the 3 fields needed per block into 16 bytes (vs ~288 bytes in
// CharBlock). For 4 blocks this fits in a single 64-byte cache line.
struct FlatBlock {
  int offset;              // word offset into state arrays
  int n_states;            // states in this block (including NA if applicable)
  uint64_t active_mask;    // active character bits
  uint8_t has_inapplicable; // 1 if NA block (state 0 = inapplicable)
  uint8_t _pad[7];         // explicit padding to 24 bytes
}; // 24 bytes: 4 blocks = 96 bytes ≈ 1.5 cache lines (vs ~1152 for CharBlock)

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

  // Hot-path indirect scoring: cache-friendly block metadata.
  // Populated by build_dataset(); mirrors blocks[] + block_word_offset[].
  std::vector<FlatBlock> flat_blocks;
  // True when all blocks have weight == 1 (common EW case).
  // When true AND upweight_mask == 0, the specialized EW indirect
  // functions can skip per-block weight multiply and upweight checks.
  bool all_weight_one = false;

  // IW metadata (per original pattern)
  int n_patterns;                      // number of unique patterns
  std::vector<int> min_steps;          // minimum steps per pattern
  std::vector<int> pattern_freq;       // original weight (for reporting)
  double concavity;                    // IW concavity constant k; HUGE_VAL = EW

  // Extended IW (XPIWE): per-pattern effective concavity and Φ-rescaling.
  // Standard IW: eff_k[p] = concavity, phi[p] = 1.0 for all p.
  // XPIWE: eff_k[p] = concavity / f[p], phi[p] = (1+eff_k[p])/(1+concavity).
  // See Goloboff (2014) "Extended implied weighting", §missing entries.
  std::vector<double> eff_k;
  std::vector<double> phi;

  // Scoring mode (derived from concavity / info_amounts at build time)
  ScoringMode scoring_mode = ScoringMode::EW;

  // Profile parsimony lookup table (populated only when scoring_mode == PROFILE).
  // Column-major layout matching R: info_amounts[(step-1) + info_max_steps * pattern]
  // where step is the total step count (1-based) for that character.
  // Row 0 (= 1 total step = min steps for binary chars) has cost 0.
  std::vector<double> info_amounts;
  int info_max_steps = 0;              // number of rows in info_amounts

  // Character simplification metadata (populated by simplify_patterns).
  // ew_offset: sum of (precomputed_steps * weight) for all patterns,
  //   including removed uninformative ones. Added to EW score in score_tree().
  int ew_offset = 0;
  // Per-pattern step offset: topology-independent steps removed during
  //   simplification. Used by profile scoring to restore correct total steps,
  //   and by IW to adjust min_steps. Index by original pattern index.
  std::vector<int> precomputed_steps;

  // State index that represents the inapplicable ("-") state, or -1 if none.
  // Populated by build_dataset(); used by HSJ scoring.
  int inapp_state = -1;

  // HSJ scoring data (populated when scoring_mode == HSJ).
  // These are set by the Rcpp bridge after build_dataset().
  std::vector<HierarchyBlock> hierarchy_blocks;
  // tip_labels: per-tip per-original-char state labels (0-based).
  //   Layout: tip_labels[tip * n_orig_chars + char]
  std::vector<int> tip_labels;
  int n_orig_chars = 0;
  double hsj_alpha = 1.0;

  // Sankoff/xform scoring data (populated when scoring_mode == XFORM).
  // Each recoded hierarchy block becomes one Sankoff character.
  // Set by the Rcpp bridge after build_dataset().
  int sankoff_n_chars = 0;
  int sankoff_max_states = 0;
  std::vector<int> sankoff_n_states;           // [n_chars]
  std::vector<double> sankoff_cost_matrices;   // [n_chars * max_states * max_states]
  std::vector<int> sankoff_forced_root;        // [n_chars]
  // Flat tip costs: tip_costs[tip * stride + ch * max_states + state]
  // where stride = n_chars * max_states.
  // 0.0 = state allowed, INF = state disallowed.
  std::vector<double> sankoff_tip_costs;

  // Diagnostic counter: total TBR/SPR-class candidate rearrangements evaluated
  // (the analogue of TNT's "Total rearrangements examined"). Accumulated by
  // tbr_search() across an entire serial driven_search. `mutable` because the
  // kernels take `const DataSet&`. Valid only single-threaded: in parallel each
  // worker copies the DataSet (ts_parallel.cpp), so per-worker counts are not
  // aggregated. Excludes NNI-warmup and annealing candidates (neither funnels
  // through tbr_search).
  mutable long long n_candidates_evaluated = 0;

  // Exact-verify optimum memoization (NA path; see exact_verify_sweep in
  // ts_tbr.cpp).  A topology certified a true unrooted-TBR optimum under the
  // current weighting regime is cached here so repeated convergences — notably
  // across the ratchet's regime excursions — skip the O(n^2) full-neighbourhood
  // sweep.  Keyed by hash(child-pairs) ^ dataset-fp ^ weight-fp; cleared when
  // the dataset fingerprint changes.
  //
  // Lives on DataSet (NOT a function-local `static thread_local`) deliberately:
  // each parallel worker owns a private `ds_local` copy for its whole lifetime,
  // so this gives the same per-thread, cross-replicate persistence the old
  // thread_local had — but WITHOUT MinGW emutls, whose thread_local teardown
  // across std::thread spawn/exit corrupted the heap (parallel-only crash).
  // `mutable` because the scorer takes `const DataSet&`.  Single-writer per
  // copy: workers touch only their own ds_local; the shared prototype's cache
  // is written solely in the post-join (single-threaded) MPT phase.
  mutable std::unordered_set<uint64_t> evs_false_cache;
  mutable uint64_t evs_last_fp = 0;

  // Per-pattern step scratch for the weighted (IW/profile) full-rescore path
  // (fitch_score_ew).  Lives on DataSet for the SAME reason as evs_false_cache
  // above, NOT a function-local `static thread_local`: MinGW tears a
  // thread_local std::vector down via emutls when each std::thread worker
  // exits, and that teardown corrupted the heap across the parallel search's
  // repeated worker spawn/exit cycles (Windows-only hang in test-ts-parallel.R;
  // Linux native TLS is unaffected).  The per-worker `ds_local` copy gives the
  // same per-thread, cross-call capacity persistence the thread_local had.
  // `mutable` because the scorer takes `const DataSet&`; single-writer per copy.
  mutable std::vector<int> char_steps_scratch;
};

// Build a DataSet from R-side data.
//
// contrast_r: n_tokens x n_states matrix (doubles, 0/1) — the phyDat contrast
// tip_data_r: n_tips x n_patterns integer matrix — phyDat token indices (1-based)
// weight_r:   n_patterns integer vector — pattern frequencies
// levels_r:   character vector of state labels; "-" marks the inapplicable state
// obs_count_r: n_patterns integer vector — number of non-missing taxa per pattern
//              (used only when xpiwe = true; nullptr otherwise)
//
// Returns a fully populated DataSet with patterns expanded by weight.
DataSet build_dataset(
    const double* contrast_r, int n_tokens, int n_states,
    const int* tip_data_r, int n_tips, int n_patterns,
    const int* weight_r,
    const char** levels_r,
    const int* min_steps_r = nullptr,
    double concavity = HUGE_VAL,
    const double* info_amounts_r = nullptr,
    int info_max_steps = 0,
    bool xpiwe = false,
    double xpiwe_r = 0.5,
    double xpiwe_max_f = 5.0,
    const int* obs_count_r = nullptr);

} // namespace ts

#endif // TS_DATA_H
