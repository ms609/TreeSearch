#include "ts_splits.h"
#include "ts_data.h"
#include <algorithm>
#include <cstring>

namespace ts {

// splitmix64-style mixer to reduce XOR cancellation
static uint64_t mix(uint64_t x) {
  x ^= x >> 30;
  x *= 0xbf58476d1ce4e5b9ULL;
  x ^= x >> 27;
  x *= 0x94d049bb133111ebULL;
  x ^= x >> 31;
  return x;
}

// Per-word prime multipliers for multi-word split hashing
static const uint64_t PRIMES[] = {
  0x9e3779b97f4a7c15ULL,
  0x517cc1b727220a95ULL,
  0x6c62272e07bb0142ULL,
  0x62b821756295c58dULL,
  0xcdb32970830fcaa1ULL,
  0xc1b6e8e4253e850fULL,
  0x3a39d80cf26f5e87ULL,
  0xf51f10c91a8e8a49ULL,
};

static uint64_t prime_for_word(int w) {
  if (w < 8) return PRIMES[w];
  // For words beyond our table, derive from mix
  return mix(static_cast<uint64_t>(w) * 0x9e3779b97f4a7c15ULL);
}

// Canonicalize: if bit 0 is set, complement all words.
// Mask the final word to n_tips bits.
static void canonicalize_split(uint64_t* s, int words_per_split, int n_tips) {
  bool flip = (s[0] & 1ULL) != 0;
  if (flip) {
    for (int w = 0; w < words_per_split; ++w) {
      s[w] = ~s[w];
    }
  }
  // Mask trailing bits in the last word
  int trailing = n_tips % 64;
  if (trailing != 0) {
    s[words_per_split - 1] &= (1ULL << trailing) - 1;
  }
}

SplitSet compute_splits(const TreeState& tree) {
  int n_tip = tree.n_tip;
  int wps = (n_tip + 63) / 64;

  // Temporary per-node bitsets: tip membership of each node's subtree
  size_t total = static_cast<size_t>(tree.n_node) * wps;
  std::vector<uint64_t> tip_bits(total, 0);

  // Initialize tips: tip i has bit i set
  for (int t = 0; t < n_tip; ++t) {
    int word = t / 64;
    int bit = t % 64;
    tip_bits[static_cast<size_t>(t) * wps + word] = 1ULL << bit;
  }

  // Postorder traversal: union children's bitsets at each internal node
  for (int pi = 0; pi < static_cast<int>(tree.postorder.size()); ++pi) {
    int node = tree.postorder[pi];
    int ni = node - n_tip;
    int lc = tree.left[ni];
    int rc = tree.right[ni];
    uint64_t* dst = &tip_bits[static_cast<size_t>(node) * wps];
    const uint64_t* lbits = &tip_bits[static_cast<size_t>(lc) * wps];
    const uint64_t* rbits = &tip_bits[static_cast<size_t>(rc) * wps];
    for (int w = 0; w < wps; ++w) {
      dst[w] = lbits[w] | rbits[w];
    }
  }

  // Collect non-trivial splits.
  // For each internal node except the root, its subtree bitset defines a split.
  // Skip trivial splits (single-tip subtrees are already excluded since we
  // only visit internal nodes). Also skip root's children — for an unrooted
  // tree the two children of the root produce complementary (redundant) splits,
  // so we skip one of them (the right child of the root).
  int root = n_tip; // root node index
  int root_right = tree.right[0]; // right child of root (root's ni = 0)

  // Count splits first
  int n_splits = 0;
  for (int pi = 0; pi < static_cast<int>(tree.postorder.size()); ++pi) {
    int node = tree.postorder[pi];
    if (node == root) continue;
    if (node == root_right) continue;
    // Check if this is a trivial split (single tip on one side)
    // A split is trivial if it separates exactly 1 tip or n_tip-1 tips.
    // Since we only visit internal nodes, each has ≥2 tips in its subtree.
    // But we need to check: does the subtree contain exactly 1 tip? No,
    // internal nodes always have ≥2 tips. So check if subtree has exactly
    // n_tip-1 tips (which would also be trivial). Actually for a fully
    // resolved binary tree, this can't happen for non-root internal nodes
    // (the complement would have ≥2 tips). So all internal node subtrees
    // except root's children complements are non-trivial.
    // Actually: an internal node could have n_tip-1 tips in its subtree
    // if it's a child of the root and the other child is a tip. In that
    // case the split separates 1 tip vs rest — trivial. We should skip those.
    const uint64_t* bits = &tip_bits[static_cast<size_t>(node) * wps];
    int count = 0;
    for (int w = 0; w < wps; ++w) {
      count += ts::popcount64(bits[w]);
    }
    if (count <= 1 || count >= n_tip - 1) continue;
    ++n_splits;
  }

  SplitSet ss;
  ss.n_tips = n_tip;
  ss.words_per_split = wps;
  ss.n_splits = n_splits;
  ss.splits.resize(static_cast<size_t>(n_splits) * wps);

  int idx = 0;
  for (int pi = 0; pi < static_cast<int>(tree.postorder.size()); ++pi) {
    int node = tree.postorder[pi];
    if (node == root) continue;
    if (node == root_right) continue;
    const uint64_t* bits = &tip_bits[static_cast<size_t>(node) * wps];
    int count = 0;
    for (int w = 0; w < wps; ++w) {
      count += ts::popcount64(bits[w]);
    }
    if (count <= 1 || count >= n_tip - 1) continue;

    uint64_t* dst = ss.split(idx);
    std::memcpy(dst, bits, sizeof(uint64_t) * wps);
    canonicalize_split(dst, wps, n_tip);
    ++idx;
  }

  return ss;
}

uint64_t hash_splits(const SplitSet& ss) {
  uint64_t h = 0;
  for (int i = 0; i < ss.n_splits; ++i) {
    const uint64_t* s = ss.split(i);
    uint64_t sh = 0;
    for (int w = 0; w < ss.words_per_split; ++w) {
      sh ^= s[w] * prime_for_word(w);
    }
    // XOR of mixed per-split hashes → order-independent
    h ^= mix(sh);
  }
  return h;
}

// Lexicographic comparison for sorting splits
static bool split_less(const uint64_t* a, const uint64_t* b, int wps) {
  for (int w = wps - 1; w >= 0; --w) {
    if (a[w] < b[w]) return true;
    if (a[w] > b[w]) return false;
  }
  return false;
}

bool splits_equal(const SplitSet& a, const SplitSet& b) {
  if (a.n_tips != b.n_tips || a.n_splits != b.n_splits) return false;
  if (a.n_splits == 0) return true;

  int wps = a.words_per_split;

  // Sort copies of both split vectors for comparison
  // Build index arrays and sort by split content
  std::vector<int> idx_a(a.n_splits), idx_b(b.n_splits);
  for (int i = 0; i < a.n_splits; ++i) idx_a[i] = i;
  for (int i = 0; i < b.n_splits; ++i) idx_b[i] = i;

  auto cmp_a = [&](int i, int j) {
    return split_less(a.split(i), a.split(j), wps);
  };
  auto cmp_b = [&](int i, int j) {
    return split_less(b.split(i), b.split(j), wps);
  };

  std::sort(idx_a.begin(), idx_a.end(), cmp_a);
  std::sort(idx_b.begin(), idx_b.end(), cmp_b);

  for (int k = 0; k < a.n_splits; ++k) {
    if (std::memcmp(a.split(idx_a[k]), b.split(idx_b[k]),
                    sizeof(uint64_t) * wps) != 0) {
      return false;
    }
  }
  return true;
}

} // namespace ts
