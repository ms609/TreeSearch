#ifndef TS_CID_H
#define TS_CID_H

// CID (Clustering Information Distance) scoring for the driven search.
//
// Scores a candidate tree against a set of precomputed input tree splits
// using mutual clustering information (MCI) with a LAP solver.
//
// MRP (Matrix Representation with Parsimony) characters provide fast
// incremental screening for TBR candidate selection; CID is computed
// via score_tree() for verification only.
//
// The MCI algorithm is adapted from TreeDist (Smith 2020).
// When TreeDist exposes LinkingTo headers, this local copy should be
// replaced with #include <TreeDist/mutual_clustering.h>.

#include "ts_data.h"
#include "ts_tree.h"
#include "ts_splits.h"
#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <unordered_map>

namespace ts {

// --------------------------------------------------------------------------
// CID split set: precomputed splits for one tree
// --------------------------------------------------------------------------
struct CidSplitSet {
  int n_splits = 0;
  int n_bins = 0;                // ceil(n_tips / 64)
  std::vector<uint64_t> data;    // n_splits * n_bins, contiguous
  std::vector<int> in_split;     // popcount per split

  // Precomputed log2 values per split (populated by prepare_cid_data)
  std::vector<double> lg2_in;    // log2(in_split[i])
  std::vector<double> lg2_out;   // log2(n_tips - in_split[i])

  // Hash index for O(1) exact-match lookup (populated by prepare_cid_data).
  // Maps split hash -> split index.  Only built for input tree split sets.
  std::unordered_map<uint64_t, int> hash_index;

  const uint64_t* split(int i) const {
    return &data[static_cast<size_t>(i) * n_bins];
  }
};

// --------------------------------------------------------------------------
// LAP cost matrix (simplified from TreeDist's CostMatrix)
// --------------------------------------------------------------------------
// Integer cost matrix for the Jonker-Volgenant LAP solver.
// Row-major: element (i,j) at data_[i * dim_padded_ + j].
using lap_cost = int_fast64_t;
static constexpr lap_cost LAP_BIG =
    std::numeric_limits<lap_cost>::max() / 4096;
static constexpr size_t LAP_BLOCK = 16;

struct LapMatrix {
  int dim_ = 0;
  int dim_padded_ = 0;
  std::vector<lap_cost> data_;

  LapMatrix() = default;

  void resize(int dim) {
    dim_ = dim;
    dim_padded_ = static_cast<int>(
        ((static_cast<size_t>(dim) + LAP_BLOCK - 1) / LAP_BLOCK) * LAP_BLOCK);
    size_t needed = static_cast<size_t>(dim_padded_) * dim_padded_;
    if (data_.size() < needed) data_.resize(needed, 0);
  }

  lap_cost& operator()(int r, int c) {
    return data_[static_cast<size_t>(r) * dim_padded_ + c];
  }
  const lap_cost& operator()(int r, int c) const {
    return data_[static_cast<size_t>(r) * dim_padded_ + c];
  }
  lap_cost* row(int r) {
    return &data_[static_cast<size_t>(r) * dim_padded_];
  }
  const lap_cost* row(int r) const {
    return &data_[static_cast<size_t>(r) * dim_padded_];
  }

  void pad_row_after(int r, int start_col, lap_cost val) {
    size_t off = static_cast<size_t>(r) * dim_padded_;
    std::fill(data_.begin() + off + start_col,
              data_.begin() + off + dim_, val);
  }
  void pad_after_row(int start_row, lap_cost val) {
    size_t off = static_cast<size_t>(start_row) * dim_padded_;
    std::fill(data_.begin() + off,
              data_.begin() + static_cast<size_t>(dim_) * dim_padded_, val);
  }
};

// Reusable scratch for LAP solver
struct LapScratch {
  std::vector<lap_cost> v;
  std::vector<int> matches;
  std::vector<int> free_rows;
  std::vector<int> col_list;
  std::vector<lap_cost> d;
  std::vector<int> pred;
  std::vector<int> rowsol;
  std::vector<int> colsol;
  LapMatrix score_pool;    // reusable cost matrix

  void ensure(int dim) {
    int padded = static_cast<int>(
        ((static_cast<size_t>(dim) + LAP_BLOCK - 1) / LAP_BLOCK) * LAP_BLOCK);
    if (static_cast<int>(v.size()) < padded)      v.resize(padded);
    if (static_cast<int>(matches.size()) < dim)    matches.resize(dim);
    if (static_cast<int>(free_rows.size()) < dim)  free_rows.resize(dim);
    if (static_cast<int>(col_list.size()) < dim)   col_list.resize(dim);
    if (static_cast<int>(d.size()) < dim)          d.resize(dim);
    if (static_cast<int>(pred.size()) < dim)       pred.resize(dim);
    if (static_cast<int>(rowsol.size()) < dim)     rowsol.resize(dim);
    if (static_cast<int>(colsol.size()) < dim)     colsol.resize(dim);
  }
};

// --------------------------------------------------------------------------
// CidData: all precomputed data for CID scoring
// --------------------------------------------------------------------------
struct CidData {
  int n_trees;
  int n_tips;
  int n_bins;                    // ceil(n_tips / 64)

  // Precomputed input tree splits
  std::vector<CidSplitSet> tree_splits;

  // Precomputed per-tree clustering entropies
  std::vector<double> tree_ce;
  double mean_tree_ce;

  // Per-tree weights (1.0 normally; modified during ratchet)
  std::vector<double> tree_weights;
  double weight_sum;

  // Saved weights for ratchet restore
  std::vector<double> saved_tree_weights;
  double saved_weight_sum;

  // Normalized scoring mode (1 - mean(MCI_i / CE_i))
  bool normalize;

  // MRP screening parameters.
  double mrp_concavity = 7.0;
  double screening_tolerance = 0.0;

  // Block boundaries: mrp_tree_block_start[i] = first CharBlock index
  // for input tree i. mrp_tree_block_start[n_trees] = total blocks.
  std::vector<int> mrp_tree_block_start;

  // --- Precomputed constants (populated by prepare_cid_data) ---
  double lg2_n = 0.0;           // log2(n_tips)
  int max_splits = 0;           // max n_splits across all input trees

  // --- Persistent candidate buffers (reused across cid_score calls) ---
  mutable std::vector<uint64_t> cand_tip_bits;
  mutable CidSplitSet cand_buf;
  mutable std::vector<bool> match_a, match_b;

  // Score budget for early termination.  Set < HUGE_VAL to enable.
  mutable double score_budget = HUGE_VAL;

  // LAP scratch (reused across cid_score calls)
  mutable LapScratch lap_scratch;
};

// --------------------------------------------------------------------------
// Public API
// --------------------------------------------------------------------------

// Score a candidate tree against input trees using CID.
double cid_score(TreeState& tree, const CidData& cd);

// Build an MRP DataSet from CidData.
DataSet build_mrp_dataset(CidData& cd);

// Prepare CidData after tree_splits are populated: build hash indices,
// precompute log2 values, presize scratch buffers.
void prepare_cid_data(CidData& cd);

// Convert a SplitSet (from ts_splits.h) to a CidSplitSet.
CidSplitSet splitset_to_cid(const SplitSet& ss, int n_tips);

// Compute splits directly into a preallocated CidSplitSet.
void compute_splits_cid(const TreeState& tree,
                        std::vector<uint64_t>& tip_bits_work,
                        CidSplitSet& out);

// Compute clustering entropy of a split set.
double clustering_entropy(const CidSplitSet& ss, int n_tips);

// Compute clustering entropy using precomputed log2 values.
double clustering_entropy_fast(const CidSplitSet& ss, int n_tips,
                               double lg2_n);

// Compute mutual clustering information between two split sets.
double mutual_clustering_info(const CidSplitSet& a, const CidSplitSet& b,
                              int n_tips, LapScratch& scratch);

// Jonker-Volgenant LAP solver.
lap_cost lap_solve(int dim, LapMatrix& cost,
                   int* rowsol, int* colsol, LapScratch& scratch);

// Save / restore CidData tree weights for ratchet perturbation.
void save_cid_weights(CidData& cd);
void restore_cid_weights(CidData& cd);

// Synchronize CidData tree weights with MRP block perturbation.
void sync_cid_weights_from_mrp(CidData& cd, const DataSet& ds);

} // namespace ts

#endif // TS_CID_H
