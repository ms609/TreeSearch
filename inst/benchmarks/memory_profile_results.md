# Phase 3D: Memory Layout Profiling Results

Date: 2026-03-16
Platform: Windows, R 4.5.2, GCC 14.2.0
CPU: Intel (L1 32 KB, L2 256 KB typical)

## 1. Baseline Measurements

### TBR pass phase breakdown

All timings in microseconds (μs), averaged over 3 random trees per dataset.

| Dataset | Tips | Blocks | Words | Clips | Candidates | Clip+Incr (μs) | Indirect (μs) | Unclip (μs) |
|---------|------|--------|-------|-------|------------|-----------------|----------------|-------------|
| Vinther2008 | 23 | 6 | 28 | 38 | 3,585 | 789 | 286 | 268 |
| Agnarsson2004 | 62 | 8 | 59 | 112 | 56,501 | 2,948 | 5,175 | 856 |
| synth_20 | 20 | 4 | 11 | 34 | 2,535 | 271 | 65 | 93 |
| synth_50 | 50 | 4 | 12 | 91 | 32,776 | 1,021 | 989 | 314 |
| synth_100 | 100 | 4 | 12 | 190 | 237,536 | 3,880 | 7,999 | 1,013 |
| synth_200 | 200 | 4 | 12 | 377 | 1,090,533 | 11,238 | 35,930 | 2,695 |

### Time fraction breakdown

| Dataset | Tips | % Clip+Incr | % Indirect | % Unclip |
|---------|------|-------------|------------|----------|
| synth_20 | 20 | 63.2 | 15.1 | 21.7 |
| synth_50 | 50 | 43.9 | 42.6 | 13.5 |
| synth_100 | 100 | 30.1 | 62.0 | 7.9 |
| synth_200 | 200 | 22.5 | 72.1 | 5.4 |

**Conclusion:** Indirect scoring dominates at scale (72% at 200 tips). The clip+incremental
phase dominates at small scales because the incremental downpass is O(depth) ≈ O(n) for
small trees (depth ≈ n), while indirect evaluation is O(n²).

### Per-candidate indirect timing

| Dataset | Tips | total_words | Candidates | ns/candidate |
|---------|------|-------------|------------|--------------|
| Vinther2008 | 23 | 28 | 3,585 | 79.9 |
| Agnarsson2004 | 62 | 59 | 56,501 | 91.6 |
| synth_20 | 20 | 11 | 2,535 | 25.6 |
| synth_50 | 50 | 12 | 32,776 | 30.2 |
| synth_100 | 100 | 12 | 237,536 | 33.7 |
| synth_200 | 200 | 12 | 1,090,533 | 32.9 |

**Conclusion:** Per-candidate cost is stable across tree sizes (~33 ns for `total_words=12`),
confirming that cache effects are not increasing per-candidate cost. The cost scales linearly
with `total_words` (28 words → 80 ns, 59 words → 92 ns).

### Scaling analysis

- Indirect time scaling exponent: **2.78** (vs expected 2.0 for O(n²))
- Candidate count scaling exponent: **2.66**
- The super-quadratic scaling is primarily from candidate count growth (2.66),
  not from per-candidate cost degradation (stable at ~33 ns).
- The extra 0.12 exponent may come from TBR rerooting generating O(k) sub-edges
  per clip, where k is subtree size.

### Snapshot overhead

| Tips | Save (μs) | Restore (μs) | Size (KB) |
|------|-----------|---------------|-----------|
| 20 | 0.3 | 0.3 | 14.6 |
| 50 | 1.1 | 1.1 | 40.2 |
| 100 | 2.5 | 2.3 | 80.8 |
| 200 | 5.4 | 5.0 | 162.1 |

**Conclusion:** Snapshot save/restore is negligible — 5 μs per operation at 200 tips,
compared to 36 ms for indirect evaluation. StateSnapshot optimization (Step 6) is not
worth pursuing.

## 2. Steps Investigated and Decisions

### Step 3: Postorder node renumbering — SKIPPED

Analysis of node-ID strides during postorder traversal (50-tip tree):
- Mean stride: 34.6 node IDs (~52 cache lines at `total_words=12`)
- Max stride: 93 node IDs (~140 cache lines)

However, the downpass is **not the hot path** — it's only 22% of time at 200 tips. The
state arrays fit comfortably in L2 (prelim for 200 tips = 37 KB; total state data ≈ 162
KB). Since the bottleneck is indirect scoring (which uses vroot_cache with linear access),
postorder renumbering would not improve the hot path.

**Decision:** Not implemented. Cost/benefit ratio unfavorable.

### Step 4: Binary-character specialization — SKIPPED

Block `n_states` values for typical datasets:
- Vinther2008: 4, 4, 5, 5, 5, 5 (total_words=28)
- Agnarsson2004: 7, 7, 7, 7, 7, 8, 8, 8 (total_words=59)
- synth_200 (binary+NA): 3, 3, 3, 3 (total_words=12)

`n_states` per block is determined by the **total number of applicable states in the
contrast matrix**, not by individual character state coverage. All standard blocks share
the same `n_states`. Binary characters contribute to blocks with the full `n_states`
because `state_remap` assigns globally consecutive indices.

**Decision:** Per-block unrolling for binary characters is not possible with the current
block structure. Changing this would require per-block state counts, which is a deep
architectural change. Not worth it for Phase 3D.

Verified: all inner loops correctly iterate `blk.n_states` (not `total_words`). No bug.

### Step 5: Block-major layout — SKIPPED

The vroot_cache (Phase 2B) already provides linear access for the indirect scoring hot
path. Per-candidate cost is stable across tree sizes, confirming no cache pressure issue.
State arrays for 200 tips fit in L2 (162 KB total).

For morphological data (the target use case), `total_words` is small (12-59) and trees
rarely exceed 500 tips. Block-major layout would add complexity without measurable benefit.

**Decision:** Not implemented. Experiment not justified by profiling data.

### Step 6: StateSnapshot reduction — SKIPPED

Snapshot overhead is <0.01% of TBR pass time at scale. Not worth optimizing.

## 3. Optimizations Applied

### Postorder save/restore in TBR (ts_tbr.cpp)

After `spr_unclip()`, the tree topology is identical to before `spr_clip()`, so the
postorder traversal is the same. Previously, `build_postorder()` (O(n) DFS with vector
allocations) was called to reconstruct it. Now the pre-clip postorder is saved and
restored via `assign()` (O(n) memcpy, no allocation).

Similarly, after `state_snap.restore()` on rejection, the postorder is already restored
by the snapshot's memcpy. The redundant `build_postorder()` calls were removed.

**Changes:**
- Save `tree.postorder` before `spr_clip()`, restore after `spr_unclip()`
- Remove 2 redundant `build_postorder()` calls after `state_snap.restore()`

**Impact:** Eliminates ~377 `build_postorder()` calls per TBR pass at 200 tips. Each call
saves O(n) DFS traversal plus 2 vector allocations. Estimated savings: 1-3% of the
unclip phase. The benefit is modest because unclip is only 5% of total TBR pass time;
the real bottleneck (indirect scoring at 72%) is addressed by Phase 3E (SIMD).

## 4. Implications for Future Phases

### Phase 3E (SIMD) — highest priority

The profiling clearly shows that the **indirect scoring inner loop** is the primary target
for optimization. At 200 tips, it consumes 72% of TBR pass time. The inner loop is:

```cpp
for (int b = 0; b < ds.n_blocks; ++b) {
    uint64_t any_hit = 0;
    for (int s = 0; s < blk.n_states; ++s) {
        any_hit |= (clip_prelim[offset+s] & vroot[offset+s]);
    }
    uint64_t needs_step = ~any_hit & blk.active_mask;
    extra_steps += blk.weight * popcount64(needs_step);
}
```

This is a textbook SIMD target: independent AND/OR operations over contiguous uint64_t
arrays. SSE2 can process 2 words per instruction, AVX2 can process 4. With `n_states`
typically 3-8 per block, even 2× throughput from SSE2 would be significant.

### Algorithmic improvements

The candidate count scaling exponent (2.66 > 2.0) suggests that TBR rerooting generates
more candidates than pure SPR. Reducing the candidate set (e.g., tighter bounds on which
rerootings to try) could reduce the constant factor.

## 5. Files Created/Modified

### Created:
- `inst/benchmarks/bench_memory.R` — profiling harness
- `inst/benchmarks/memory_profile_results.md` — this file
- `tests/testthat/test-ts-memory-layout.R` — 32 regression tests

### Modified:
- `src/ts_rcpp.cpp` — added `ts_bench_tbr_phases` diagnostic (append only), added
  `#include <chrono>` and `#include <random>`
- `src/TreeSearch-init.c` — registered `ts_bench_tbr_phases` (7 args)
- `src/ts_tbr.cpp` — postorder save/restore optimization (3 changes)
- `R/RcppExports.R` — regenerated via `Rcpp::compileAttributes()`
- `src/RcppExports.cpp` — regenerated

### Test status:
- memory-layout: 32/32 passing
- driven: 53/53 passing
- tbr-bench: 26/26 passing
- fuse: 16/16 passing (1 skip)
- sector: 32/32 passing
