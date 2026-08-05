# Per-candidate kernel baseline on REAL large matrices (near-optimal regime)

**Date:** 2026-07-14 · Harness `dev/profiling/reeval/kern_real.R` (load .nex, fitch-convert `-`→`?`,
converge to a kernel local optimum, then time a near-optimal round via committed `TS_IW_TIMING`,
min-of-runs). Build `.agent-disc` (worktree tip). Matrices pulled from Hamilton neotrans to
`.agent-disc/matrices/`. This is the number the kernel work must beat.

| dataset | nTip | nChar | nStates | SPR ns | REROOT ns | **TOT/cand ns** |
|---|---|---|---|---|---|---|
| project4204 | 163 | 37   | 2 | 8.6  | 5.6  | **6.1** |
| project175  | 165 | 71   | 6 | 9.6  | 6.5  | **7.0** |
| project4327 | 197 | 823  | 6 | 26.4 | 16.5 | **17.6** |
| project2668 | 196 | 1227 | 6 | 44.3 | 40.0 | **40.6** |
| project970  | 157 | 1844 | 6 | 79.7 | 75.3 | **76.0** |

## Key finding — cost ∝ total_words (candidates do NOT bail early here)

Per-candidate cost scales strongly with nChar (37c→6ns … 1844c→76ns; super-linear at the top =
cache). If near-optimal candidates bailed after ~1 block (Goloboff's "first few characters exceed
the bound"), cost would be ~flat in nChar. It is not. **Cause:** on a converged, plateau-rich
morphological tree, most reroot neighbours are EQUAL length; with `acceptEqual=FALSE` + strict-`<`
cutoff, an equal-length candidate never *exceeds* the cutoff, so it processes **all blocks** — a
full `edge_set[below]` gather + reduce over all `total_words = n_blocks × n_states`. REROOT is
77–81% of candidates and dominates (project970 REROOT 75 ns of the 76).

**Implication for levers:**
- A block-0-hot cache re-layout helps LESS than a shallow-bail model predicts — the reroot reads
  *all* blocks, not just block-0. The gather that hurts is the *whole* `edge_set[below]`.
- The high-payoff targets are **(a) reducing words gathered/reduced per candidate** (incremental
  cross-candidate update so successive rerootings don't re-gather the full edge set; or a compacter
  edge_set footprint to cut cache misses) and **(b) the equal-length plateau processing** — an
  equal-length candidate could perhaps be recognised/short-circuited without a full scan.
- ns/word grows with total_words (0.23 → 0.34 → 0.44 ns/word for 823 → 1227 → 1844 chars) =
  cache-miss pressure as the edge_set buffer outgrows cache. Footprint reduction cuts this directly.

The `kernel-reslice-design` workflow was launched with a block-0-bail framing; reconcile its output
against this ∝-total_words / full-gather reality (the equal-length plateau, not shallow bail, is the
governing regime for the dominant reroot path).
