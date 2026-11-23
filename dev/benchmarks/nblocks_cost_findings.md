# Per-Candidate Cost vs Number of Character Blocks

**Task:** T-075  
**Date:** 2026-03-18  
**Agent:** A

## Setup

- 9 neotrans matrices selected from the 100–130 tip range
- All have inapplicable characters (NA-aware scoring)
- 5 random tree seeds per matrix
- Measured via `ts_bench_tbr_phases()` (one full TBR clip–evaluate–unclip pass)

## Key finding

Per-candidate indirect scoring cost is **linear** in both `n_blocks` and
`total_words`, with no significant nonlinearity (quadratic term p = 0.41).

### Model: `ns_per_cand ~ n_blocks + total_words`

| Term | Coefficient | SE | Interpretation |
|------|------------|-----|----------------|
| intercept | 2.4 ns | 0.7 | Base overhead per candidate |
| n_blocks | 3.3 ns | 0.2 | Per-block overhead (loop, function call) |
| total_words | 0.29 ns | 0.02 | Per-word cost (bit-parallel ops) |

R² = 0.990 (45 observations from 9 datasets × 5 seeds)

### Predicted cost at range extremes

| n_blocks | total_words | Predicted ns/candidate | Observed mean |
|----------|-------------|----------------------|---------------|
| 3 | 16 | 17.1 | 17.8 |
| 11 | 76 | 61.2 | 61.2 |

Ratio: 3.6× cost increase from simplest to most complex dataset.

### Standalone models

- `n_blocks` alone: R² = 0.931, slope ≈ 5.4 ns/block
- `total_words` alone: R² = 0.885, slope ≈ 0.62 ns/word

## Practical implications

1. **No threshold effect**: Cost scales linearly — there's no critical
   n_blocks value after which performance degrades sharply.

2. **Block overhead dominates**: At typical total_words (30–80), the per-block
   overhead (3.3 ns × n_blocks) contributes more than per-word cost
   (0.29 ns × total_words) for datasets with many state-count groups.

3. **Optimisation opportunity**: Merging blocks with adjacent state counts
   (e.g., 5-state and 6-state characters into a single padded block) could
   reduce n_blocks by 2–4, saving ~7–13 ns/candidate. At 300k candidates
   per clip, this would save ~2–4 ms per clip pass, or ~100–200 ms across
   a full TBR sweep with 50 clips. Meaningful for large datasets but not
   critical — this is a low-priority micro-optimisation.

4. **For strategy selection**: n_blocks can be computed cheaply at dataset
   load time. Datasets with n_blocks ≥ 10 will have ~3× higher per-candidate
   cost than datasets with n_blocks ≤ 4, which affects expected search
   duration. This could inform time estimates in the Shiny app.

## Data

Raw results: `nblocks_cost_bench.csv` (45 rows: 9 datasets × 5 seeds)

### Datasets used

| File | n_tips | n_char | n_blocks | total_words | Mean ns/cand |
|------|--------|--------|----------|-------------|-------------|
| project2144.nex | 109 | 123 | 3 | 16 | 17.8 |
| project987.nex | 108 | 114 | 4 | 30 | 24.2 |
| project2191.nex | 105 | 215 | 5 | 46 | 32.6 |
| project3422.nex | 110 | 277 | 6 | 42 | 34.0 |
| project4264 (1).nex | 112 | 441 | 7 | 50 | 39.4 |
| project1157.nex | 110 | 138 | 8 | 26 | 36.1 |
| project691.nex | 103 | 443 | 9 | 64 | 53.6 |
| project625.nex | 106 | 236 | 10 | 86 | 59.9 |
| project2292.nex | 114 | 493 | 11 | 76 | 61.2 |
