# T-069 Stress Test Findings — 150–225 taxa
Agent F, 2026-03-18

## Datasets

| File | Taxa | Chars | NA blocks | Inapplicable |
|------|------|-------|-----------|--------------|
| project175.nex  | 165 | 71  | 2  | 0%    |
| project3763.nex | 205 | 103 | 3  | 50.1% |
| syab07204.nex   | 225 | 748 | 12 | 25.1% |

## Key Findings

### 1. Scaling exponents (synthetic series, n=20–225)

| Metric | Exponent | Expected |
|--------|----------|---------|
| `n_candidates` | **n^2.86** | O(n^2) = 2.0 |
| `indirect_us`  | **n^2.73** | — |
| `clip_incr_us` | **n^1.50** | — |

Candidate count scales slightly super-quadratically (larger pruned subtrees give more valid regraft positions). Indirect scoring tracks candidates closely. Clip/incremental is sub-linear relative to candidates — incremental state amortises well.

Both exponents are consistent with the existing AGENTS.md note (~n^2.8 TBR cost).

### 2. NA block count drives per-candidate cost

| Dataset | n_tips | n_blocks | ns/candidate |
|---------|--------|----------|--------------|
| project175  | 165 | 2  | 12.6 ns |
| project3763 | 205 | 3  | 19.2 ns |
| syab07204   | 225 | 12 | **57.5 ns** |

syab07204's 12 NA character blocks cause ~4.6× higher per-candidate cost than the 2-block case, and 3× higher than 3-block. The NA three-pass scoring cost is proportional to n_blocks, not just n_tips. This is a real bottleneck for large, character-rich matrices with many inapplicable characters.

The existing baseline in AGENTS.md (`~23 ns at 75 tips`) was measured on small inapplicable.phyData sets. Large real matrices with many NA blocks can be 2–3× slower per candidate.

### 3. TBR fraction surpasses ratchet+drift at 200+ taxa

| Dataset | TBR% | Ratchet% | Drift% |
|---------|------|----------|--------|
| project175 (165t, thorough)  | 17% | 38% | 42% |
| project3763 (205t, default)  | **57%** | 13% | 28% |
| syab07204 (225t, default)    | **49%** | 13% | 27% |

At ≤100 taxa, ratchet+drift dominate (~65–70%). At 200+ taxa, TBR itself becomes the largest single cost (49–57%). This crossover happens around 150–175 taxa. The phase distribution shift is driven by the super-quadratic TBR cost overwhelming the approximately-linear perturbation overhead.

### 4. Pool collapse at large n with many characters

syab07204 (225t, 748 chars) produced pool sizes of **8 and 2** from 2 replicates (2 reps each, nThreads=2). In contrast, project3763 (205t, 103 chars) filled the 100-tree pool even from 2 reps.

The near-empty pool for syab07204 means:
- Tree fusing has almost no material to work with
- MPT enumeration from the pool will be from very few seeds
- Users may get poor solutions without many more replicates

This is expected behaviour (each TBR pass takes ~150ms, so a 2-rep run completes very few TBR iterations), but it highlights that **recommended replicates should scale with taxa × chars**. At 225t / 748 chars, users need 10–20+ replicates for reliable results.

### 5. Score variability at large n

| Dataset | Score seed1 | Score seed2 | Δ |
|---------|------------|------------|---|
| project175  | 419  | 424  | 5  (1.2%) |
| project3763 | 1643 | 1513 | 130 (7.9%) |
| syab07204   | 11785 | 11933 | 148 (1.3%) |

project3763 shows high variability (7.9%) despite only 205 taxa — likely because the 50% inapplicable data creates a very complex landscape. High inapplicable fractions interact with the NA three-pass to create many near-equal plateau trees.

### 6. Memory (snapshot bytes per TBR pass)

| Dataset | Snapshot KB |
|---------|------------|
| project175 (165t, 2 blocks)  | 66.8 KB |
| project3763 (205t, 3 blocks) | 290.8 KB |
| syab07204 (225t, 12 blocks)  | **547.2 KB** |

Snapshot memory is manageable (well under 1 MB per pass), but the 547 KB for syab07204 means that with nThreads=2 each thread carries ~1 MB of snapshot state. Not a memory problem, but cache pressure contributes to the elevated per-candidate cost.

## Suggested Follow-up Tasks

- **T-073 (potential)**: Benchmark per-candidate cost as a function of `n_blocks` (hold n_tips fixed). Determine whether there's a block-count threshold beyond which a different NA scoring strategy would help.
- **T-074 (potential)**: Auto-scale `maxReplicates` recommendation in `SearchControl()` based on n_tips × n_chars × n_blocks.
- Revisit `thorough` strategy for large char-dense matrices: at 225t/748 chars, the ratchet+drift overhead is proportionally small (40%), so increasing ratchet/drift cycles is cheap relative to per-pass TBR cost.
