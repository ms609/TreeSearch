# Benchmark Results Analysis (Agent A, T-005)

## Dataset

8 datasets × 6 strategies × 3 reps = 144 planned runs.
55/144 succeeded (38%) due to T-025 optimization-dependent UB segfault.
Aria2015 (35 tips) and Dikow2009 (88 tips) had highest crash rates.

## Key Findings

### 1. All strategies find optimal on small datasets (≤43 tips)
- Longrich (20 tips), Vinther (23 tips), Griswold (43 tips): 100% optimal
- Strategy choice doesn't matter much for small datasets

### 2. Thorough and ratchet_heavy win on large datasets
- Zhu2013 (75 tips): `thorough` found best-known (649), `sprint` failed (652)
- Giles2015 (78 tips): `ratchet_heavy` found best (714), others 716-720
- Dikow2009 (88 tips): `ratchet_heavy` and `drift_heavy` both found 1612 (vs best-known 1614)

### 3. Sprint is fastest but loses quality at scale
- Sprint uses 3 ratchet cycles, no drift, minimal sectorial
- At ≤43 tips: optimal quality, 2-10× faster wall time
- At 75+ tips: fails to find optimal within 20s timeout

### 4. Phase time distribution depends strongly on strategy
| Strategy | TBR | Ratchet | Drift | Sectorial | Fuse |
|----------|-----|---------|-------|-----------|------|
| sprint | 43% | 42% | 0% | 9% | 1% |
| default | 11% | 37% | 39% | 11% | 0% |
| ratchet_heavy | 6% | 87% | 5% | 1% | 0% |
| sectorial_heavy | 13% | 20% | 21% | 38% | 7% |
| drift_heavy | 7% | 12% | 74% | 4% | 3% |

### 5. Replicates-to-convergence varies by strategy
- Sprint: 16-43 reps (many cheap reps)
- Thorough: 6-10 reps (few expensive reps)
- At 20s timeout, sprint completes 35-100 reps; thorough completes 6-10

## Recommendations for Adaptive Strategy

1. **Size-based switching**: Use sprint for ≤30 tips, default for 30-60,
   thorough or ratchet_heavy for 60+.
2. **Phase timing feedback**: If ratchet/drift phases dominate but scores
   aren't improving, switch to more replicates with lighter per-replicate effort.
3. **Time budget**: With short timeouts, sprint covers more replicates.
   With longer timeouts, thorough explores deeper per replicate.

## Limitations

- Only 38% of runs succeeded due to T-025 bug
- 20s timeout limits large-dataset exploration
- No IW or profile parsimony benchmarks (EW only)
