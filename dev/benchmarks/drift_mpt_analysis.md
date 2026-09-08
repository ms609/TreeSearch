# T-254: Drift MPT Diversity Experiment

## Question

Drift search consumes 15–19% of wall time but contributes <1% of score
improvement (T-251). Before reducing it, we need to check whether drift
helps **MPT enumeration** — finding topologically distinct optimal trees
that the post-search TBR plateau walk uses as seeds.

## Design

- **Datasets**: Wortley2006 (37t), Zhu2013 (75t), Geisler2001 (68t)
- **Conditions**: `driftCycles=0` vs `driftCycles=2` (default preset value)
- **Seeds**: 1, 2, 3
- **Budgets**: 30s (primary, equal-budget), 120s (with consensus stopping)
- **Other params**: All match `default` preset (ratchet 12 cycles, 25%
  perturbation, XSS 3 rounds, etc.)
- **Metrics**: best score, pool tree count, n_topologies, replicates
  completed, mean pairwise Robinson-Foulds distance

### Equal-budget design

The primary comparison uses `consensusStableReps=0` to disable
consensus-stability early stopping. This ensures both conditions use the
full 30s budget, avoiding the confound that no-drift converges to consensus
stability faster (fewer replicates needed to stabilize the strict consensus).

## Results (30s, equal budget)

| Dataset     | Drift | Med score | Med trees | Med reps | Med RF | Drift % |
|-------------|:-----:|:---------:|:---------:|:--------:|:------:|:-------:|
| Geisler2001 |   0   |   1295    |    100    |    27    |  7.3   |    0    |
| Geisler2001 |   2   |   1295    |    100    |    25    |  7.4   |   18    |
| Wortley2006 |   0   |    482    |     4     |    75    |  17.3  |    0    |
| Wortley2006 |   2   |    482    |     2     |    62    |  10.0  |   15    |
| Zhu2013     |   0   |    638    |    100    |    26    |  11.6  |    0    |
| Zhu2013     |   2   |    638    |    100    |    19    |  10.2  |   17    |

### Replicate cost

| Dataset     | Reps (d=0) | Reps (d=2) | Loss |
|-------------|:----------:|:----------:|:----:|
| Geisler2001 |     27     |     24     | 10%  |
| Wortley2006 |     76     |     61     | 20%  |
| Zhu2013     |     25     |     20     | 22%  |

### Key findings

1. **Score quality**: Identical. Both conditions find the same best score
   on all datasets at all seeds.

2. **MPT count**: On Wortley2006, no-drift consistently finds 4 MPTs
   (all 3 seeds) while drift finds 1–3 (median 2). On larger datasets,
   both fill the 100-tree pool. Drift does NOT help MPT enumeration.

3. **Topological diversity**: Mean pairwise RF distances are essentially
   identical on Geisler2001 (7.3 vs 7.4 out of max 132). On Zhu2013,
   no-drift shows slightly higher RF (11.6 vs 10.2 out of max 146).
   On Wortley2006, no-drift has higher RF (17.3 vs 10.0 out of max 70).
   **Drift does not improve topological diversity.**

4. **Replicate throughput**: No-drift completes 10–22% more replicates
   in the same wall time. Each independent replicate starts from a random
   Wagner tree, providing more diverse initial basins than drift's local
   perturbation within a single basin.

5. **Consensus stability confound**: With consensus stopping enabled
   (120s budget), no-drift reaches consensus stability 2–3× faster and
   stops early. Drift prevents early stabilization (by perturbing into
   slightly different topologies) but the extra time produces no better
   scores or more MPTs. This means drift actively delays convergence
   without adding value.

## Conclusion

**Drift can be safely eliminated from the default preset.** It provides:
- Zero score benefit (confirmed both here and in T-251)
- Zero MPT enumeration benefit (fewer MPTs on Wortley2006)
- Zero topological diversity benefit
- Negative throughput impact (10–22% fewer replicates)

The time saved should be reallocated to additional replicates (which
provide genuinely independent basin sampling via random Wagner starts).

## Recommendation for T-255

- **default**: `driftCycles = 0` (was 2)
- **sprint**: already 0 (no change)
- **thorough**: reduce from 12 to 0 or 1. The thorough preset has many
  other escape mechanisms (NNI-perturbation, adaptive ratchet, outer
  cycles) that make drift redundant.
- **large**: already 0 (no change)

## Scripts and data

- `dev/benchmarks/bench_drift_mpt.R` — full experiment script
- `dev/benchmarks/results_drift_mpt_30s.csv` — 30s with consensus stopping
- `dev/benchmarks/results_drift_mpt_120s.csv` — 120s with consensus stopping
- `dev/benchmarks/results_drift_mpt_30s_nostop.csv` — 30s equal-budget (primary)
