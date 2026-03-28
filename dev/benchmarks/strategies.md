# Driven Search Strategy Space

Last updated: 2026-03-17

This document defines all tunable parameters of the C++ driven search
engine (`MaximizeParsimony()`) and proposes named strategy presets for
benchmarking (Phase 6D) and adaptive search (Phase 6F).

## Pipeline Overview

Each replicate executes this fixed phase sequence:

```
Wagner → TBR → XSS → RSS → CSS → Ratchet → Drift → Final TBR
```

Phases may be skipped by setting their cycle/round counts to 0.
Sectorial phases (XSS, RSS, CSS) only run when the tree has
≥ 2 × `sectorMinSize` tips.

Between replicates, the pool collects the best tree(s) and tree
fusing may run (every `fuseInterval` replicates).

---

## Parameter Categories

### A. Strategy Parameters (per-replicate search behavior)

These control how each replicate explores tree space. They are the
primary targets for strategy tuning in Phase 6D.

#### A1. Wagner Start

| R parameter | C++ field | Default | Description |
|-------------|-----------|---------|-------------|
| `wagnerStarts` | `wagner_starts` | 1 | Random Wagner trees built per replicate; best-scoring one used as TBR starting point. Higher values improve starting topology at low cost for small datasets. |

#### A2. TBR

| R parameter | C++ field | Default | Description |
|-------------|-----------|---------|-------------|
| `tbrMaxHits` | `tbr_max_hits` | 1 | Equal-score hits before TBR declares convergence. Higher values explore the plateau more thoroughly. |
| `tabuSize` | `tabu_size` | 100 | Tabu list capacity for TBR. Prevents revisiting recently-explored topologies on plateaus. 0 = disabled. |

#### A3. Ratchet

| R parameter | C++ field | Default | Description |
|-------------|-----------|---------|-------------|
| `ratchetCycles` | `ratchet_cycles` | 10 | Perturbation-then-search cycles per replicate. Primary knob for ratchet intensity. 0 = skip ratchet. |
| `ratchetPerturbProb` | `ratchet_perturb_prob` | 0.04 | Per-character probability of perturbation. Higher = more disruptive. |
| `ratchetPerturbMode` | `ratchet_perturb_mode` | 0 | 0 = zero (silence characters), 1 = upweight (double weight), 2 = mixed (zero some, double others). |
| `ratchetPerturbMaxMoves` | `ratchet_perturb_max_moves` | 0 (auto) | Max TBR moves during perturbation phase. 0 = `max(20, min(200, n_tip/8))`. |
| `ratchetAdaptive` | `ratchet_adaptive` | FALSE | Auto-tune `perturbProb` to target a ~30% escape rate. |

#### A4. Drift

| R parameter | C++ field | Default | Description |
|-------------|-----------|---------|-------------|
| `driftCycles` | `drift_cycles` | 6 | Suboptimal-exploration cycles per replicate. 0 = skip drift. |
| `driftAfdLimit` | `drift_afd_limit` | 3 | Max absolute fit difference (steps) for accepting suboptimal moves. |
| `driftRfdLimit` | `drift_rfd_limit` | 0.1 | Max relative fit difference for accepting suboptimal moves. |

#### A5. Sectorial Search

| R parameter | C++ field | Default | Description |
|-------------|-----------|---------|-------------|
| `xssRounds` | `xss_rounds` | 3 | Exclusive Sectorial Search (systematic partition) rounds. 0 = skip XSS. |
| `xssPartitions` | `xss_partitions` | 4 | Number of non-overlapping sectors per XSS round. |
| `rssRounds` | `rss_rounds` | 1 | Random Sectorial Search rounds after XSS. 0 = skip RSS. |
| `cssRounds` | `css_rounds` | 1 | Constrained Sectorial Search (full-tree exact scoring) rounds. 0 = skip CSS. |
| `cssPartitions` | `css_partitions` | 4 | Partitions for CSS. |
| `sectorMinSize` | `sector_min_size` | 6 | Minimum sector clade size (tips). |
| `sectorMaxSize` | `sector_max_size` | 50 | Maximum sector clade size (tips). |

#### A6. Tree Fusing

| R parameter | C++ field | Default | Description |
|-------------|-----------|---------|-------------|
| `fuseInterval` | `fuse_interval` | 3 | Fuse best tree against pool every N replicates. |
| `fuseAcceptEqual` | `fuse_accept_equal` | FALSE | Accept equal-score fusions (increases pool diversity). |

### B. Convergence Parameters (when to stop)

These control total search effort across replicates. Independent of
per-replicate strategy — benchmarking should generally fix these.

| R parameter | C++ field | Default | Description |
|-------------|-----------|---------|-------------|
| `maxReplicates` | `max_replicates` | 100 | Hard cap on replicates. |
| `targetHits` | `target_hits` | `max(10, n_tip/5)` | Stop after this many independent hits to the best score. |

### C. Pool Parameters

| R parameter | C++ field | Default | Description |
|-------------|-----------|---------|-------------|
| `poolMaxSize` | `pool_max_size` | 100 | Maximum trees retained in the pool. |
| `poolSuboptimal` | `pool_suboptimal` | 0.0 | Score tolerance for retaining suboptimal trees. |

### D. Infrastructure Parameters (not strategy-relevant)

| R parameter | C++ field | Default | Description |
|-------------|-----------|---------|-------------|
| `concavity` | — | Inf | Scoring mode: Inf = EW, finite = IW, "profile" = profile parsimony. |
| `nThreads` | — | 1 | Worker threads. |
| `verbosity` | `verbosity` | 1 | 0 = silent, 1 = per-replicate, 2 = per-phase. |
| `progressCallback` | — | NULL (auto) | Custom progress reporting function. |
| `constraint` | — | (none) | Topology constraint (splits). |
| — | `max_seconds` | 0 | Timeout in seconds (available in C++ bridge, not exposed in R-level `MaximizeParsimony`). |

### E. Not Yet Implemented (noted in production plan)

| Parameter | Description | Status |
|-----------|-------------|--------|
| SPR vs TBR phase choice | Use SPR first, escalate to TBR only where SPR plateaus | Not implemented (T-012) |
| NNI pre-pass | Quick NNI before TBR | Not implemented |

---

## Strategy Vector

For Phase 6D benchmarking, the **strategy vector** consists of the 20
Category A parameters. Each preset specifies values for all 20.

---

## Named Strategy Presets

### 1. `sprint`

Minimal effort for fast interactive exploration. Skips expensive phases.
Suitable as a quick-look default or for very small datasets where a
single TBR pass is often sufficient.

```
wagnerStarts       = 1
tbrMaxHits         = 1
tabuSize           = 0
ratchetCycles      = 3
ratchetPerturbProb = 0.04
ratchetPerturbMode = 0
ratchetPerturbMaxMoves = 0
ratchetAdaptive    = FALSE
driftCycles        = 0        # skip drift
driftAfdLimit      = 3
driftRfdLimit      = 0.1
xssRounds          = 1
xssPartitions      = 4
rssRounds          = 0        # skip RSS
cssRounds          = 0        # skip CSS
cssPartitions      = 4
sectorMinSize      = 6
sectorMaxSize      = 50
fuseInterval       = 5
fuseAcceptEqual    = FALSE
```

**Rationale**: 3 ratchet cycles (vs 10) provides some escape from local
optima without large time cost. No drift (most expensive phase per cycle).
Minimal sectorial (1 XSS round, no RSS/CSS). No tabu (saves memory and
TBR overhead for quick passes).

### 2. `default`

Current production defaults. Balanced for general use.

```
wagnerStarts       = 1
tbrMaxHits         = 1
tabuSize           = 100
ratchetCycles      = 5
ratchetPerturbProb = 0.04
ratchetPerturbMode = 0
ratchetPerturbMaxMoves = 0
ratchetAdaptive    = FALSE
driftCycles        = 2
driftAfdLimit      = 3
driftRfdLimit      = 0.1
xssRounds          = 3
xssPartitions      = 4
rssRounds          = 1
cssRounds          = 0
cssPartitions      = 4
sectorMinSize      = 6
sectorMaxSize      = 50
fuseInterval       = 3
fuseAcceptEqual    = FALSE
```

### 3. `thorough`

More exhaustive exploration. More cycles of everything, adaptive ratchet,
multiple Wagner starts, wider plateau exploration.

```
wagnerStarts       = 3
tbrMaxHits         = 3
tabuSize           = 200
ratchetCycles      = 20
ratchetPerturbProb = 0.04
ratchetPerturbMode = 2        # mixed
ratchetPerturbMaxMoves = 0
ratchetAdaptive    = TRUE
driftCycles        = 12
driftAfdLimit      = 5
driftRfdLimit      = 0.15
xssRounds          = 5
xssPartitions      = 6
rssRounds          = 3
cssRounds          = 2
cssPartitions      = 6
sectorMinSize      = 6
sectorMaxSize      = 80
fuseInterval       = 2
fuseAcceptEqual    = TRUE
```

**Rationale**: Doubles most cycle counts. Adaptive ratchet tunes perturbation
intensity automatically. Mixed perturbation mode (zero + upweight) provides
more diverse perturbation landscapes. More Wagner starts improve starting
point quality. Higher `tbrMaxHits` + `tabuSize` explore plateaus better.
`fuseAcceptEqual` increases pool diversity for fusing.

### 4. `ratchet_heavy`

Emphasize ratchet perturbation for escaping deep local optima. Useful
when the fitness landscape has many local optima separated by large
barriers (common in large datasets with many inapplicable characters).

```
wagnerStarts       = 1
tbrMaxHits         = 1
tabuSize           = 100
ratchetCycles      = 30
ratchetPerturbProb = 0.08
ratchetPerturbMode = 2        # mixed
ratchetPerturbMaxMoves = 0
ratchetAdaptive    = TRUE
driftCycles        = 2        # reduced
driftAfdLimit      = 3
driftRfdLimit      = 0.1
xssRounds          = 1        # reduced
xssPartitions      = 4
rssRounds          = 0        # skip
cssRounds          = 0        # skip
cssPartitions      = 4
sectorMinSize      = 6
sectorMaxSize      = 50
fuseInterval       = 3
fuseAcceptEqual    = FALSE
```

**Rationale**: 3× ratchet cycles, 2× perturbation probability, adaptive
tuning + mixed mode. Drift and sectorial reduced to leave time budget
for ratchet. Most time goes to perturbation-escape cycles.

### 5. `sectorial_heavy`

Emphasize sectorial search for large trees where full-tree TBR is
expensive. Decompose the problem into cheaper subproblems.

```
wagnerStarts       = 1
tbrMaxHits         = 1
tabuSize           = 100
ratchetCycles      = 5        # reduced
ratchetPerturbProb = 0.04
ratchetPerturbMode = 0
ratchetPerturbMaxMoves = 0
ratchetAdaptive    = FALSE
driftCycles        = 3        # reduced
driftAfdLimit      = 3
driftRfdLimit      = 0.1
xssRounds          = 8        # increased
xssPartitions      = 6        # more partitions
rssRounds          = 4        # increased
cssRounds          = 3        # increased
cssPartitions      = 6
sectorMinSize      = 6
sectorMaxSize      = 80       # larger sectors
fuseInterval       = 2
fuseAcceptEqual    = TRUE
```

**Rationale**: Heavy sectorial search (XSS + RSS + CSS) with more
partitions and larger max sector size. Ratchet and drift reduced.
For large trees (60+ tips), sectorial search per-step cost is lower
than full-tree TBR, so more sectorial rounds may yield better
time-to-optimal.

### 6. `drift_heavy`

Emphasize tree drifting for exploring the near-optimal landscape.
Useful when the fitness landscape has broad plateaus or many
near-optimal trees.

```
wagnerStarts       = 1
tbrMaxHits         = 1
tabuSize           = 100
ratchetCycles      = 5        # reduced
ratchetPerturbProb = 0.04
ratchetPerturbMode = 0
ratchetPerturbMaxMoves = 0
ratchetAdaptive    = FALSE
driftCycles        = 20       # increased
driftAfdLimit      = 5        # wider
driftRfdLimit      = 0.2      # wider
xssRounds          = 2        # reduced
xssPartitions      = 4
rssRounds          = 1
cssRounds          = 0        # skip
cssPartitions      = 4
sectorMinSize      = 6
sectorMaxSize      = 50
fuseInterval       = 3
fuseAcceptEqual    = TRUE
```

**Rationale**: 3× drift cycles with relaxed acceptance criteria
(AFD 5, RFD 0.2) allow the search to wander farther from local
optima via incremental suboptimal moves. Ratchet and sectorial
reduced. `fuseAcceptEqual` helps propagate diverse drifted topologies.

---

## Preset Summary Table

| Preset | Wagner | TBR hits | Ratchet | Drift | XSS | RSS | CSS | Fuse int |
|--------|--------|----------|---------|-------|-----|-----|-----|----------|
| sprint | 1 | 1 | 3 cyc | off | 1 rnd | off | off | 5 |
| default | 1 | 1 | 10 cyc | 6 cyc | 3 rnd | 1 rnd | 1 rnd | 3 |
| thorough | 3 | 3 | 20 cyc adaptive | 12 cyc | 5 rnd | 3 rnd | 2 rnd | 2 |
| ratchet_heavy | 1 | 1 | 30 cyc adaptive | 2 cyc | 1 rnd | off | off | 3 |
| sectorial_heavy | 1 | 1 | 5 cyc | 3 cyc | 8 rnd | 4 rnd | 3 rnd | 2 |
| drift_heavy | 1 | 1 | 5 cyc | 20 cyc | 2 rnd | 1 rnd | off | 3 |

---

## Usage in Benchmarking (Phase 6D)

The benchmarking framework should:

1. Fix convergence parameters (`maxReplicates`, `targetHits`) identically
   across presets to make wall-clock comparisons fair.
2. For each benchmark dataset × preset combination, measure:
   - Time to find the best-known score (from `datasets.md`)
   - Total time for convergence or timeout
   - Number of replicates to convergence
   - Phase-level timing breakdown (from `timings` attribute)
3. The results matrix (datasets × presets → metrics) feeds Phase 6E
   (predictive model) and Phase 6F (adaptive search).

## Usage in Adaptive Search (Phase 6F)

The warmup-then-switch approach:
1. Run 2–3 replicates with `default` preset while collecting phase timings.
2. Compute dataset features + phase yield metrics (e.g., "ratchet improved
   score in 80% of cycles" → ratchet-heavy might help).
3. Select the best preset for remaining replicates.

Alternatively, online adaptation could smoothly interpolate between presets
based on per-phase improvement rates.

---

## R Helper Function

The `dev/benchmarks/bench_datasets.R` benchmark utility can use a
`get_strategy(name)` helper. Example:

```r
get_strategy <- function(name = c("sprint", "default", "thorough",
                                   "ratchet_heavy", "sectorial_heavy",
                                   "drift_heavy")) {
  name <- match.arg(name)
  strategies <- list(
    sprint = list(
      wagnerStarts = 1L, tbrMaxHits = 1L, tabuSize = 0L,
      ratchetCycles = 3L, ratchetPerturbProb = 0.04,
      ratchetPerturbMode = 0L, ratchetPerturbMaxMoves = 0L,
      ratchetAdaptive = FALSE,
      driftCycles = 0L, driftAfdLimit = 3L, driftRfdLimit = 0.1,
      xssRounds = 1L, xssPartitions = 4L, rssRounds = 0L,
      cssRounds = 0L, cssPartitions = 4L,
      sectorMinSize = 6L, sectorMaxSize = 50L,
      fuseInterval = 5L, fuseAcceptEqual = FALSE
    ),
    default = list(
      wagnerStarts = 1L, tbrMaxHits = 1L, tabuSize = 100L,
      ratchetCycles = 5L, ratchetPerturbProb = 0.04,
      ratchetPerturbMode = 0L, ratchetPerturbMaxMoves = 0L,
      ratchetAdaptive = FALSE,
      driftCycles = 2L, driftAfdLimit = 3L, driftRfdLimit = 0.1,
      xssRounds = 3L, xssPartitions = 4L, rssRounds = 1L,
      cssRounds = 0L, cssPartitions = 4L,
      sectorMinSize = 6L, sectorMaxSize = 50L,
      fuseInterval = 3L, fuseAcceptEqual = FALSE
    ),
    thorough = list(
      wagnerStarts = 3L, tbrMaxHits = 3L, tabuSize = 200L,
      ratchetCycles = 20L, ratchetPerturbProb = 0.04,
      ratchetPerturbMode = 2L, ratchetPerturbMaxMoves = 0L,
      ratchetAdaptive = TRUE,
      driftCycles = 12L, driftAfdLimit = 5L, driftRfdLimit = 0.15,
      xssRounds = 5L, xssPartitions = 6L, rssRounds = 3L,
      cssRounds = 2L, cssPartitions = 6L,
      sectorMinSize = 6L, sectorMaxSize = 80L,
      fuseInterval = 2L, fuseAcceptEqual = TRUE
    ),
    ratchet_heavy = list(
      wagnerStarts = 1L, tbrMaxHits = 1L, tabuSize = 100L,
      ratchetCycles = 30L, ratchetPerturbProb = 0.08,
      ratchetPerturbMode = 2L, ratchetPerturbMaxMoves = 0L,
      ratchetAdaptive = TRUE,
      driftCycles = 2L, driftAfdLimit = 3L, driftRfdLimit = 0.1,
      xssRounds = 1L, xssPartitions = 4L, rssRounds = 0L,
      cssRounds = 0L, cssPartitions = 4L,
      sectorMinSize = 6L, sectorMaxSize = 50L,
      fuseInterval = 3L, fuseAcceptEqual = FALSE
    ),
    sectorial_heavy = list(
      wagnerStarts = 1L, tbrMaxHits = 1L, tabuSize = 100L,
      ratchetCycles = 5L, ratchetPerturbProb = 0.04,
      ratchetPerturbMode = 0L, ratchetPerturbMaxMoves = 0L,
      ratchetAdaptive = FALSE,
      driftCycles = 3L, driftAfdLimit = 3L, driftRfdLimit = 0.1,
      xssRounds = 8L, xssPartitions = 6L, rssRounds = 4L,
      cssRounds = 3L, cssPartitions = 6L,
      sectorMinSize = 6L, sectorMaxSize = 80L,
      fuseInterval = 2L, fuseAcceptEqual = TRUE
    ),
    drift_heavy = list(
      wagnerStarts = 1L, tbrMaxHits = 1L, tabuSize = 100L,
      ratchetCycles = 5L, ratchetPerturbProb = 0.04,
      ratchetPerturbMode = 0L, ratchetPerturbMaxMoves = 0L,
      ratchetAdaptive = FALSE,
      driftCycles = 20L, driftAfdLimit = 5L, driftRfdLimit = 0.2,
      xssRounds = 2L, xssPartitions = 4L, rssRounds = 1L,
      cssRounds = 0L, cssPartitions = 4L,
      sectorMinSize = 6L, sectorMaxSize = 50L,
      fuseInterval = 3L, fuseAcceptEqual = TRUE
    )
  )
  strategies[[name]]
}
```

This helper will be formalized in the benchmarking framework (T-004).

---

## External Benchmark Datasets (MorphoBank corpus)

### Train/validation split

The `neotrans/inst/matrices/` directory contains ~800 MorphoBank phylogenetic
matrices. These supplement the 14 bundled datasets for broader, less
overfitting-prone benchmarking.

**Split rule:** A matrix belongs to the **validation** set if its MorphoBank
project number is divisible by 5 (i.e., `project_id %% 5 == 0`); all others
are **training**. The 7 `syab*` files (non-MorphoBank) are always training.

After filtering (ntax ≥ 20, parse OK, dedup): 535 training, 124 validation.

**Usage rules:**
- **Training** matrices may be used freely during development and tuning.
- **Validation** matrices are a **one-way door**: run once to confirm that
  improvements generalize. Results must **never** inform strategy tuning.
- If validation is ever used for tuning, the split is compromised and must
  be rebuilt with a new rule.

### Dedup

Multi-file projects (same MorphoBank project, separate `.nex` files) often
contain the same character matrix with minor taxon-sampling variations. These
are flagged as `dedup_drop = TRUE` in the catalogue. The dedup uses pairwise
character identity ≥ 95% on shared taxa (requiring ≥ 80% taxon overlap),
keeping the largest matrix per redundancy cluster.

24 near-duplicates are excluded, leaving 659 usable matrices.

### Fixed 25-matrix training sample

For routine benchmarking, a fixed sample of 25 matrices is used
(`MBANK_FIXED_SAMPLE` in `bench_datasets.R`). Selected via max-min distance
on standardized (ntax, nchar, pct_missing, pct_inapp) within each tier:

| Tier | Count | Keys |
|------|-------|------|
| Small (20–30) | 7 | project532, project2346, project2451, project4501, project944, project971_(1), project2762 |
| Medium (31–60) | 7 | project826, project561, project571, project4146_(3), project3688, project4049, project423 |
| Large (61–120) | 7 | project4286, project4359, project4397, project2084_(1), project2771, project2184, project3938 |
| XLarge (121+) | 4 | syab07201, project4133, project804, project4284 |

**Do not modify this list.** Benchmark comparisons require the same sample.

### Fixed 20-matrix Brazeau-track sample

For benchmarking under Brazeau et al. (2019) NA-aware scoring, a separate
fixed sample is used (`MBANK_BRAZEAU_SAMPLE` in `bench_datasets.R`).
Selected from training matrices with **pct_inapp ≥ 4%** (where the
three-pass algorithm materially differs from Fitch scoring). Max-min
distance selection on (ntax, nchar, pct_inapp):

| Tier | Count | Keys |
|------|-------|------|
| Small (20–30) | 5 | project4182, project2346, project906, project4112, project537 |
| Medium (31–60) | 6 | project709, project561, project2359, project4761, project4146_(3), project4867 |
| Large (61–120) | 6 | project4286, project3512_(2), project2084_(1), project2086, project2771, project3938 |
| XLarge (121+) | 3 | project4103, project804, syab07204 |

**Do not modify this list.** Benchmark comparisons require the same sample.

---

## Benchmark Tracks

Strategy tuning and performance comparison use two distinct tracks,
each run under both equal weights (EW) and implied weights (IW, k=10).

### Track 1: Fitch (TNT comparison and core search quality)

- **Scoring:** `fitch_mode()` (inapplicable treated as missing)
- **Datasets:** 14 bundled + `MBANK_FIXED_SAMPLE` (25 matrices, all scorings)
- **Purpose:** Direct comparison with TNT; measures raw search quality on
  the same objective function
- **TNT comparison:** `bench_tnt_compare.R` with `use_fitch = TRUE`

### Track 2: Brazeau (NA-algorithm-specific strategy tuning)

- **Scoring:** Default Brazeau et al. (2019) three-pass algorithm
- **Datasets:** `MBANK_BRAZEAU_SAMPLE` (20 matrices, pct_inapp ≥ 4%)
  + relevant bundled datasets
- **Purpose:** Tests whether the different per-evaluation cost and fitness
  landscape of NA-aware scoring warrant different strategy presets
- **No TNT comparison:** TNT does not implement Brazeau scoring

### Weighting dimension

Both tracks should be run under:
- **EW** (equal weights, `concavity = Inf`)
- **IW k=10** (implied weights, `concavity = 10`)

IW reshapes the fitness landscape (homoplasy penalties change which local
optima exist), so the optimal search strategy may differ between EW and IW.
k=10 is chosen as representative: moderate weighting, commonly used in
practice.

### Priority order for new benchmark runs

1. **Brazeau × EW** on `MBANK_BRAZEAU_SAMPLE` — biggest current gap
2. **Fitch × IW(k=10)** on `MBANK_FIXED_SAMPLE` — detects IW strategy needs
3. **Brazeau × IW(k=10)** on `MBANK_BRAZEAU_SAMPLE` — full interaction
   (only if 1 or 2 reveals meaningful differences)

### Sample-size validation protocol

After the first Brazeau × EW benchmark round, run a strategy-ranking
stability analysis to confirm 20 matrices are sufficient:

1. For k = 5, 8, 10, 12, 14, 16, 18, 20: bootstrap 500 subsamples of
   size k from the 20 matrices.
2. Within each subsample, rank strategies by mean score gap (strategy
   score − best-of-strategies, averaged across k datasets).
3. Compute Kendall's W between each bootstrap ranking and the full-20
   ranking.
4. Plot median W vs k. If W has plateaued at k=20 (ΔW < 0.02 for last
   2 increments), 20 is sufficient. If still climbing, expand the sample.

---

## Brazeau-track Phase Profiling (T-290b, 2026-03-28)

**Method:** 6 Brazeau-sample datasets (project2346/23t, project2359/42t,
project4146\_(3)/59t, project2086/91t, project2771/94t, project804/173t)
× 4 scoring conditions ({Brazeau, Fitch} × {EW, IW k=10}) × 2 strategy
presets × 3 seeds. 30s timeout. Timings from `attr(result, "timings")`.

### Brazeau / Fitch cost ratio by phase

| Phase | EW/default | EW/thorough | IW10/default | IW10/thorough |
|-------|:----------:|:-----------:|:------------:|:-------------:|
| **Wagner** | **3.58×** | **3.85×** | **4.99×** | **5.20×** |
| TBR | 0.89× | 0.89× | 0.69× | 0.89× |
| XSS | 1.00× | 1.14× | 0.95× | 0.94× |
| RSS | 1.38× | 1.29× | 1.23× | 1.06× |
| CSS | — | 1.30× | — | 1.08× |
| Ratchet | **1.33×** | 1.29× | 1.23× | 1.09× |
| Final TBR | 1.23× | 1.19× | 1.02× | 0.96× |

**Key finding:** Wagner is the extreme outlier at 3.6–5.2× under Brazeau.
All other phases are within 0.7–1.4× of Fitch cost. The IW penalty for
Wagner is larger than EW (4.99–5.20× vs 3.58–3.85×), likely because IW
weighting changes which states are ambiguous during NA-aware scoring.

### Phase distribution (% of wall time)

| Phase | Brazeau/EW/default | Brazeau/EW/thorough | Fitch/EW/default | Fitch/EW/thorough |
|-------|:------------------:|:-------------------:|:----------------:|:-----------------:|
| Ratchet | **74.3%** | **63.0%** | **75.6%** | **65.0%** |
| Wagner | 8.9% | 6.7% | 4.5% | 2.8% |
| TBR | 6.7% | 4.3% | 8.3% | 5.2% |
| XSS | 5.3% | 6.4% | 6.4% | 6.8% |
| RSS | 2.9% | 9.6% | 3.1% | 9.7% |
| CSS | — | 7.2% | — | 7.0% |
| Final TBR | 1.9% | 2.1% | 2.1% | 2.3% |

Under Brazeau, Wagner is nearly double its Fitch share (8.9% vs 4.5%
for default; 6.7% vs 2.8% for thorough). All other phases are very close
between scoring modes. The landscape structure (ratchet dominance, sectorial
share) is essentially the same under Brazeau as under Fitch.

### Replicate rates

| Condition | Mean reps/30s |
|-----------|:------------:|
| Brazeau/EW/default | 26.8 |
| Brazeau/EW/thorough | 23.3 |
| Brazeau/IW10/default | 23.2 |
| Brazeau/IW10/thorough | 19.3 |
| Fitch/EW/default | 27.7 |
| Fitch/EW/thorough | 24.8 |
| Fitch/IW10/default | 24.3 |
| Fitch/IW10/thorough | 22.0 |

Brazeau runs at 95–97% of Fitch rep rate for the same preset. The landscape
is not meaningfully harder to navigate under Brazeau.

### Default vs thorough score comparison (Brazeau, 30s, 20-matrix sample)

| Weighting | Datasets where thorough better | Mean best delta | Notes |
|-----------|:---:|:---:|---|
| EW | 2/20 (10%) | +12.6 | Driven by 86t (+113) and 225t (+139) datasets |
| IW10 | 2/20 (10%) | +0.03 | Very small; 3/20 datasets show near-zero gap |

The median improvement is 0 in all conditions. Thorough helps only at
≥86 tips (where the extra ratchet cycles and sectorial rounds have room
to operate). For ≤72 tips, default and thorough produce identical results
at 30s.

### Implications for Brazeau preset tuning

1. **Wagner cost:** `wagnerStarts = 3` (thorough) costs 3.85–5.20× more
   per start under Brazeau vs Fitch, but Wagner is only 6–9% of wall time.
   Reducing to `wagnerStarts = 1` saves ~4% of replicate time, enabling
   ~1 extra replicate per 30s run. Whether that extra replicate outweighs
   better starting topology depends on dataset size. See T-290c for
   wagnerStarts = 1 vs 3 experiment.

2. **Ratchet:** The 1.1–1.3× Brazeau overhead per ratchet cycle is modest,
   and the landscape difficulty is unchanged (same rep rate, same ratchet
   dominance). No evidence that Brazeau requires more ratchet cycles than
   Fitch. Current settings (12 cycles/default, 20/thorough) are appropriate.

3. **Overall:** Fitch-tuned presets are appropriate for Brazeau. The only
   parameter worth revisiting is `wagnerStarts` in the thorough preset,
   and only for larger datasets where Wagner per-start cost is non-trivial.

### wagnerStarts = 1 vs 3 under Brazeau (analytical, T-290)

The default vs thorough score comparison reveals the mechanism behind the
wagnerStarts benefit:

| Dataset | n_tax | n_char | med_reps (default) | med_reps (thorough) | best_delta (EW) |
|---------|:-----:|:------:|:------------------:|:-------------------:|:---------------:|
| project2084\_(1) | 86 | 3660 | **0** | **0** | +113 (thorough better) |
| syab07204 | 225 | ~200 | **0** | **0** | +139 (thorough better) |
| project2086 | 91 | 453 | 7 | 4 | −1 (default slightly better) |
| project2771 | 94 | ~450 | 9 | 6 | 0 (tied) |
| ≤72 tips (all) | ≤72 | — | ≥11 | ≥10 | 0 (all tied) |

The two datasets where thorough is significantly better (86t/3660c and
225t) **complete zero replicates in 30s** for both presets. The score
difference is therefore NOT from extra replicates or more ratchet cycles —
it is entirely from **better starting topology quality**.

With wagnerStarts=3 + nniFirst=TRUE, the search evaluates 3 Wagner trees,
runs NNI from each, and selects the best NNI-local optimum as the TBR
starting point. At 86–225 tips with hundreds of characters, TBR convergence
takes 30+ seconds per run. In that regime:
- The extra Wagner+NNI setup cost (3.6–5.2× per start, but still <500ms
  total for 3 starts at 86t) is negligible relative to TBR convergence time.
- The better starting point found by 3-start selection meaningfully reduces
  the remaining distance to the TBR local optimum.

For datasets where multiple replicates complete (≤94 tips, ≤450 chars),
wagnerStarts=1 and wagnerStarts=3 produce equivalent scores. The extra
Wagner starts provide no marginal benefit once ratchet escape drives the
search across multiple replicates.

**Conclusion:** `wagnerStarts = 3` in the thorough preset is the correct
choice under Brazeau, for exactly the same reason as under Fitch: it
dominates at large/complex datasets where initial descent quality is the
binding constraint, and is negligible at small/medium datasets where
ratchet escape drives quality. The 3.6–5.2× per-start Wagner cost under
Brazeau is irrelevant in practice — Wagner overhead is ~7% of wall time,
and TBR convergence is orders of magnitude more expensive at the scales
where wagnerStarts matters.

### T-290c confirmatory: wagnerStarts = 1 vs 3, empirical (2026-03-28)

Quick benchmark (EW, 5 seeds, 30s + 60s) on two contrasting datasets to
validate the T-290b analytical conclusion. Both under Brazeau scoring.

| Dataset | n_tax | n_char | timeout | med_reps | w1_best | w3_best | delta_best |
|---------|:-----:|:------:|:-------:|:--------:|:-------:|:-------:|:----------:|
| project2084\_(1) | 86 | 3660 | 30s | 0 | **29258** | 29385 | **−127 (w1 better)** |
| project2084\_(1) | 86 | 3660 | 60s | 1 | 29196 | **28632** | **+564 (w3 better)** |
| project2086 | 91 | 453 | 30s | 8 | 2180 | **2133** | +47 (w3 better) |
| project2086 | 91 | 453 | 60s | 17 | 2169 | **2133** | +36 (w3 better) |

**The 30s / 0-rep result is surprising:** wagnerStarts=1 is better for
project2084_(1) at 30s. The reason: at 86t with 3660 chars, Brazeau Wagner
construction is expensive (~4× more than Fitch). Three starts consume
~20-25s of a 30s budget, leaving almost no TBR time. With one start,
25s+ is available for TBR — and TBR progress from a mediocre starting
topology beats no TBR at all.

**At 60s (≥1 rep completed):** wagnerStarts=3 wins by a large margin (+564
steps). The better starting topology matters at large-dataset scale.

**Reconciling with T-290b:** The T-290b reasoning estimated Wagner overhead
as "<500ms total for 3 starts at 86t". This holds for small/medium datasets
(≤500 chars). At 3660 chars with Brazeau, Wagner is far more expensive. The
30s/0-rep case is a **large-preset** scenario, not a thorough-preset scenario.

**Final conclusion (T-290b + T-290c combined):**

| Regime | Recommended wagnerStarts | Reasoning |
|--------|:------------------------:|-----------|
| budget >> replicate_time (multiple reps) | 3 | Small benefit; ratchet dominates |
| budget ≈ replicate_time (∼1 rep) | 3 | Better topology → large improvement |
| budget < replicate_time (0 reps, large preset) | 1 | Less Wagner → more TBR → better score |

The current preset assignments are correct:
- **thorough** (65–119t, ≤500c): wagnerStarts=3 ✓ (gets multiple reps or ~1 rep)
- **large** (≥120t, high char count): wagnerStarts=1 ✓ (0–1 reps; TBR time > topology quality)
MDEOF 2>&1
