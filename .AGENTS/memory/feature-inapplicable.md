# Alternative Inapplicable-Handling Algorithms

Load this when: working on HSJ scoring, x-transform recoding, Sankoff engine,
`inapplicable=` parameter, or character hierarchy specification.

Plan: `.positai/plans/2026-03-19-0643-alternative-inapplicable-handling-algorithms.md`

Adding HSJ (Hopkins & St. John 2021) and step-matrix/x-transformation
(Goloboff et al. 2021) scoring as alternatives to the existing Brazeau
et al. (2019) three-pass algorithm. Both require an explicit character
hierarchy specification.

---

## New files

| File | Purpose | Status |
|------|---------|--------|
| `R/CharacterHierarchy.R` | `CharacterHierarchy` S3 class, `validate_hierarchy()`, `hierarchy_from_names()`, `hierarchy_chars()`, `hierarchy_controlling()`, `non_hierarchy_weights()` | Complete, 34 tests passing |
| `tests/testthat/test-CharacterHierarchy.R` | Unit tests for hierarchy specification + weight partitioning | Complete |
| `src/ts_hsj.h` | `HierarchyBlock` struct (with `absent_state`), `hsj_score()` declaration, `partition_weights()` | Complete |
| `src/ts_hsj.cpp` | `partition_weights()`, `fitch_label_char()` (with uppass), `score_hierarchy_block()`, `hsj_score()` | Complete (full-rescore only; not wired to search pipeline) |
| `src/ts_sankoff.h` | `SankoffChar`, `SankoffData` structs, `sankoff_score()`, `sankoff_score_char()`, `sankoff_uppass()` | Complete |
| `src/ts_sankoff.cpp` | Sankoff downpass, uppass, root forcing | Complete |
| `R/recode_hierarchy.R` | `recode_hierarchy()`: x-transformation recoding (Goloboff et al. 2021) | Complete, 49 tests |
| `tests/testthat/test-recode-hierarchy.R` | Unit tests for recode_hierarchy() | Complete |
| `inst/REFERENCES.bib` | Added `Goloboff2021` entry | Complete |

---

## Modified files

| File | Change |
|------|--------|
| `DESCRIPTION` | Added `CharacterHierarchy.R` to Collate field |
| `R/MaximizeParsimony.R` | Added `hierarchy`, `inapplicable`, `hsj_alpha` params with validation |
| `src/ts_data.h` | Added `inapp_state` field to `DataSet` (for HSJ) |
| `src/ts_data.cpp` | Populate `inapp_state` in `build_dataset()` |

---

## Design decisions

- `hierarchy` is a **separate argument** to `MaximizeParsimony()` (not a phyDat attribute)
- `inapplicable` and `hsj_alpha` are **top-level args** alongside `concavity`
- Default `hsj_alpha = 1.0`
- IW + hierarchy and Profile + hierarchy: **deferred**
- Constraint interaction: **ignored** for now
- Resampling: **hierarchical** — resample top-level chars; when a controlling primary is sampled, also resample within its block; recurse for nested hierarchies

---

## Resampling with hierarchy (T-124)

`Resample()` now accepts `hierarchy`, `inapplicable`, and `hsj_alpha`
parameters. When `inapplicable != "brazeau"`, resampling is hierarchy-aware:

- **Resampling units**: each non-hierarchy character = 1 unit; each
  top-level hierarchy block (primary + all dependents) = 1 atomic unit.
- **Jackknife**: retain `proportion` of units without replacement.
- **Bootstrap**: sample `n_units` units with replacement (blocks can be
  duplicated).
- Per replicate: `.HierarchicalResampleWeights()` computes pattern weights
  for non-hierarchy chars and per-block sample counts. `.ResampleHierarchy()`
  calls `ts_driven_search` per replicate with filtered HSJ blocks or xform
  chars.
- **No C++ changes**: reuses existing `ts_driven_search` HSJ/xform infrastructure.
- **Parallelism**: serial R loop over replicates (C++ inter-search parallelism
  via `nThreads` still available within each replicate).

---

## Key algorithm notes (HSJ)

- Paper's Algorithm 1 initializes `a(l) = p(l) = 0` for all leaves. This is
  incorrect for enforcing observed leaf states. Correct initialization:
  leaf with primary absent → `a(l) = 0, p(l) = INF`; primary present →
  `a(l) = INF, p(l) = 0`. Verified against hand-computed example.
- `score_hierarchy_block()` operates per hierarchy block. Non-hierarchy
  characters use standard Fitch. Total = Fitch(non-hierarchy) + Σ HSJ(blocks).
- Secondary character labels at internal nodes from Fitch first-pass
  (inapplicable treated as a separate state).
- HSJ is full-rescore only (no incremental variant). Performance mitigation:
  candidate screening via Fitch, full HSJ only for promising candidates.

---

## Phase 2 (step-matrix/x-transform) — Complete

Sankoff engine (`ts_sankoff.h/.cpp`) implements downpass, uppass, root forcing.
R-level `recode_hierarchy()` combines primary + secondaries into composite
step-matrix character with asymmetric costs (gain:loss = n+1:1). Multistate
secondaries supported (state count = ∏k_i + 1). Nested hierarchies deferred.
Integration complete: `ScoringMode::XFORM` in `score_tree()` dispatches
Fitch(non-hierarchy) + Sankoff(recoded). `MaximizeParsimony()` accepts
`inapplicable = "xform"`. End-to-end search verified.
