# MorphyLib Deprecation Plan

Last updated: 2026-03-17

## Overview

MorphyLib is a vendored C library (~3,700 LOC across 10 files) that provides
Fitch parsimony scoring with support for inapplicable characters (Brazeau et al.
2017). The C++ engine (`ts_fitch.cpp` + `ts_fitch_na.inc`) now provides
equivalent or superior scoring, making MorphyLib redundant for most use cases.

This document catalogs all remaining MorphyLib dependencies and proposes a
tiered migration plan.

---

## MorphyLib C source files (candidates for removal)

| File | Lines | Contents |
|------|-------|----------|
| `src/morphy.c` | 1,126 | Core Morphy implementation |
| `src/morphy.h` | 63 | Core header |
| `src/mpl.c` | 1,038 | MorphyLib public API |
| `src/mpl.h` | 827 | MorphyLib public API header |
| `src/morphydefs.h` | 213 | Type definitions |
| `src/mplerror.h` | 54 | Error codes |
| `src/RMorphy.c` | 338 | R-MorphyLib bridge (`.Call` wrappers) |
| `src/RMorphy.h` | 35 | Bridge header |
| `src/RMorphyUtils.c` | 40 | Utilities (`preorder_morphy`, `morphy_iw`, etc.) |
| `src/RMorphyUtils.h` | 9 | Utilities header |
| `src/build_postorder.h` | 184 | Contains `RANDOM_TREE_SCORE` + `MORPHYLENGTH` |
| **Total** | **~3,930** | |

---

## R functions that call MorphyLib

### Tier 0: Already migrated to C++ (no MorphyLib dependency)

| Function | File | Notes |
|----------|------|-------|
| `MaximizeParsimony()` | `R/MaximizeParsimony.R` | C++ driven search (EW/IW/profile/constraint) |
| `AdditionTree()` | `R/AdditionTree.R` | C++ `ts_wagner_tree` |
| `Resample()` | `R/Morphy.R` | C++ `ts_resample_search` |
| `SuccessiveApproximations()` | `R/SuccessiveApproximations.R` | C++ `ts_successive_approx` |

### Tier 1: Easy to replace (C++ equivalent exists, just needs R glue)

| Function | File | MorphyLib calls | C++ replacement |
|----------|------|-----------------|-----------------|
| `TreeLength.phylo()` (EW branch) | `R/tree_length.R` | `PhyDat2Morphy`, `MorphyTreeLength` | `ts_fitch_score()` handles EW/IW/profile |
| `CharacterLength()` / `FastCharacterLength()` | `R/tree_length.R` | Per-character morphy scoring loop | `ts_na_char_steps()` returns per-pattern step counts |
| `Consistency()` / `IConsistency()` | `R/Consistency.R` | Via `CharacterLength()` / `FastCharacterLength()` | Migrates automatically when CharacterLength migrates |
| `RandomTreeScore()` | `R/RandomTreeScore.R` | `RANDOM_TREE_SCORE` C function | `ts_fitch_score()` on a `RandomTree()` |

**Estimated effort**: ~1 day. These are straightforward rewrites where the C++
function already exists and just needs to be called with the right arguments.

### Tier 2: Moderate effort (R-loop search functions)

| Function | File | MorphyLib calls | Notes |
|----------|------|-----------------|-------|
| `Morphy()` | `R/Morphy.R` | `PhyDat2Morphy`, `preorder_morphy`, `morphy_iw`, `morphy_profile`, `SingleCharMorphy`, `UnloadMorphy` | Core R-loop search. Users who need per-iteration callbacks or custom stopping criteria must use this. Could be reimplemented using `ts_fitch_score()` for scoring + R-level tree rearrangement, but the R-level rearrangement functions (SPR, TBR, NNI in R) are themselves MorphyLib-dependent (Tier 3). |
| `MorphyBootstrap()` | `R/Bootstrap.R` | `mpl_set_charac_weight`, `mpl_apply_tipdata` | Weight perturbation via MorphyLib API. Could be replaced by `ts_resample_search()`. |
| `Jackknife()` | `R/Jackknife.R` | Same weight manipulation | Same as MorphyBootstrap. `Resample()` already provides C++ jackknife; `Jackknife()` is the legacy version. |
| `Ratchet()` | `R/Ratchet.R` | `PhyDat2Morphy`, `UnloadMorphy` | Legacy ratchet. C++ `ts_ratchet_search` exists. |
| `CustomSearch()` | `R/CustomSearch.R` | `PhyDat2Morphy`, `UnloadMorphy` | Generic search framework. Powers `Morphy()`. |

**Estimated effort**: ~2-3 days. Most can be deprecated in favor of existing
C++ equivalents. `Morphy()` is the only one that provides unique functionality
(R-loop with callbacks).

### Tier 3: Hard / low priority (R-level tree rearrangement)

| Function | File | MorphyLib calls | Notes |
|----------|------|-----------------|-------|
| `RootedSPR*()` functions | `R/SPR.R` | Via scoring callbacks | R-level SPR implementation. C++ `ts_spr_search` replaces this. |
| `RootedTBR*()` functions | `R/TBR.R` | Via scoring callbacks | R-level TBR. C++ `ts_tbr_search` replaces this. |
| `RootedNNI*()` functions | `R/NNI.R` | Via scoring callbacks | R-level NNI. C++ `ts_nni_search` replaces this. |
| `Sectorial*()` functions | `R/Sectorial.R` | Commented out (dead code) | Already dead. |

**Estimated effort**: Low priority. These R-level rearrangement functions are
legacy code. The C++ search engine replaces their functionality entirely.
They could be deprecated without replacement (users should call
`MaximizeParsimony()` instead).

### Tier 4: MorphyLib API wrappers (remove last)

| Function | File | Notes |
|----------|------|-------|
| `PhyDat2Morphy()` | `R/mpl_morphy_objects.R` | Creates morphy objects from phyDat |
| `SingleCharMorphy()` | `R/mpl_morphy_objects.R` | Creates single-character morphy objects |
| `UnloadMorphy()` | `R/mpl_morphy_objects.R` | Frees morphy objects |
| `MorphyWeights()` | `R/mpl_morphy_objects.R` | Get/set character weights |
| `mpl_new_Morphy()`, etc. | `R/mpl_morphyex.R` | Thin wrappers around 20+ MorphyLib C functions |
| `MorphyTreeLength()` | `R/tree_length.R` | Direct morphy scoring |
| `MorphyLength()` | `R/tree_length.R` | Low-level morphy scoring |

These are all internal or semi-internal functions. Remove them after all
Tier 1-3 callers are migrated.

---

## Recommended migration sequence

### Phase A: Tier 1 migration (immediate, no API changes)

1. **`TreeLength.phylo()` EW branch**: Use `ts_fitch_score()` instead of
   `MorphyTreeLength()`. The IW and profile branches already avoid MorphyLib.

2. **`CharacterLength()` / `FastCharacterLength()`**: Use `ts_na_char_steps()`
   for per-pattern step counts. Expand by `attr(dataset, "index")` for
   per-character output.

3. **`RandomTreeScore()`**: Replace with `ts_fitch_score(RandomTree(nTip), ...)`.
   Function is rarely used; consider deprecating entirely.

4. **`Consistency()` / `IConsistency()`**: Migrate automatically once
   `CharacterLength()` is migrated.

### Phase B: ~~Deprecate legacy search functions~~ — RETAINED

**Decision (2026-03-19):** `Ratchet()`, `Jackknife()`, `TreeSearch()`,
`EdgeListSearch()`, `MultiRatchet()`, and `MorphyBootstrap()` are **not
deprecated**. They provide a custom search framework with pluggable
`TreeScorer` and `EdgeSwapper` functions, used by Hopkins & St John (2021).
`MaximizeParsimony()` and `Resample()` are faster for standard parsimony
but cannot accommodate arbitrary scoring functions.

Docs updated to direct standard-parsimony users to `MaximizeParsimony()` /
`Resample()` while keeping these functions available for custom criteria.

### Phase C: Decouple custom search from MorphyLib

**Assessment (2026-03-19, T-095):** The default MorphyLib scorers
(`PhyDat2Morphy`, `MorphyLength`, `UnloadMorphy`) **can** be replaced with C++
equivalents while preserving the custom search framework's pluggable API.

#### Design

The custom search framework uses dependency injection: `TreeSearch()`,
`Ratchet()`, `Jackknife()`, etc. accept `InitializeData`, `TreeScorer`, and
`CleanUpData` as function parameters. The default chain is:

```
PhyDat2Morphy(dataset) → morphyPtr
MorphyLength(parent, child, morphyPtr) → score
UnloadMorphy(morphyPtr) → NULL
```

Replace with:

```
CppScorerData(dataset) → list(contrast, tip_data, weight, levels, ...)
CppTreeLength(parent, child, data) → score  [wraps ts_fitch_score]
identity → no cleanup needed (no external pointer)
```

#### Implementation sketch

```r
CppScorerData <- function(dataset, concavity = Inf) {
  at <- attributes(dataset)
  tip_data <- matrix(unlist(dataset, use.names = FALSE),
                     nrow = length(dataset), byrow = TRUE)
  list(
    contrast = at$contrast,
    tip_data = tip_data,
    weight = at$weight,
    levels = at$levels,
    min_steps = if (is.finite(concavity)) at$min.length else integer(),
    concavity = if (is.finite(concavity)) concavity else -1.0,
    infoAmounts = if (.UseProfile(concavity)) at$info.amounts else NULL,
    original_weight = at$weight  # for bootstrap/jackknife restore
  )
}

CppTreeLength <- function(parent, child, dataset, ...) {
  edge <- cbind(parent, child)
  ts_fitch_score(edge, dataset$contrast, dataset$tip_data,
                 dataset$weight, dataset$levels,
                 dataset$min_steps, dataset$concavity,
                 dataset$infoAmounts)
}

CppCleanUp <- function(dataset) invisible(NULL)
```

#### Weight manipulation for bootstrap/jackknife

`MorphyBootstrap()` and `Jackknife()` perturb character weights via
`mpl_set_charac_weight()`. With C++ scoring, weight perturbation is simpler:
modify `dataset$weight` directly (it's an R vector, no FFI calls needed).
A replacement `CppBootstrap()` would:

1. Save `dataset$original_weight`
2. Resample: `dataset$weight <- tabulate(sample.int(nChar, replace = TRUE), nChar)`
3. Run `EdgeListSearch()`
4. Restore: `dataset$weight <- dataset$original_weight`

#### Edge ordering

`MorphyLength()` converts to postorder internally. `ts_fitch_score()` also
handles edge ordering in C++. No ordering concerns with R-level EdgeSwappers
(they output renumbered edges that both scorers accept).

#### What stays the same

- The `TreeScorer(parent, child, dataset, ...)` contract is unchanged.
- All `EdgeSwapper` functions (R-level TBR/SPR/NNI) are unaffected — they
  only manipulate edge vectors, never call MorphyLib.
- Custom user-supplied `TreeScorer` functions continue to work.
- Per-iteration R callbacks, custom stopping criteria all preserved.

#### Scoring equivalence

`ts_fitch_score()` handles all three modes (EW, IW, profile) and inapplicable
characters. It should produce identical scores to `MORPHYLENGTH` for standard
datasets. A validation step should compare scores on a test suite before
switching the default.

#### Effort and risk

- **Effort**: ~0.5 day. Three small functions + update default arguments in
  `TreeSearch()`, `Ratchet()`, `Jackknife()`, `MorphyBootstrap()`.
- **Risk**: Low. The custom search framework's API is preserved unchanged.
  Scoring delegation is the only change. A compatibility test comparing
  MorphyLib vs C++ scores on the existing test suite provides full coverage.
- **Backward compatibility**: Users who explicitly pass `TreeScorer = MorphyLength`
  continue to work (MorphyLib is still present). Only the *defaults* change.

#### Recommendation

**Option 2 (rewrite defaults)** is straightforward and low-risk. Implement it
as the next step. This eliminates MorphyLib as a *required* dependency for the
custom search framework while preserving it as an option.

After one release cycle with C++ defaults, MorphyLib can be moved to Suggests
or removed entirely (Phase D).

### Phase D: Remove MorphyLib source

Once all R functions are migrated:
1. Remove `src/morphy.c`, `src/morphy.h`, `src/mpl.c`, `src/mpl.h`,
   `src/morphydefs.h`, `src/mplerror.h`
2. Remove `src/RMorphy.c`, `src/RMorphy.h`, `src/RMorphyUtils.c`,
   `src/RMorphyUtils.h`
3. Remove `src/build_postorder.h` (contains `RANDOM_TREE_SCORE` +
   `MORPHYLENGTH`)
4. Remove all `_R_wrap_mpl_*` entries from `src/TreeSearch-init.c`
5. Remove `R/mpl_morphyex.R`, `R/mpl_morphy_objects.R`
6. Remove MorphyLib-specific tests

**Impact**: ~3,930 lines of C code removed, ~500 lines of R wrapper code
removed. Reduces compile time, binary size, and maintenance burden.

---

## Test files with MorphyLib dependencies

| Test file | MorphyLib usage | Migration notes |
|-----------|-----------------|-----------------|
| `test-RMorphy.R` | `preorder_morphy()` | Remove after Phase A |
| `test-tree_length.R` | `morphy_profile()`, `MorphyTreeLength()` | Rewrite scoring calls |
| `test-Morphy.R` | `Morphy()` search tests | Keep if Morphy() kept |
| `test-CustomSearch.R` | `Morphy()` via CustomSearch | Deprecate with CustomSearch |
| `test-pp-info_extra_step.R` | Indirect via `CharacterLength()` | Migrates automatically |
| `test-ts-iw.R` | `morphy_iw()` as reference scorer | Replace with hard-coded scores |
| `test-ts-ratchet-search.R` | `SingleCharMorphy()`, `UnloadMorphy()` | Clean up references |

---

## Risk assessment

- **Low risk**: Tier 1 migration. These are straightforward function-for-function
  replacements with existing C++ equivalents.
- **Medium risk**: Tier 2 deprecation. Legacy users of `Jackknife()`, `Ratchet()`,
  etc. need migration guidance. Deprecation warnings + one release cycle.
- **Higher risk**: `Morphy()` deprecation. Some users may depend on the R-loop
  search architecture for custom workflows. Needs a clear migration path.

## Timeline estimate

| Phase | Effort | Prerequisite |
|-------|--------|-------------|
| Phase A (Tier 1) | 1 day | None |
| Phase B (Tier 2 deprecation) | 1 day | None |
| Phase C (`Morphy()` decision) | 2-3 days | Phase A |
| Phase D (source removal) | 1 day | Phases A-C + 1 release cycle |
