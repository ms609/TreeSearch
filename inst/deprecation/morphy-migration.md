# MorphyLib Deprecation Plan — status

Last updated: 2026-07-08

## Overview

MorphyLib was a vendored C library (~3,930 LOC across 11 files) providing
Fitch parsimony scoring with support for inapplicable characters (Brazeau
et al. 2017). The native C++ engine (`ts_fitch.cpp`/`ts_fitch_na.inc` +
friends) now handles all scoring, including a correctness fix MorphyLib
lacked (`{1-}` ambiguous-with-inapplicable tokens; see
`memory/native-morphy-bgs-discrepancy.md`).

**Status: the plan below is complete on `feature/native-morphy-helpers`
(PR #251), including the full Phase D source deletion. It has not yet been
merged into `cpp-search`** — this checkout's `src/` and `R/` still carry the
full vendored library and the old `Morphy()`/`CustomSearch` API until #251
lands. Everything else this document originally tracked (Tiers 0-4, Phases
A/C/D, the R-level rearrangement functions) is done; see divergences below
for where the actual outcome differs from what was planned.

---

## Decisions that materially diverge from the original plan

1. **`Morphy()` was deleted outright, not migrated.**
   This doc's Tier 2 (and the follow-up `dev/plans/2026-06-15-morphylib-removal-stages.md`)
   treated `Morphy()` as the one function with genuinely unique functionality
   (R-loop search with per-iteration callbacks) and posed two migration
   options: (a) swap its internal scoring closures to native, or (b) redirect
   it to `MaximizeParsimony()`. Neither was taken. `Morphy()` was deleted
   wholesale (commit `d324ae149`): it was never released on CRAN, and its
   EW/IW/profile/ratchet/constraint feature set is fully covered by
   `MaximizeParsimony()` (standard criteria) or `TreeSearch()`/`Ratchet()`
   (custom criteria). `MaximizeParsimony()`'s legacy-parameter delegation now
   raises an informative error instead of calling it.

2. **Rare gap treatments (`ambiguous`/`missing`/`extra`) were dropped, not
   reproduced.** The staged plan proposed reproducing these natively by
   recoding `-` before scoring. That recode path was never built — only the
   default `inapplicable` gap treatment is supported now. The deprecated
   shims (`PhyDat2Morphy()` etc., in `R/morphy-deprecated.R`) explicitly
   `stop()` if any other `gap=` is requested.

3. **`InitializeData`/`CleanUpData` are deprecated, not just re-pointed at
   native defaults.** Phase C's implementation sketch (`CppScorerData()`/
   `CppTreeLength()`/`CppCleanUp()`) kept the injection-point parameters and
   only changed what the defaults did. The actual implementation goes
   further: `concavity` moved onto `PrepareData()` itself (a prepared
   dataset now carries its own scoring mode), and `InitializeData`/
   `CleanUpData` on `TreeSearch()`/`Ratchet()`/`Jackknife()` are deprecated
   (default `NULL`; `.Deprecated()` warning if supplied) rather than
   silently swapped. Maintainer rationale (this session): guaranteed
   on-exit cleanup is "vanishingly rare and replaceable" — a user who needs
   it can use their own `on.exit()`. The old hooks still work when supplied,
   with a warning, so this is backward compatible for one release cycle.

4. **DESCRIPTION attribution resolved alongside this work**, not as a
   separate follow-up: MorphyLib copyright dropped, Brazeau recast as `ctb`
   for the Brazeau-Guillerme-Smith inapplicable-data *algorithm* (not for
   MorphyLib code, which no longer exists in the tree).

---

## What's actually left

- **Merge PR #251 into `cpp-search`.** This is the only remaining action —
  every migration tier and phase below is complete on the feature branch.
- Nothing else identified. Tiers 0-4 and Phases A/C/D are done; Tier 3
  (`RootedSPR`/`RootedTBR`/`RootedNNI`) turned out to have no live MorphyLib
  dependency at all (they only ever manipulated edge vectors); `Sectorial*.R`
  is fully deleted rather than merely dead code; the MorphyLib-adjacent test
  files (`test-tree_length.R`, `test-CustomSearch.R`,
  `test-pp-info_extra_step.R`, `test-ts-iw.R`, `test-ts-ratchet-search.R`)
  were updated and verified free of live MorphyLib calls, while
  `test-RMorphy.R`, `test-Morphy.R` and `test-mpl_morphy_objects.R` were
  deleted along with the code they tested.

## Reference

`dev/plans/2026-06-15-morphylib-removal-stages.md` has the more granular
Stage 1/2/3 breakdown this final push followed (superseding the Tier/Phase
structure above); it is retained for history but its "NEEDS REVIEW" items
are resolved per the divergences above.
