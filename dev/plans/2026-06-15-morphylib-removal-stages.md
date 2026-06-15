# Dropping MorphyLib entirely from TreeSearch 2.0

Goal (maintainer decision, 2026-06-15): route ALL scoring through the native
Fitch kernel and **delete vendored MorphyLib**. Native scores are the correct
v2.0 behaviour — we are NOT preserving MorphyLib equivalence (MorphyLib mis-
scores `{1-}` ambiguous-with-inapplicable tokens; see
`memory/native-morphy-bgs-discrepancy.md`).

## Stage 1 — native `MorphyLength` for the inapplicable mode  ✅ (PR #251)
`MorphyLength()`/`MorphyTreeLength()` score the default (inapplicable) gap mode
via the native kernel (`FastCharacterLength`). Pointer retained so `Morphy()`'s
C-loop + introspection are untouched. Fixes `{1-}` for `TreeSearch`/`Ratchet`/
`Jackknife`/custom.Rmd. Verified `== TreeLength` across Lobo + CL datasets.

## Stage 2 — axe `Morphy()` entirely  ✅ (done on this branch)

Originally scoped as "migrate `Morphy()`'s C-loop", but `Morphy()` turned out to
be redundant and unreleased, so it was deleted instead. For reference, it scored
through three closures in its hot loop (former R/Morphy.R:317-345):
`Morphy()` (R/Morphy.R) scores through three closures called in the hot search
loop (R/Morphy.R:317-345):
- `.EWScore`   -> `preorder_morphy(edge, morphyObj)`
- `.IWScore`   -> `morphy_iw(edge, morphyObjs, weight, minLength, charSeq, concavity, target)`
- `.ProfileScore` -> `morphy_profile(edge, morphyObjs, weight, charSeq, profiles, target)`

**DECISION (MRS, 2026-06-15): axe `Morphy()` entirely** rather than migrate it.
It was never on CRAN (not in 1.8.0) and is the redundant middle between
`MaximizeParsimony()` (native EW/IW/profile + constraints + ratchet) and
`TreeSearch()` (custom-criteria R-loop). Functionally nothing is lost.

Done on this branch: deleted `Morphy()` + its private helpers
(`.EWScore`/`.IWScore`/`.ProfileScore`/`.TBRSearch`/`.CombineResults`/
`.ReplaceResults`/`.Time`/`.DateTime`) from `R/Morphy.R`; dropped
`export(Morphy)`; replaced `MaximizeParsimony()`'s legacy-param delegation
(`do.call(Morphy, …)`) with an informative error; deleted `test-Morphy.R`
(+ snapshots) and the `Morphy()` calls in `test-CustomSearch.R`; fixed the
dangling `[Morphy()]` doc links; regenerated NAMESPACE/man.

Consequence: `morphy_iw`/`morphy_profile` are now dead (Morphy() was their only
caller), and the rare `ambiguous`/`extra` gap modes no longer need a native
reimplementation (nothing in production uses them; `.NativeMorphyData` still
falls back to MorphyLib for them via `PhyDat2Morphy`, removed in Stage 3).

## Stage 3 — delete vendored MorphyLib  ⚠ LARGE, all-or-nothing
Once nothing calls MorphyLib:
- Make `PhyDat2Morphy()`/`SingleCharMorphy()` native-only (return the native
  handle, no pointer); `UnloadMorphy()` a no-op; `is.morphyPtr()` accepts it.
- Migrate `RandomTreeScore()` (uses `preorder_morphy`/`RandomMorphyTree`) to native.
- Remove: `preorder_morphy`/`morphy_iw`/`morphy_profile`, the `MORPHYLENGTH`
  C routine + `GetMorphyLength`/`C_MorphyLength`/`MorphyTreeLength` raw variants,
  all ~28 `mpl_*` bindings (R/mpl_morphyex.R), the introspection helpers
  (`MorphyWeights`/`SetMorphyWeights`/`GapHandler`/`summary.morphyPtr`/
  `MorphyErrorCheck`), and `RandomMorphyTree`.
- Delete the vendored MorphyLib C/C++ sources under `src/`; update Makevars /
  `R_init` registration; drop the MorphyLib copyright/attribution from
  DESCRIPTION (keep historical credit per GPL as appropriate).
- Update/replace tests: `test-RMorphy.R`, `test-mpl_morphy_objects.R`,
  `test-pp-random-tree.R`, `test-RandomTreeScore.R`, and any `Morphy()` test
  asserting exact trees/scores on inapplicable data.

### Stage 3 findings (resolved)
- **Wagner is NOT a MorphyLib dependency.** The native engine has its own
  self-contained addition-tree builder (`src/ts_wagner.cpp`:
  `wagner_tree`/`random_wagner_tree`/`biased_wagner_tree`, used by
  `ts_driven`). MorphyLib's `src/wagner.c` (`mpl_wagner_*`, `MPLstate`) is
  internal MorphyLib plumbing, reachable only through `morphy.c` → the deleted
  R bridge. `MPLstate`/`MPLndsets`/`MPLpartition` appear ONLY in MorphyLib
  files; no `ts_*` source uses them. So MorphyLib deletes wholesale — no
  `MPLstate` extraction, no Wagner rebuild.
- **Resampling weight coupling (the real gating issue).** `Jackknife()` and
  `MorphyBootstrap()` resampled by mutating character weights in place on the
  MorphyLib pointer (`mpl_set_charac_weight`) and re-scoring. Stage 1's native
  `MorphyLength` reads weights frozen into `attr(morphyObj,"native")` at
  construction, so those mutations were ignored — a latent regression.
  Fixed by **functional native resampling**: `.SetMorphyWeight()` returns a
  copy of the handle carrying the resampled weight vector; the custom-search
  framework (`EdgeListSearch` + `TreeScorer`/`EdgeSwapper`) is unchanged.
  Regression covered by a new test (doubling weights doubles the score).
  Custom-search resampling is preserved; `Resample()` (no custom TreeScorer)
  is a separate path and was not a substitute.

## Status
Stage 3 IMPLEMENTED on this branch: `PhyDat2Morphy`/`SingleCharMorphy`
native-only; `UnloadMorphy` a no-op; `RandomTreeScore` native (kept
`RandomMorphyTree` + the pure C `RANDOM_TREE`); removed `mpl_*` bindings,
introspection helpers, `MORPHYLENGTH`/`RANDOM_TREE_SCORE`/`GetMorphyLength`/
`C_MorphyLength`, `morphy_iw`/`morphy_profile`/`preorder_morphy`; deleted all
17 vendored MorphyLib C/C++ sources; trimmed `TreeSearch-init.c` +
`build_postorder.h`; regenerated RcppExports/NAMESPACE/man. Package compiles,
installs, and the test suite is green. Outstanding: DESCRIPTION MorphyLib
copyright/attribution (sensitive — Brazeau/BGS credit) left for maintainer
decision.
