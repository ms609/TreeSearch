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

## Stage 2 — migrate `Morphy()`'s C-loop + the rare gap modes  ⚠ NEEDS REVIEW
`Morphy()` (R/Morphy.R) scores through three closures called in the hot search
loop (R/Morphy.R:317-345):
- `.EWScore`   -> `preorder_morphy(edge, morphyObj)`
- `.IWScore`   -> `morphy_iw(edge, morphyObjs, weight, minLength, charSeq, concavity, target)`
- `.ProfileScore` -> `morphy_profile(edge, morphyObjs, weight, charSeq, profiles, target)`

Migration options (DECISION NEEDED — both change superseded-path behaviour, and
native!=MorphyLib on `{1-}` data means `Morphy()` may return different trees):
  (a) **Swap the closures to native** keeping the existing search loop:
      - EW: trivial — `morphyObj` already carries the native attr; call
        `.NativeMorphyLength`.
      - IW: `ts_fitch_score(edge, contrast, tip_data, weight, levels, min_steps,
        concavity = k)`. **Loses** the `target` early-abort bound (efficiency,
        not correctness) and the per-char `charSeq`/`morphyObjs` machinery.
      - Profile: `ts_fitch_score(..., infoAmounts = profiles)` similarly.
        `morphy_profile`'s `target` early-abort is lost.
  (b) **Redirect `Morphy()` to `MaximizeParsimony()`** — `Morphy()` is fully
      superseded for EW/IW/profile. Cleanest end-state, but must map args
      (`ratchIter`/`tbrIter`/`maxHits`/`concavity`/constraint) and accept that
      the search path (and trees found) change.
Recommendation: (a) if `Morphy()`'s search API must be preserved verbatim;
(b) if a behaviour change to this superseded function is acceptable.

Rare gap modes (currently fall back to MorphyLib via `.NativeMorphyData` ->
NULL): reproduce natively by recoding the dataset before scoring —
  - `ambiguous`/`missing`: matrix round-trip, `-` -> `?` (drops the `-` level).
  - `extra`: `-` -> a fresh ordinary state symbol.
Then `.NativeMorphyData`/`.SingleCharNative` return data for all modes.

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

## Why Stages 2-3 are not done unattended
They change behaviour of superseded public functions, span the C/C++ build, and
require many test rewrites + CI cycles — they warrant a maintainer decision on
the `Morphy()` approach (a vs b) and review before landing.
