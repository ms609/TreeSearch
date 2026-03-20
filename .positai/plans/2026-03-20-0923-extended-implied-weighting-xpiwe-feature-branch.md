# Plan: Extended Implied Weighting (XPIWE) — feature/xpiwe

Date: 2026-03-20
Branch: `feature/xpiwe` (worktree `../TS-Xpiwe`, branched from `cpp-search`)

---

## Background

Goloboff (2014) "Extended implied weighting" (*Cladistics* 30: 260–272,
doi:10.1111/cla.12047) proposes a modification to the IW fit function that
explicitly accounts for the minimum possible steps for each character.

**Standard IW** (Goloboff 1993) minimises:

    Σ_c  w_c · e_c / (k + e_c)

where `e_c = s_c − m_c` (extra steps = actual − minimum), `k` = concavity.

**XPIWE** replaces `k` with a per-character effective concavity `k + m_c`:

    Σ_c  w_c · e_c / (k + m_c + e_c)

Equivalently, maximise `(k + m_c) / (k + m_c + e_c)` instead of
`k / (k + e_c)`.

Rationale: a character that requires many steps on any tree (`m_c` large) is
intrinsically "harder" and should have its extra homoplasy penalised less
steeply. Standard IW ignores `m_c` in the denominator, which also means
characters with lots of missing data (which artificially compress `m_c` toward
zero) are penalised more than they should be.

> **⚠ Formula confirmation needed:** The formula above is consistent with
> secondary sources and with the TNT implementation, but should be verified
> against the paper (doi:10.1111/cla.12047) before coding begins.

---

## Performance cost

The formula change (`k` → `k + m_c` per pattern) touches three places:

| Function | Change | Overhead |
|---|---|---|
| `compute_iw()` | `k + extra` → `eff_k[p] + extra` | ~0 — one extra add per pattern |
| `precompute_iw_delta()` | `k + e` → `eff_k[p] + e` | ~0 — same |
| `indirect_iw_length_*()`  | **no change** — consumes pre-computed delta | zero |

The `eff_k[p]` vector is precomputed once when the `DataSet` is built.
Because the indirect scoring functions only consume `iw_delta[pattern]`
(which is precomputed from `eff_k[p]`), the hot TBR/drift inner loops need
no modification.

**Conclusion: XPIWE has essentially zero runtime overhead vs standard IW.**
Any benchmark difference will be in the noise.

---

## Decision 1 — Should XPIWE be the default?

**✅ DECIDED: Yes — XPIWE is the default for 2.0.0.**

2.0.0 is already a declared major version bump, so a breaking change in IW
behaviour is appropriate. Standard IW (`extended_iw = FALSE`) remains
available for users who need to reproduce pre-2.0.0 results.

Impact:
- `concavity = Inf` (EW): unaffected — `extended_iw` is irrelevant when
  there is no minimum-step normalisation.
- `concavity = "profile"`: unaffected.
- `concavity = <numeric>` (IW): **behaviour changes.** Users get XPIWE by
  default; standard IW requires explicit `extended_iw = FALSE`.
- NEWS entry must clearly flag this as a breaking change.

---

## Decision 2 — API parameter design

**✅ DECIDED: Option A — boolean flag, default `TRUE`.**

```r
MaximizeParsimony(dataset, concavity = 10, extended_iw = TRUE)
TreeLength(tree,  dataset, concavity = 10, extended_iw = TRUE)
Resample(dataset,           concavity = 10, extended_iw = TRUE)   # Phase 2
SuccessiveApproximations(dataset, concavity = 10, extended_iw = TRUE)  # Phase 2
```

Rules:
- `extended_iw = TRUE` activates XPIWE when `is.finite(concavity)`.
- When `concavity = Inf` or `concavity = "profile"`, `extended_iw` is
  ignored (no message — it is simply inapplicable, not wrong).
- `extended_iw = FALSE` gives the legacy Goloboff (1993) formula.
- C++ `make_dataset()` receives `bool xpiwe` (always explicit from R; C++
  default stays `false` to keep the C++ API conservative).

`SearchControl()` exposure: **include in Phase 1** — low cost, and expert
users calling `SearchControl()` deserve the same access.

---

## Decision 3 — Scope of this feature branch

### In scope (Phase 1 — this branch)

1. `ScoringMode::XPIWE` added to `ts_data.h` enum.
2. `DataSet::eff_k` (per-pattern double vector) added; populated in
   `build_dataset()` as `k + global_min[p]` for XPIWE or `k` for IW.
   - `global_min[p] = min_steps[p] + precomputed_steps[p]`
     (adjusted min steps plus the topology-independent offset removed during
     simplification).
3. `compute_xpiwe()` and `precompute_xpiwe_delta()` in `ts_fitch.cpp`
   (or unify with IW via `eff_k` — see C++ design below).
4. Dispatch in `compute_weighted_score()` and `precompute_weighted_delta()`.
5. `make_dataset()` in `ts_rcpp.cpp`: new `bool xpiwe = false` parameter;
   sets `ScoringMode::XPIWE` and builds `eff_k`.
6. **R functions updated**: `MaximizeParsimony()`, `TreeLength()`.
7. Tests: `test-ts-xpiwe.R` (Tier 2), formula unit tests, search end-to-end.
8. Docs: parameter documented in `@param`, vignette note.

### Deferred (Phase 2 — same or follow-on branch)

- `Resample()` + `SuccessiveApproximations()` — these call `ts_driven_search`
  and already pass `concavity`; adding `xpiwe` is straightforward but
  lower priority.
- `SearchControl()` — expose `extended_iw` for expert parameter bundles.
- Profile + XPIWE — remains incompatible (profile has its own correction).
- HSJ/xform + XPIWE — non-hierarchy characters could use XPIWE, but deferred.
- ParsSim — uses `TreeLength()` so Phase 1 already covers simulation scoring.

---

## Decision 4 — C++ implementation strategy

### Recommended: per-pattern `eff_k` vector in DataSet

Add to `DataSet`:
```cpp
std::vector<double> eff_k;  // per-pattern effective k; k for IW, k+m for XPIWE
```

Populate in `build_dataset()` after setting `scoring_mode`:
```cpp
ds.eff_k.resize(n_patterns);
for (int p = 0; p < n_patterns; ++p) {
  ds.eff_k[p] = ds.concavity;
  if (ds.scoring_mode == ScoringMode::XPIWE) {
    ds.eff_k[p] += ds.min_steps[p] + ds.precomputed_steps[p];
  }
}
```

Then unify `compute_iw` to use `eff_k[p]`:
```cpp
double compute_iw(const DataSet& ds, const std::vector<int>& char_steps) {
  double score = 0.0;
  for (int p = 0; p < ds.n_patterns; ++p) {
    int extra = char_steps[p] - ds.min_steps[p];
    if (extra > 0) {
      score += ds.pattern_freq[p] *
               (static_cast<double>(extra) / (ds.eff_k[p] + extra));
    }
  }
  return score;
}
```

Similarly for `precompute_iw_delta`. The indirect functions remain
**unchanged** — they only consume `iw_delta[p]`, which is now implicitly
XPIWE-aware via `eff_k`.

**Advantage**: no new functions needed; IW and XPIWE share a single code
path; the `ScoringMode::XPIWE` enum value signals that `eff_k` was
populated with `k + m[p]`.

**Alternative**: separate `compute_xpiwe` / `precompute_xpiwe_delta`
functions with branching dispatch in `compute_weighted_score`. Slightly
more code but avoids touching the IW functions directly. The per-pattern
`eff_k` approach is cleaner.

---

## Decision 5 — Interaction with `precomputed_steps` / offset

`ds.min_steps[p]` is the adjusted minimum (= global_min − precomputed_steps).
For XPIWE we need `global_min[p]` in the denominator.

`global_min[p] = ds.min_steps[p] + ds.precomputed_steps[p]`

Both vectors are already present in `DataSet`. No new data needs to be
passed from R; `eff_k` is computed entirely within `build_dataset()`.

---

## Open questions / decisions for the user

1. **Formula confirmation** *(still pending)*: Please verify from Goloboff
   (2014) §2 (doi:10.1111/cla.12047) that the XPIWE fit function is
   `(k + m) / (k + m + e)`, i.e. the loss being minimised is
   `e / (k + m + e)`.  Secondary sources and TNT's implementation both
   point to this, but confirming from the primary source before committing
   the formula to code removes any risk.

2. **Scope of `Resample()` / `SuccessiveApproximations()` in Phase 1**:
   These functions call `ts_driven_search` (already updated in Phase 1);
   the R-side change is straightforward — add `extended_iw = TRUE` and
   pass it through.  Worth doing now to ship a consistent API, or defer?

3. **Validation for non-IW modes**: When `concavity = Inf` and
   `extended_iw = TRUE`, silently ignore (recommended) or emit a
   `message()`? The parameter is vacuously correct (EW is unchanged
   regardless of m), so a warning would be noise for users who set
   `extended_iw = TRUE` globally and then switch to EW.

---

## Implementation checklist (Phase 1)

### C++ (`src/`)

- [ ] `ts_data.h`: Add `std::vector<double> eff_k;` to `DataSet`;
  add `XPIWE` to `ScoringMode` enum.
- [ ] `ts_data.cpp`: Populate `ds.eff_k` in `build_dataset()` per above.
- [ ] `ts_fitch.cpp`: Change `compute_iw` and `precompute_iw_delta` to use
  `ds.eff_k[p]` instead of `ds.concavity`.
- [ ] `ts_fitch.h`: Update declarations / comments.
- [ ] `ts_rcpp.cpp`: Add `bool xpiwe = false` param to `make_dataset()` and
  to `ts_fitch_score`, `ts_driven_search`, `ts_successive_approx`,
  `ts_resample_search` (at minimum the first two for Phase 1).
- [ ] `TreeSearch-init.c`: Re-run `Rscript check_init.R` after any signature
  change; update arg counts as needed.

### R (`R/`)

- [ ] `MaximizeParsimony.R`: Add `extended_iw = TRUE` parameter; ignore
  silently when `!is.finite(concavity)` or profile; pass to C++.
- [ ] `tree_length.R`: Add `extended_iw = TRUE` to `TreeLength()` and all
  S3 methods; pass through to `ts_fitch_score`.
- [ ] `SearchControl()`: add `extended_iw = TRUE` (Phase 1).
- [ ] `Bootstrap.R` / `Resample.R`: add `extended_iw = TRUE` (Phase 1 if
  scope confirmed, Phase 2 otherwise).

### Tests (`tests/testthat/`)

- [ ] `test-ts-xpiwe.R` (Tier 2, `skip_on_cran()` at file level):
  - Unit test: `compute_iw` vs `compute_xpiwe` differ when `min_steps > 0`.
  - Hand-computed example: 5-taxon dataset, known `m` per character,
    verify XPIWE score = EW score when `m = 0`, verify XPIWE ≠ IW when `m > 0`.
  - Formula test: XPIWE with `m = 0` equals IW (same formula).
  - End-to-end: `MaximizeParsimony(dataset, concavity = 10, extended_iw = TRUE)`
    runs and returns valid trees.
  - Score comparison: XPIWE score ≠ IW score on dataset with informative
    multistate characters (non-trivial min_steps).
  - `TreeLength` with `extended_iw = TRUE` gives same score as
    `MaximizeParsimony` on the best tree.
- [ ] Extend `test-iw-scoring.R` to include cross-check that IW score is
  unchanged by the `eff_k` refactor (regression).

### Docs

- [ ] `@param extended_iw` entries in `MaximizeParsimony.Rd` and
  `TreeLength.Rd`.
- [ ] Short paragraph in `vignettes/profile-scores.Rmd` explaining when to
  prefer XPIWE over standard IW.
- [ ] NEWS entry.
- [ ] Run `devtools::check_man()` and `spelling::spell_check_package()`.

---

## Worktree setup

```bash
# From TreeSearch-a (on cpp-search)
git worktree add ../TS-Xpiwe feature/xpiwe
# or if branch doesn't exist yet:
git worktree add -b feature/xpiwe ../TS-Xpiwe cpp-search
```

Build into `.agent-X/` (rename X to appropriate agent letter for this
worktree, e.g. `.agent-pman/` for the primary user):

```bash
SRC=$(pwd) && TMPBUILD=$(mktemp -d) && \
  rm -f src/*.o src/*.dll && \
  (cd "$TMPBUILD" && R CMD build --no-build-vignettes --no-manual \
     --no-resave-data "$SRC") && \
  R CMD INSTALL --library=.agent-pman "$TMPBUILD"/TreeSearch_*.tar.gz && \
  rm -rf "$TMPBUILD"
```

---

## Summary of key up-front decisions

| Decision | Choice | Status |
|---|---|---|
| Formula: `e/(k+m+e)` | Yes — verify from paper before committing | **Formula confirmation pending** |
| Make XPIWE the default? | **Yes — `extended_iw = TRUE`** | ✅ Decided |
| Parameter name | **`extended_iw`** (Option A boolean flag) | ✅ Decided |
| C++ strategy | `eff_k` vector; unify IW/XPIWE path | ✅ Decided |
| Phase 1 scope | `MaximizeParsimony` + `TreeLength` + `SearchControl` | ✅ Decided |
| `Resample` / `SA` in Phase 1? | TBD | **Pending user input** |
| Silent-ignore vs message for EW + `extended_iw` | TBD | **Pending user input** |
