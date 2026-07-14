# Lever #7 — plain-EW scorer monomorphization (template-by-criterion-flavour)

**Status:** SPEC, ready to build. Off `cpp-search`. Land #7 BEFORE #6 (this is the tractable one).
**Written 2026-07-14** from the per-move investigation (see memory node `tnt-per-move-kernel-gap`).

## 0. Mission (and the honest ceiling)
Land an EXACT (byte-identical) per-candidate speedup in the TBR scorer by lifting per-candidate
**criterion-flavour dispatch out of the hot loop** — instantiate one monomorphic candidate-scan loop
per flavour, chosen once at dispatch, so the dead branches compile away and the compiler can inline /
pick a specialised kernel.

- **Measured ceiling:** a *runtime* branch-strip alone ≈ **5% whole-search**; the **template** version
  additionally unlocks a weight-aware "flat" kernel worth **~17% of the SPR scan** on unit-weight data
  → realistic target **~5–10% whole-search, EXACT**. We take the 5%.
- **Orthogonal to 5432 reach** — this is a general wall-clock-to-optimum lever; it does NOT change which
  trees are found. Do not sell it as a 5432 fix.

## 1. Read first (durable evidence — all on disk, survive compaction)
- Memory `tnt-per-move-kernel-gap` — the whole per-move investigation + the "#7 MEASURED" block. Key
  facts: the ~24× per-move gap is NOT arithmetic (T-250 + 64-bit disasm: TS SSE2 ≥ TNT SSE2), NOT
  strippable orchestration alone (~5%); the reroot x4 path is ALREADY branch-free unit-weight (why
  77–81% of candidates "had nothing to strip").
- `dev/profiling/s7-fastpath-sizing.md` (worktree branch `worktree-agent-af64d954c25f78a8e`) — the sizing:
  dead-flavour-branch strip = ~4–6% of SPR; the unit-weight **flat** scorer (drops CharBlock deref +
  weight-multiply + active_mask) = ~17% of SPR ≈ 4.7% whole-search but **WEIGHT-BLIND (unsafe as a global
  toggle)**. A first A/B falsely read ~0% because the fast path never fired (the collapsed-vector is
  always sized, so the `!use_collapsed` gate was BASE-vs-BASE) — a fire-counter + a deliberately-wrong
  mode caught it. **Replicate that positive control.**
- Prototype to generalise from: the `TS_EW_STRIP` runtime probe (default-OFF, uncommitted) in
  `src/ts_tbr.cpp` on that worktree branch — a specialised branch-free plain-EW SPR loop, proven
  byte-identical (BASE==STRIP==307). #7 turns this from a runtime toggle into compile-time specialisation.

## 2. The lever (design — this is the user's framing, and it's the right one)
The per-candidate loop re-tests, EVERY candidate, checks that are **constant for the whole search**:
- **Criterion flavour (constant per search):** `has_na` (standard Fitch vs NA three-pass), `use_iw`
  (EW vs implied weights), `all_weight_one`/`use_flat` (unit vs weighted blocks). You never switch
  EW↔IW mid-search.
- **Search-mode (constant within a clip loop):** `use_collapsed`/`collapsed`, `sector_mask`,
  `constrained`, tabu, pool, `screen_probe`, `b2_ceiling`, `iw_timing`.

**Monomorphize:** template the candidate-scan loop on these as `bool` template params, instantiate once
at dispatch (a small dispatch switch selects the instantiation), so the compiler eliminates the dead
branches and better inlines/allocates. "A template, called once, per flavour of optimality criterion."

**The weight-aware kernel is where #7 beats the bare 5%** — and why the *template* form matters, not a
runtime toggle: the flat unit-weight kernel is unsafe as a global switch (weight-blind), but SAFE as an
instantiation gated on `all_weight_one` — the compiler picks the flat kernel only for unit-weight data,
the general (weighted) kernel otherwise. That safely captures the ~17%-of-SPR piece.

**Note TS already does a slice of this:** the reroot x4 kernel is a branch-free unit-weight path. #7 =
extend the same monomorphization to (a) the **SPR scan** and (b) the remaining flavour branches on the
reroot path. Cover **both** candidate paths.

Scope: **EW first** (the hot/common case, the win). IW/NA instantiations can follow but are not the prize.

## 3. Correctness gate (MANDATORY — before any wall claim)
Monomorphization is a pure reorganization ⇒ final `score` AND `attr(result,"candidates_evaluated")` AND
per-pass `n_candidates_evaluated` MUST be **BYTE-IDENTICAL** vs baseline on Wortley2006 / Zhu2013 /
Zanol2014 × ≥2 seeds (EW: gaps `-`→`?`), via BOTH `TreeSearch:::ts_tbr_diagnostics` (kernel-level) AND
full `MaximizeParsimony` (exercises the ratchet `use_flat=false` + `sector_mask` regimes the kernel
oracle can't reach). **Replicate the positive control:** a fire-counter proving the specialised path
actually executed, plus a deliberately-wrong instantiation that DIVERGES — otherwise a no-op passes the
gate falsely (the agent's ~0% trap). Any drift = bug: STOP, do not proceed to timing.

## 4. Measurement
Per-candidate ns via an in-process `std::chrono` around the scoring region (reuse the `TS_IW_TIMING`
hook), single-thread, **min-of-runs** (load-robust), baseline vs templated, local panel (Zhu 75t /
Zanol 74t / Dikow 88t). Report SPR-scan % and whole-search %. Honest expectation ~5–10% whole-search;
the flat-kernel-on-unit-weight is the bulk. 5432/482t confirmation is Hamilton (optional; the sign is
established on the panel).

## 5. Build & standing constraints
- **Own git worktree** off cpp-search; commit NAMED files only; **NEVER push**; **NEVER** touch/merge/
  checkout `cpp-search` or `main`; no `git add -A` (cpp-search is shared by parallel sessions).
- Per-agent `R CMD INSTALL --library=.agent-X` or `dev/build-fast.R` (ccache); **never** `load_all` for
  anything measured.
- Templates touch headers ⇒ ccache can serve **stale .o** ⇒ build `CCACHE_DISABLE=1 … --preclean` and
  validate (the stale-object ABI gotcha).
- This machine's `~/.R/Makevars.win` zeroes PKG_CXXFLAGS ⇒ put any `-D` in **PKG_CPPFLAGS** and grep the
  build log to confirm the flag is live.
- No `thread_local` on the hot path (MinGW emutls); file-scope static for single-thread timing only.
- R conventions: camelCase vars, BigCamelCase funcs, no `<<-`/right-assign, base-R over tidyverse,
  TreeTools over ape.
- **Anchors drift — verify against the current tip** (`git show HEAD:src/ts_tbr.cpp`): the SPR scan and
  reroot loop live in `tbr_search` (src/ts_tbr.cpp); scorers `fitch_indirect_length_cached` + the
  flat/x4 kernels (src/ts_fitch.cpp); flavour flags on `DataSet` (`all_weight_one`, `n_blocks`,
  `flat_blocks`, ts_data.h).

## 6. Deliverable
The templated EW scorer (either as the default path — it's exact — or flag-gated for a staged land),
the bit-identity gate result + positive control, the measured per-candidate + whole-search delta, and a
durable `dev/profiling/s7-land.md`. Committed to your worktree branch, NOT pushed — supervisor merges.
