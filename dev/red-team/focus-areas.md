# Red-team focus areas — TreeSearch

Rotation table for the `/red-team` skill. Built once, edited rarely. Each `/red-team`
invocation reviews **one** area (the next in rotation, see `last_focus:` at the bottom of
`log.md`) at its earned tier, then records the round in `log.md`. Verified non-trivial
findings are filed in `findings.md`. Durable lessons (bug patterns, fragile areas) live in
`../expertise/red-team.md`.

## start_tier

`start_tier` is the tier a **never-visited or freshly-rotated** area starts at. Unlike the
skill's default (everything `sonnet`), these tiers **encode measured maturity** from the
ported round history (`log.md`): areas whose seams have only ever yielded subtle,
opus-class bugs start higher; immature seams that still bleed cheap bugs start at `sonnet`.
The rotation still adjusts per recorded yield — a dry round escalates one tier, a yielding
round re-visits at the same tier with a fresh agent, a high-severity signal escalates
immediately. Treat these as the starting point, not a ceiling.

| # | Area | Files | start_tier | Key questions |
|---|------|-------|-----------|---------------|
| 1 | **Fitch scoring correctness** | `src/ts_fitch.h/.cpp`, `src/ts_fitch_na.h`, `src/ts_fitch_na_incr.h`, `src/ts_fitch_na_dirty.h` | **opus** | Does incremental / dirty-set scoring match full `score_tree()`? Bounded variants bail correctly? NA three-pass edge cases? Write a targeted test if you find a gap. |
| 2 | **Search topology invariants** | `src/ts_tbr.cpp`, `src/ts_drift.cpp`, `src/ts_search.cpp` | **opus** | After every rejected move, is topology fully restored? Undo stack correct? No stale `postorder`? Constraint metadata re-synced on *all* reject paths (incl. tabu)? Symmetry-breaking hash collisions? |
| 3 | **Ratchet & perturbation** | `src/ts_ratchet.cpp`, `src/ts_sector.cpp`, `src/ts_fuse.cpp`, `src/ts_prune_reinsert.cpp` | **opus** | `active_mask`/`upweight_mask`/`flat_blocks` fully restored after perturbation? Sectorial reinsertion reverts on worse score? `build_reduced_dataset` copies all needed fields? Fuse handles tied scores? |
| 4 | **Parallelism & RNG** | `src/ts_parallel.cpp`, `src/ts_rng.h/.cpp`, `src/ts_driven.cpp`, `src/ts_resample.cpp` | **opus** | Thread-local RNG set before any search call? **No R API (incl. `unif_rand`/`Get/PutRNGstate`) from worker threads** — note the resample path. Pool mutex correct? Atomic stop-flag races? Seeds drawn from R RNG before spawn? |
| 5 | **Data pipeline & simplification** | `src/ts_data.h/.cpp`, `src/ts_simplify.h/.cpp` | **opus** | `build_dataset` handles edge cases (all-ambiguous, single-state, zero-weight, `n_states==32` UBSAN)? `build_reduced_dataset` copies all fields? XPIWE `obs==0` division? |
| 6 | **R ↔ C++ interface** | `src/ts_rcpp.cpp`, `src/TreeSearch-init.c`, `R/RcppExports.R`, `R/MaximizeParsimony.R`, `R/SearchControl.R` | **sonnet** | Arg counts match? Concavity sentinel translated? Edge-matrix conventions? Return value attributes/types set (frozen-API `logical` vs `integer`)? Parameter validation in R layer? |
| 7 | **Shiny module wiring** | `inst/Parsimony/server.R`, `inst/Parsimony/server/mod_*.R`, `inst/Parsimony/server/events.R` | **sonnet** | Forward-ref callbacks resolve? Cross-module `updateXxxInput` namespaces correct? Re-entrancy / double-launch guards? Stale dataset-hash on async tasks? `onStop` cleanup (cancel signal + temp files)? Orphaned observers? |
| 8 | **Test suite health** | `tests/testthat/test-ts-*.R`, `tests/testthat/helper-ts.R` | **sonnet** | Tier guards correct? Vacuous (always-pass) assertions? Missing `TreeSearch:::` prefixes? `set.seed()` before `sample()`? Edge-case coverage gaps (3-tip, single-char, all-NA)? Enduring regression for incremental-rescore? |
| 9 | **Wagner & addition trees** | `src/ts_wagner.h/.cpp`, `R/AdditionTree.R`, `R/PolEscapa.R` | **opus** | NA-incremental scoring staleness acceptable? Constraint mapping (LCA-based) correct? Retry loop fires? 3-taxon base case handles all orderings? R-layer index/`sequence` validation (OOB-write guard)? |
| 10 | **Profile & IW scoring** | `src/ts_fitch.cpp` (IW/profile paths), `src/ts_data.cpp` (precompute) | **opus** | `e/(k+e)` delta correct? Profile `info_amounts` lookup + capping matches? `concavity = 1.0` sentinel activates weighted path? `precompute_profile_delta` includes `precomputed_steps` offset? Clipped-subtree homoplasy in screening? |
| 11 | **Zero-length-branch collapse (MPT set)** | `src/ts_collapsed.cpp/.h`, `src/ts_splits.cpp` (`compute_collapsed_splits`), `src/ts_rcpp.cpp` (`ts_collapse_flags_batch`), `src/ts_tbr.cpp` (enum `add_collapsed` sites), `R/MaximizeParsimony.R` (collapse block) | **opus** | DEFAULT-ON since 2026-06-24, so every `MaximizeParsimony` call exercises it. Does `compute_collapsed_flags_aggressive` flag the *correct* min-length-0 branches under **IW / profile / NA**, not just EW (verified)? Is it really rooting-invariant, or does tip-rooting+`RenumberTips(labs)` alignment break on constraint trees / user start trees / `RenumberTips` permutations (cf. [[na-validation-alignment-gotcha]])? Can the dedup key `write.tree(SortTree(unroot(t)))` over-merge (two distinct collapsed topologies → same key) or under-merge across rootings? `result$scores == best_score` float-equality safe under IW/profile? Degenerate inputs: star tree, single MPT, 3–4 tips, all-resolved (must be exact no-op), fully-unresolved? Does collapse ever produce a tree that violates an active `constraint`? |
| 12 | **Red-team process meta-review** | `dev/red-team/focus-areas.md`, `dev/red-team/log.md`, `dev/red-team/findings.md` | **sonnet** | Are any areas too broad — spanning multiple distinct seams such that a finder concentrating on one file family misses another? Are any too narrow — a single-feature scope that would be better merged into a neighbour? Do any areas overlap (same source files audited under two different area headings)? Has any area gone persistently dry (≥ 3 consecutive rounds with zero confirmed findings) — should it be retired, merged, or downtiered? Are there new code seams (recently merged features, new source files) not covered by any existing area? Are tier assignments calibrated to actual yield recorded in `log.md` — any area that keeps surprising at its current tier and should escalate, or one that has been consistently empty and should drop? Propose concrete restructuring actions (split, merge, retire, add, re-tier) with rationale tied to `log.md` yield history. |
| 13 | **Constrained search correctness** | `src/ts_constraint.h/.cpp`, `src/ts_nni_perturb.cpp`, constraint integration points in `src/ts_driven.cpp` (fuse), `src/ts_parallel.cpp` (parallel-fuse), `src/ts_wagner.cpp`/`src/ts_sector.cpp` (posthoc retry), `src/ts_tbr.cpp` (`regraft_violates_constraint`) | **opus** | Does every `impose_constraint()` caller verify-before-capture, not just trust an improved score (T-213 gap, fixed d9a4f827: `nni_perturb_search` was the one caller that didn't re-check `constraint_node[]` after repair — fuse/parallel-fuse already did)? Any other heuristic-repair or posthoc-retry caller (Wagner build retry, sector) that skips discard-on-failure? Is `impose_one_pass`'s `best_node` reference stale after its own move-out loop's `topology_spr()` calls relocate a node — traced mechanism, produced one `std::bad_alloc` crash under experimental code, did NOT reproduce in 600 stress-test seeds against shipped code; needs a targeted adversarial tree construction, not more random seeds, to confirm either way. Is `map_constraint_nodes`/DFS-timestamp resync correct on every topology-mutation path, including reject paths (cross-check vs area 2's tabu-reject question)? Are nested/overlapping constraint splits handled consistently across TBR clip-gating, Wagner retry, and sector/fuse posthoc paths? |

### Maturity / tier rationale (one line each)

- **1 Fitch correctness — opus.** Crown jewel; T-300 (systematic delta=−3) and T-306 were
  opus-class subtle bugs. Prime **fable**-escalation target the moment opus runs dry.
- **2 Topology invariants — opus.** Deep state-restore subtleties; T-235 (SPR stale state),
  T-316 (P1 stale constraint metadata after tabu rejection).
- **3 Ratchet & perturbation — opus.** Mature, but the `build_reduced_dataset` /
  mask-sync class (T-275, T-303) keeps recurring.
- **4 Parallelism & RNG — opus.** Looked mature until T-309 (P1: R RNG API on worker
  threads via the parallel `Resample()` path, 2026-06-15). Concurrency = P1-capable.
- **5 Data pipeline — opus.** DAT-001 (`1u<<32` UBSAN at `n_states==32`), DAT-002 (XPIWE
  `obs==0` division). Edge inputs still bite.
- **6 R ↔ C++ interface — sonnet.** Mostly mechanical arg-count / sentinel audits; cheap to
  run. T-310 (frozen-API `pruneReinsertNni` type) shows it still occasionally yields —
  escalate if a sonnet pass goes dry.
- **7 Shiny wiring — sonnet.** **Immature seam:** a Sonnet pass found 5 bugs on 2026-06-16
  (T-309…T-313). Keep mining cheap until it runs dry.
- **8 Test suite health — sonnet.** Reliably yields inline fixes (`set.seed`, vacuous
  asserts) and test-gap notes (T-304).
- **9 Wagner & addition — opus.** Kernel code; WGN-01 (P1 OOB write via `AdditionTree(sequence=)`).
- **10 Profile & IW — opus.** Numerical delta algebra; subtle conservative bugs (profile
  delta capping). Secondary **fable**-candidate alongside area 1.
- **11 Zero-length collapse — opus, NEVER REVIEWED.** New + default-on; implementation hit
  three subtle traps in one sitting (conservative flags rooting-sensitive; aggressive flags
  need tip-rooting; tip-data alignment via `RenumberTips`). Cross-mode correctness
  (IW/profile/NA) and dedup canonicalization are the under-verified seams. Start opus.
- **12 Red-team meta-review — sonnet, NEVER REVIEWED.** Pure document review; no code to
  trace. Cheap per round; findings are restructuring proposals (split/merge/retire/add/re-tier)
  that improve every future round. Escalate to opus only if evaluating a proposed split
  requires reading source files to assess scope boundaries.
- **13 Constrained search — opus, NEVER REVIEWED.** Split out 2026-07-02 after fixing T-213
  (`nni_perturb_search` missing verify-before-capture on `impose_constraint()`'s heuristic
  repair — cf. [[impose-constraint-verify-gap]]). Same session traced but did NOT confirm a
  second, more severe failure mode (stale `best_node` across `topology_spr()` relocation
  inside `impose_one_pass`, hand-derived from a crash under experimental code, absent after
  a 600-seed stress test against shipped code) — exactly the kind of subtle, heuristic +
  topology-mutation state-invariant bug that areas 1/2/9/10/11 keep finding at opus tier.
  Start opus; a confirming repro of the stale-`best_node` mechanism would be the highest-value
  first finding.
