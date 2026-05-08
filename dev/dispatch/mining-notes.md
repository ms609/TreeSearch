# Mining Notes: Legacy PositAI Artifacts

## Surviving facts from agent-*.md

### In-progress / Parked Tasks

- **agent-c.md: T-214 PARKED on GHA 23536512228** — Multi-split constraint enforcement bug during TBR search. Root cause identified: `classify_clip_constraints()` marks clips as UNCONSTRAINED incorrectly when constraint tips and extras straddle attachment edge. Two-part fix implemented (post-hoc `map_constraint_nodes()` + FORBIDDEN clip zone). Added test-ts-constraint-multi.R (806 assertions). Needs GHA result.

- **agent-e.md: T-289f PARKED — GHA 23690338955 (feature/tbr-batch); Hamilton down** — Prune-reinsert PR NNI polish cost reduction. Stage 5 submitted as SLURM 16622224. Root cause of Stage 4 failure: full TBR convergence after each PR cycle (~7s per 5 cycles). New SearchControl() params added: `pruneReinsertNni` (NNI vs TBR polish) and `pruneReinsertFullMoves` (limit full-tree TBR). Stage 5 results indicate pr_nni wins 7/10 conditions; benefit dataset-dependent, reverses at >=206t. Feature not enabled in large preset, available via SearchControl().

### Critical Findings

- **agent-a.md: TS-PruneRI directory orphaned** — After T-266 completion and branch deletion, local git metadata removed but directory remains (manual cleanup needed).

- **agent-a.md: T-204 fix complexity** — GHA 23641482723 failed due to T-204's `.Deprecated()` addition to `PhyDat2Morphy`/`UnloadMorphy` causing warnings in examples (PhyDat2Morphy.Rd, MorphyWeights.Rd, GapHandler.Rd, SingleCharMorphy.Rd, Morphy.R constraint example). Fixed via WORDLIST updates and `\donttest{}`/`suppressWarnings()` wrappers.

- **agent-a.md: S-RED focus 10 bug fixed** — precompute_profile_delta had old_cost=0 when s>info_max_steps. Fixed in commit 7cff7870 (15 tests pass).

- **agent-a.md: PR #213 (cid-consensus) aborted** — GHA conflict: ts_tbr.cpp between CID and T-263 snapshot. Needs E/human review.

- **agent-d.md: S-RED focus 4 — consensus stability bug in parallel path** — Idle polls incorrectly increment unchanged counter → premature termination. Identified and fixed.

- **agent-g.md: G-006 filed** — nni_search in ts_prune_reinsert.h/.cpp lacks ConstraintData* parameter (found during S-RED Focus 30-31).

### Completed & Merged

- **agent-a.md: T-266 (PR #235)** — Taxon pruning-reinsertion perturbation strategy. Commit afbf531f. Phase distribution: Ratchet 46.3%, NNI-perturb 34.3%, RSS 7.4%, CSS 4.4%, XSS 3.2%, TBR 3.2%. T-274 filed for benchmarking nniPerturbCycles=0 vs 5.

- **agent-a.md: T-270 (vignette docs)** — Completed; updated vignettes/search-algorithm.Rmd (new pipeline step 5a, post-ratchet sectorial subsection). Commit d8f3c769.

- **agent-b.md: T-277 (PR #236 open)** — ScoreSpectrum() Chao1 landscape coverage estimator. Awaiting human review/merge.

- **agent-b.md: T-275, T-230, T-235, T-226 completed** — Prune-reinsert EW guard, replicate-count warning gate, full_rescore after rejected SPR regraft, remove "Trees in sequence" option.

- **agent-f.md: F-030 (PR #239, merged)** — TBR clip-ordering Phase 2. Feature/weighted-clip-order deleted; worktree TS-WeightClip pending manual deletion.

- **agent-f.md: T-245 (PR #238, merged)** — TBR 4-wide candidate batching.

- **agent-g.md: T-289f Stage 5 complete** — Prune-Reinsert NNI vs TBR Polish benchmark (SLURM 16622421, 7h). Five large-tree datasets (131-206t), 20 seeds, EW scoring. pr_nni wins 7/10 conditions. Not enabled in large preset (benefit dataset-dependent). Strategies.md updated.

- **agent-g.md: T-290c** — wagnerStarts=1 vs 3 under Brazeau scoring (2 datasets, 86-91t). Preset assignments confirmed correct.

## Notes from .positai/

### Expertise files copied to dev/expertise/

All 6 expertise files copied:
- **coordination.md** — (copy of existing coordination.md reference; kept for legacy context)
- **fitch-scoring.md** — Technical reference on Fitch scoring implementation
- **profiling.md** — R package profiling techniques and tools
- **red-team.md** — Code review and correctness verification checklists
- **shiny-app.md** — Shiny app architecture and development notes
- **tnt.md** — TNT algorithm comparison and benchmarking notes

### Plan files copied to dev/plans/

- **2026-03-22-1348-full-polytomy-search-for-treesearch-c-engine.md** — In-depth design for polytomy-search (collapsed-edge optimization). Approach B chosen (binary internals + collapsed-edge flags, ~16–24 agent-days estimated vs Approach A ~9-13 weeks). No C++ changes needed beyond Phase 1–10 (regions, TBR/SPR/drift, pool dedup, ratchet, sectorial, Wagner, testing, benchmarking). TNT benchmark re-run planned to validate score parity.

### Briefing files reviewed

- **briefing-multistate-profile.md** — T-101 done; T-102–T-107 open. Extends profile parsimony from 2 to multi-state (3+). Reuses MaddisonSlatkin() from concordance-FitchInfo branch for multi-state information content. Recommends MC-calibrated normal approximation for >5 states with exact anchor at s_min. ~3.5-hour implementation effort estimated.
  - **Decision: KEPT.** Contains non-derivable mathematical theory and prototype R code for multi-state profile parsimony. Survival value: guides T-102–T-107 task execution and performance tuning.

- **briefing-progressive-results.md** — T-129. Recommends progress-file polling using existing C++ callback infrastructure (no C++ changes needed; TREESEARCH_PROGRESS_FILE env var). Mirrors cancel-file pattern. ~2–3 hours estimated. Max-rep/best-score/hits display during search.
  - **Decision: KEPT.** Contains implementation guidance and correctness rationale (why NOT to stream partial trees mid-search). Survival value: prevents re-analysis of the rejected alternatives (partial-tree streaming, R-level chunking).

### .positai/settings.json reviewed

**Content:** PositAI-era Sonnet 4.6 model config + permission allowlist (edit *.md/*.h/*.cpp/*.R, bash commands, git, Hamilton-HPC + r-package-profiling skills, TreeDist/TS-MadSlat external dirs).

**Decision:** NOT copied to dev/. Current `.claude/settings.json` supersedes this entirely (Claude Code replaces PositAI). The model ID, thinking effort, and skill references are no longer applicable (Claude Code doesn't use PositAI providers). Permission allowlist is project-specific but `.claude/settings.json` will be maintained as the canonical config.

### .positai/skills/ directory noted

- **hamilton-hpc/SKILL.md** — Hamilton HPC integration skill. Deferred to separate Claude Code skill setup (not copied to dev/). These become `.claude/commands/` or Claude Code integrations separately.

---

## Summary of archival decisions

| Category | Files | Action | Justification |
|----------|-------|--------|----------------|
| **expertise** | 6 files | → dev/expertise/ | Still load-bearing technical references; decoupled from PositAI |
| **plans** | 1 file | → dev/plans/ | Polytomy-search plan (16–24 agent-days) needs full context; referenced in to-do |
| **briefings** | 2 files | → dev/briefings/ | Contain non-derivable theory + implementation guidance for open tasks |
| **settings.json** | — | Discard | Superseded by `.claude/settings.json`; PositAI config no longer applicable |
| **skills/** | 1 file | Note only | Hamilton-HPC → Claude Code skill (separate setup); not duplicated |

---

## Word count: 732 words (this document)
