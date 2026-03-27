# TreeSearch Task Queue

## How this works

- Tasks are sorted by priority (highest first within each status group).
- An agent claims a task by changing its status to `ASSIGNED (X)`.
- When a task is being developed in a **git worktree**, set its status to
  `WORKTREE (name)` where *name* is the worktree directory (e.g.
  `WORKTREE (TS-CID-cons)`). This distinguishes human/long-running worktree
  work from agent assignments and prevents double-claiming.
- On completion, **delete** the row from this file and append a summary row
  to `completed-tasks.md` (see workflow in AGENTS.md).
- Tasks awaiting GHA results: `PARKED (<Letter>, GHA <run_id>)`.
- Tasks with an open PR awaiting human merge: `PR #N (<Letter>)`.
  S-COORD cleans these up after merge.
- Standing tasks (S-RED, S-PROF, S-COORD) are always present. When one is
  completed, reset it to OPEN. Their effective priority is dynamic:
  - ≥6 OPEN specific tasks → standing tasks are P3
  - 3–5 OPEN specific tasks → standing tasks are P2
  - <3 OPEN specific tasks → standing tasks are P1

---

## Active Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-150 | P2 | PR #213 (F) | — | **CID-optimal consensus tree search** | PR #213. Vignette fix (TreeTools::Consensus) commit f8bfee49. GHA 23650002703. |
| T-204 | P2 | PR #216 (F) | — | **Decouple R-loop search from MorphyLib.** Native C++ scorer defaults for `TreeSearch()`, `Ratchet()`, `Jackknife()`; `concavity` param; MorphyLib soft-deprecated. | GHA 23649607006 PASSED. Ready for merge. |


### Bugs

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-279 | P2 | PR #237 (F) | — | **Constrained drift: stale `cd->constraint_node` after RFD rejection + redundant full_rescore.** In `ts_drift.cpp::drift_phase()`, when a constrained suboptimal move passes the constraint check (so `cd` is updated), then is rejected by the RFD test, `update_constraint()` is not called on the restored topology. Bug in both IW reject (lines 647–649) and EW reject (line 689). Same class as T-278. Also eliminates a redundant `drift_full_rescore()` call in the EW reject path (return value at line 657 was discarded, causing a second rescore at line 689). Only affects constrained drift (`driftCycles > 0`); drift disabled in all presets since T-255. | Found by S-RED focus 8 (F, 2026-03-27). PR TBD. feature/drift-constraint-fix commit e85ec84f. |

### Performance Optimization (180+ tips)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-245 | P3 | OPEN | — | **TBR candidate batching.** Restructure TBR rerooting inner loop to evaluate 4 regraft candidates in lockstep, exploiting memory-level parallelism (while one candidate's data transits L2→L1, ALU works on another). Phase profiling shows TBR+enumeration = 86% of 180-tip wall time; estimated ~13% overall gain. | New branch `feature/tbr-batch`. Validate on Hamilton with same benchmark setup. Most invasive change — needs careful correctness testing. |

### Alternative Homologies (Goloboff 2026) — `feature/alt-homology` / `TS-AltHom`

Ref: Goloboff (2026) *Cladistics* doi:10.1111/cla.70033.
Plan: `.positai/plans/2026-03-27-1415-implement-goloboff-2026-alternative-homologies-with-step-matrix-recoding.md`

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-280 | P3 | OPEN | — | **AltHom Phase 1: `AlternativeHomology` S3 class & core recoding (MVP).** Create `R/AlternativeHomology.R` (constructor, validation, print), `R/recode_alt_homology.R` (correspondence enumeration, morphotype states, cost matrix, tip assignment). Wire into `TreeLength()` for scoring on a fixed tree. Reproduce paper's Definition 1 cost matrix + Table 1 as tests. | WORKTREE (TS-AltHom). Invertible, no external constraints, two part-types only. |
| T-281 | P3 | OPEN | T-280 | **AltHom Phase 2: Constraints & options.** Non-invertible (`>`), adjacent (`>>`), restricted homology (`!`), configurable part transformation costs, adjacent-loss merging (`<`). Reproduce Definitions 2–3 and their cost matrices. | WORKTREE (TS-AltHom). |
| T-282 | P3 | OPEN | T-280 | **AltHom Phase 3: Wire into `MaximizeParsimony()` search pipeline.** Accept `AlternativeHomology` in `hierarchy` param, prepare xformArgs, end-to-end search. Also wire `Resample()` and `SuccessiveApproximations()`. | WORKTREE (TS-AltHom). |
| T-283 | P3 | OPEN | T-280 | **AltHom Phase 4: External inapplicability.** An external character can make individual characters, parts, or entire part sets inapplicable. Expand state enumeration for externally-disabled states. | WORKTREE (TS-AltHom). |
| T-284 | P3 | OPEN | T-280 | **AltHom Phase 5: Combination pruning.** Implement `xlinks&` (pairwise compatibility), `xlinks!` (observed-state-only), `xlinks@` (uninformative-state restriction) to reduce supercharacter state count. Verify same optimal trees as unpruned. | WORKTREE (TS-AltHom). |
| T-285 | P3 | OPEN | T-280 | **AltHom Phase 6: Implied weighting support.** Compute combined minimum steps across all valid alignments (not sum of per-char minima). Required for correct IW homoplasy counts. | WORKTREE (TS-AltHom). |
| T-286 | P3 | OPEN | T-280 | **AltHom Phase 7: Mixed `AlternativeHomology` + `CharacterHierarchy`.** Support datasets with both simple hierarchy blocks and alternative homology blocks in one analysis. | WORKTREE (TS-AltHom). |
| T-287 | P3 | OPEN | T-284 | **AltHom Phase 8: Static alignment fallback.** For datasets where supercharacter exceeds practical state limit, generate alternative static datasets (one per alignment) and search each. | WORKTREE (TS-AltHom). |
| T-288 | P3 | OPEN | T-282 | **AltHom Phase 9: Documentation & vignette.** `vignettes/alternative-homologies.Rmd`, roxygen docs for all new exports, `inst/REFERENCES.bib` entry. | WORKTREE (TS-AltHom). |

### TNT Comparison & Strategy Learning

### Strategy Tuning

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-289 | P2 | PARKED (F, SLURM 16606222) | — | **Prune-reinsert benchmark (Hamilton).** Stage 1: sweep cycles (1/3/5) × drop fraction (5/10/20/30%) × RANDOM on 5 datasets (37–180t), 5 seeds, 30s. Stage 2: best configs + INSTABILITY selection + ratchet-replacement at 30s/60s. Scripts: `bench_prune_reinsert.R`, `t289_hamilton.sh`. Goal: decide whether to enable in presets and at what settings. | Fixed Rscript invocation bug (R_LIBS_USER pattern). SLURM 16606222 submitted 2026-03-27. |
| T-269 | P3 | ASSIGNED (F) | — | **Fine-grained sectorial interleaving benchmark.** Compare current coarse `outerCycles=2` (all ratchet cycles batched per pass) against fine-grained interleaving (e.g. `outerCycles=12, ratchetCycles=12` → one ratchet cycle per sectorial pass), approximating TNT's per-iteration pattern. T-256 showed extra sectorial *rounds* don't help, but *timing* of sectorial relative to perturbation hasn't been tested. **Design note (u.005):** Full TNT-style interleaving IS architecturally supported now — setting `outerCycles = ratchetCycles` achieves one sector pass per ratchet cycle. T-257 first validated that any post-ratchet sectorial helps at all (merged). T-269 benchmarks the fine-grained variant to determine whether the per-cycle overhead is worth enabling in presets. | Low priority. Use 3–5 gap datasets at 30s/60s budgets. |


### Housekeeping

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|





### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | OPEN | — | **Standing: Red-team review** | Last run: 2026-03-27 focus 16 by F (ts_ratchet.cpp, 259+61 lines). Bug found and fixed: F-015 (stale constraint_node after best_tree revert in ratchet_search() non-escape path — commit ae6a3528). GHA 23653228247. Same class as T-278/T-279/E-003. All other invariants correct (save/restore, FlatBlock sync, perturb modes, adaptive tuning). Next: ts_resample.cpp or ts_nni_perturb.cpp. |
| S-PROF | dyn | OPEN | — | **Standing: Performance profiling** | Last run: 2026-03-27 by A (round 6: thorough-preset phase distribution at 75t; NNI-perturb 34% time / 14% hit rate; T-274 filed). |
| S-COORD | dyn | OPEN | — | **Standing: Coordination review** | Last run: 2026-03-27 round 37 by F. T-289 dispatched (SLURM 16606222, ~2.7h). F-015 fixed (ratchet constraint staleness, GHA 23653228247). PRs #213/#216/#237 still open awaiting merge. OPEN specific (non-parked): T-245, T-269 → standing at P2 (with T-289 parked, 2 OPEN + 1 PARKED = effective ≥3 → P2). |
| S-PR | dyn | OPEN | — | **Standing: PR maintenance** | Last run: 2026-03-27 round 40 by F. All 3 ready PRs unmerged: #213 (T-150, GHA 23650002703 PASS), #216 (T-204, GHA 23649607006 PASS), #237 (T-279, GHA 23650290962 PASS). cpp-search F-015 ratchet fix: GHA 23653228247 running. No new PRs needed (F-015 is direct cpp-search fix). Open PRs: #213 (T-150), #216 (T-204), #237 (T-279), #210 (DRAFT). |
