# TreeSearch Task Queue

## How this works

- Tasks are sorted by priority (highest first within each status group).
- An agent claims a task by changing its status to `ASSIGNED (d1)` (or `d2`,
  `d3`, … — ephemeral dispatcher IDs issued by `dispatch.sh`).
- When a task is being developed in a **git worktree**, set its status to
  `WORKTREE (name)` where *name* is the worktree directory (e.g.
  `WORKTREE (TS-CID-cons)`). This distinguishes human/long-running worktree
  work from agent assignments and prevents double-claiming.
- On completion, **delete** the row from this file and append a summary row
  to `completed-tasks.md` (see workflow in AGENTS.md).
- Tasks awaiting GHA results: `PARKED (d1, GHA <run_id>)`.
- Tasks with an open PR awaiting human merge: `PR #N (d1)`.
  S-COORD cleans these up after merge.
- The `Notes` column may include a bracketed model/effort hint, e.g.
  `[m:haiku e:low]`, `[m:sonnet e:medium]`, `[m:opus e:high]`. The
  dispatcher's ranker honours these hints and they override its default choice.
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


### Shiny App

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|


### Alternative Homologies (Goloboff 2026) — `feature/alt-homology` / `TS-AltHom`

Ref: Goloboff (2026) *Cladistics* doi:10.1111/cla.70033.
Plan: `dev/plans/2026-03-27-1415-implement-goloboff-2026-alternative-homologies-with-step-matrix-recoding.md`

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

### Deferred / Future Directions

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-290 | — | DEFERRED | — | **GPU-accelerated batch tree scoring.** Evaluate many TBR/SPR candidate rearrangements in a single GPU kernel launch (parallelism across *trees*, not within one tree). For a 180-leaf tree the TBR neighborhood is O(n³) ≈ millions of candidates — enough to saturate GPU hardware. Main challenges: (1) per-candidate work is tiny for Fitch+bitwise (~50 word ops), so GPU arithmetic intensity is very low; (2) tree data structures need flat-array redesign for coalesced GPU memory access; (3) for morphological data sizes (≤500 chars, k ≤ 10) CPU OpenMP parallelism across candidates likely captures most of the win with far less effort. GPU becomes more compelling for Sankoff/implied-weights scoring (O(k²) per node per char) or phylogenomic-scale data (10k+ chars). A hybrid design (CPU manages search logic, GPU batch-scores candidates) is more practical than porting the full search engine to CUDA. **References:** Santander-Jiménez et al. (2020) *J Supercomput* 76:9827 (GPU Fitch parsimony, Kepler→Turing); Santander-Jiménez & Vega-Rodríguez (2025) *Future Gen Comput Syst* (OpenMP/OpenACC/SYCL multi-platform parsimony scoring); Ayres et al. (2019) *Syst Biol* 68:1052 (BEAGLE 3 — GPU likelihood, architectural lessons). | Research: MkPrime `.agent-d.md` 2026-03-29. |
| T-291 | — | DEFERRED | — | **GPU-parallel independent search replicates.** Run 100+ search replicates simultaneously on GPU SMs (one replicate per SM; modern GPUs have 60–128 SMs). Shared read-only character matrix fits in GPU L2 cache. Main obstacle: tree search has highly irregular, data-dependent control flow (rearrangement selection, acceptance decisions, ratchet perturbation) which causes warp divergence and poor GPU utilization. Branch-and-bound in sectorial search has the same problem. CPU multicore parallelism (8–16 cores via `future`/`parallel::mclapply`, or 100+ via HPC SLURM array jobs) is far simpler and more efficient per-replicate. GPU replicates only become attractive if per-replicate arithmetic is heavy enough to dominate over control flow overhead (e.g., large Sankoff matrices). **References:** same as T-290. | Research: MkPrime `.agent-d.md` 2026-03-29. |

### TNT Comparison & Strategy Learning

### Strategy Tuning


### Housekeeping

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| T-298 | P3 | PR #242 (d8) | — | **Profile and optimize `quartet_concordance.cpp` matrix allocation** | GHA 25777319791. Resize-hoist committed, benchmarked, PR open. |
| T-300 | P3 | PARKED (d7, other ) | — | **Lazy `apply_tbr_move` rescore in `tbr_search`.** After the `score_fresh` flag was wired (companion to T-187/PSF work), the trailing `full_rescore` at function exit is now skipped when states are coherent. The remaining redundancy is the `full_rescore(tree, ds)` call at `ts_tbr.cpp:1134`, run after **every** successful `apply_tbr_move` to obtain the authoritative score for the acceptance check. Each call is O(n_node × total_words). Since the move is local (clip + reroot + regraft), the indirect-evaluation pre-check at `ts_tbr.cpp:767-772` already shows that `fitch_incremental_downpass/uppass` from the join point gives the correct score in O(affected_subtree_depth × total_words) instead. **Plan:** make `apply_tbr_move` push touched nodes onto the prealloc_undo stack, return the join node, and replace line 1134's `full_rescore` with `fitch_incremental_downpass` from that node. Estimated savings: O(n_char) per accepted move × ~10–100 accepted moves per replicate. **Risk:** medium — `apply_tbr_move` is the hot correctness-critical path; need careful unit tests covering NA/non-NA, IW/EW, constrained/unconstrained, equal-accept paths. Validate by comparing scores against current unconditional rescore on a battery of datasets before committing. |






### Standing Tasks

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| S-RED | dyn | PARKED (d4, other ) | — | **Standing: Red-team review** | Last run: 2026-03-28 focus 31 by G. ts_prune_reinsert.h/.cpp (583 lines): G-006 found + now fixed. Next: ts_search.cpp (NNI/SPR, 421 lines) and ts_nni_perturb.h/.cpp (unreviewed). |
| S-PROF | dyn | OPEN | — | **Standing: Performance profiling** | Last run: 2026-05-12 round 7 by d6. T-260 hotspot audit: std::fill (9.1%) fixed T-261 ✓; StateSnapshot per-candidate save (14.6%) mitigated by opt #7 (once-per-pass) ✓; full_rescore line 1137 (~28%) → T-300 in progress (d7). No new tasks; re-profile with VTune after T-300 lands. |
| S-COORD | dyn | PARKED (d6, other ) | — | **Standing: Coordination review** | Last run: 2026-05-12 round 47 by d5. u.118 triaged → T-301 (progress ticker multi-thread). PR #210 (cpp-search→main): 2 CI failures are infra (ASAN vignettes: missing pkgdown; Windows: code-coverage only — R CMD check passes). PR #213 (T-150): CONFLICTING, no recent CI — needs human to resolve merge conflicts. PR #216 (T-204): agent-check 23649607006 PASSED; full R-CMD-check had failures Mar 2026 on ASAN/Windows/ubuntu-old — needs re-trig or human review. Active: d1(T-294), d2(T-298), d3(T-299), d4(S-RED parked). T-280–288 WORKTREE awaiting. S-PROF/S-PR OPEN. |
| S-PR | dyn | OPEN | — | **Standing: PR maintenance** | Last run: 2026-05-12 round 48 by d5. PR #210 (cpp-search→main): MERGEABLE, fresh checks 2026-05-12 confirm 2 infra-only failures (Windows=code-coverage step, ASAN=vignettes infra) — all R-CMD-check PASS; ready for human to un-DRAFT and merge. PR #216 (T-204, feature/native-search→cpp-search): CONFLICTING — cpp-search has ~10 new commits since last merge; needs rebase then re-trig. PR #213 (T-150, feature/cid-consensus→cpp-search): CONFLICTING, no CI — needs rebase onto cpp-search. |
