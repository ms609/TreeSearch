# Fuse / Drift component isolation (task #52, component-isolation slot 4)

**Mission:** close TreeSearch's per-iteration wall gap to TNT 1.6 on equal-weights
(EW) Fitch. Quality (gapB) closed; throughput residual ~1.3-2.3× concentrated on
heavy multistate (Zanol/Zhu). NA/IW out of scope. Fitch via `-`→`?`.

**Where this sits:** TBR is CLOSED on every gate + every algorithmic thread (kernel
at-limit T-P5l / precompute dead M46/T-P5j / scaffolding below-floor T-P5m /
middle-level at-best T-P5p / lever-c dead-by-proof T-P5q). Scoring CLOSED (T-P5l).
Sectorial gate-1 done by the sectorial agent (T-S6a-e: ~96% of isolated sectorial
wall IS the inner `tbr_search` ⇒ at-limit-by-inheritance; ~2.8% byte-identical
micro-levers banked). Ratchet recipe banked (12→6, T-P5d). **Fuse + Drift are the
last untouched component slot** (user-ordered after lever-c, 2026-06-20). Composition
#40 is GATED LAST (user: recompose-from-scratch on any step-cost change ⇒ finish all
pieces first).

## Two gates per component (the program's contract)
1. **AT-LIMIT?** isolate the component's hot path; is the fuse/drift-SPECIFIC work at
   the AVX2/compiler/bandwidth ceiling, or is there a real lever? (chrono decomposition,
   like the sectorial agent's `TS_SECT_TIMING`.)
2. **vs TNT per-iteration?** shared-start race: score@opt, count-examined, wall
   (32-bit local = lower bound; Hamilton 64-bit = authoritative).

## Recon findings (2026-06-20, this session)
- **No `getenv` / hot-path CRT in ts_fuse.cpp or ts_drift.cpp** (Grep). The
  getenv-class win (13-26% TBR T-P5n; ~22% sectorial T-S6d) is ABSENT here — clean.
- **Wiring (ts_driven.cpp):** drift = outer-cycle step 5 (`drift_search`, gated
  `drift_per>0`, ts_driven.cpp:435-445) followed by a TBR polish (:544). Fuse = two
  paths: intra-replicate `tree_fuse(result.tree, ds, *pool, fp)` gated
  `params.intra_fuse && pool->size()>=1` (:568-587) + inter-replicate every
  `fuse_interval` (:954-981). `tree_fuse` runs `tbr_search` internally after each
  improvement round (the at-limit kernel) ⇒ STRONG at-limit-by-inheritance prior.
- **Default-share is SMALL:** T-P5c phase table had fuse 0-2.5%, drift not prominent;
  both are opt-in / preset-specific, NOT default-mission-wall hogs. So the gate-1
  prize is byte-identical micro-levers + the crash fix, not a big throughput lever.
- **CRASH-FIX PREREQUISITE (verified still ABSENT on cpp-search):** `intraFuse=TRUE`
  SEGFAULTS on >64-tip data (Zanol/Zhu/Giles — the heavy mission datasets).
  `reroot_at_tip0(recipient)` is called ONCE pre-loop (ts_fuse.cpp:332), NOT per round
  (round loop :357 has no re-root); `replace_subtree` (:221) has no
  `r_rest.size()!=d_rest.size()` guard. Root cause + tested fix in worktree
  `TreeSearch-nonclade` (feature/nonclade-sectors): per-round re-root + size guard,
  9/9 fuse tests, regression test at 80 tips. See MEMORY [[fuse-reroot-segfault]].
  ⇒ Fuse can only RUN on ≤64t (Wortley 37t) today; gate-2 on heavy data + fuse being a
  real mission lever both REQUIRE this port.

## Plan (4 phases)
- **P1 — gate-1 code analysis (no build):** characterize fuse/drift-specific vs
  inherited-TBR split; hunt byte-identical micro-levers (allocs-in-loops, redundant
  rebuilds, no-op round-trips à la T-S6c sectorial `ras_starts==1`); verify the
  crash-fix soundness + produce a port spec. [workflow]
- **P2 — crash-fix port:** apply per-round re-root + `replace_subtree` size guard on a
  fresh worktree off cpp-search (named files only, NO shared-branch commit); build
  per-agent; run `test-ts-fuse.R` + a >64t no-crash repro (Zanol intraFuse reps2).
- **P3 — empirical gate-1:** chrono decomposition (gated instrument) on small data
  (Wortley ≤64t for fuse; drift wherever it runs) confirming the at-limit-by-
  inheritance split + any micro-lever wall delta (local iterate-tier, ~seconds).
- **P4 — gate-2 race:** shared-start vs TNT `tfuse` / drift flags → Hamilton 64-bit
  (needs P2). score@opt + count + wall as seed distributions.

## Guardrails (unattended)
Per-agent `R CMD INSTALL` only (never load_all for perf/correctness); builds + ~30s
targeted tests stay LOCAL; heavy/parallel + 64-bit races → Hamilton SLURM (durable
per-cell output); code changes on a WORKTREE, stage named files only, NO commit unless
asked, never touch cpp-search/main broadly ([[concurrent-session-git-hazard]]); no
thread_local hot-path scratch (MinGW emutls); UK 'ize'; `return x;` C++; TreeTools over
ape. ts_sector.cpp belongs to the sectorial agent — do not edit.

## PROGRESS LOG
<!-- newest at bottom; ### <ISO> — <phase> — <one-line state> -->

### 2026-06-20 — P0 — Oriented; recon done; gate-1 analysis workflow launched; crash fix confirmed absent
Recon above. Next: P1 analysis workflow + P2 crash-fix port worktree.

### 2026-06-20 — P2 DONE — crash fix ported + validated + committed (worktree branch)
Diffed cpp-search ts_fuse.cpp vs nonclade: ONLY the 2-hunk fix differs (per-round
`reroot_at_tip0` after `++result.n_rounds`; `replace_subtree` size guard) — confirms
ts_fuse.cpp otherwise byte-identical. test-ts-fuse.R diff = purely additive (the
80-tip >64 regression test). Ported on isolated worktree
`C:/Users/pjjg18/GitHub/worktrees/TreeSearch/fuse-reroot-port` (branch
`claude/fuse-reroot-port` off cpp-search 8c57c2ec); built per-agent (`.agent-fuse`,
ccache, exit 0); **22/22 fuse tests PASS incl. the 80-tip regression** (NOT_CRAN=true).
Committed `da21f5dc` (2 named files, NOT pushed, NOT on cpp-search). **The fix is a
real CORRECTNESS fix (intraFuse segfaults on all >64t mission data) and should be
merged to cpp-search by the supervisor** (cherry-pick da21f5dc or the 2 hunks).

### 2026-06-20 — P1 (drift independent read) — drift is at-limit-by-inheritance, CONFIRMED by structure
Read ts_drift.cpp end-to-end: it is the TBR kernel DUPLICATED — `drift_collect_main_edges`
/`drift_collect_subtree_edges`/`drift_fitch_join_states`/`drift_compute_from_above`/
`drift_apply_tbr_move` are all "mirrored from ts_tbr.cpp", and `drift_phase` is a
tbr_search inner loop with an AFD/RFD accept rule; `drift_search` (cycles) calls
`tbr_search` DIRECTLY for the equal-score (:760) and convergence (:772) phases. So
drift = TBR + accept-rule ⇒ throughput rides the now-closed TBR kernel ⇒
**at-limit-by-inheritance** (same verdict as sectorial T-S6a). Drift-SPECIFIC code =
the AFD/RFD accept logic + `drift_full_rescore` (full O(N) rescore on accept/decision,
bounded by ≤max_drift_changes accepts) — no getenv, no obvious non-inherited lever.
NB the drift_* duplication is a MAINTAINABILITY smell (copies of TBR code), not a perf
lever. Awaiting P1 workflow for the byte-identical micro-lever hunt + fuse split.

### 2026-06-20 — P1 DONE + a STANDOUT QUALITY FINDING (bigger than the whole perf surface)
Gate-1 workflow verdict: **fuse + drift = AT-LIMIT-BY-INHERITANCE** (~85-90% inherited
kernel by structure) ⇒ gate-1 now CLOSED across ALL FOUR components (scoring, TBR,
sectorial, fuse/drift). Byte-identical micro-levers found are all sub-0.1% mission
(below sectorial's 2.8%): replace_subtree unordered_map→flat-vector, new_local_cost
hoist, drift_compute_from_above scratch hoist — bank as hygiene, not the point. Tier-2
(build_postorder_prealloc + elide-triple-rescore) = byte-identical-BY-ARGUMENT only
(touches the RFD-accept local_cost MASK) ⇒ DEFER (sub-floor, needs mask-equality A/B).

**THE FINDING (confirmed by code, NOT a perf lever — a SEARCH-QUALITY defect):** the
perturbation/secondary engines score candidate moves with `fitch_indirect_length_bounded`
& friends = the **union-of-finals (`final_[A]|final_[D]`) approximation that UNDERCOUNTS**
(ts_fitch.h:118-126 explicitly: edge_set directional `_cached` is "the CORRECT
replacement ... which undercounts"; ts_tbr.cpp:1608-1612 "Exact directional cost ...
replacing the union-of-finals approximation"). The directional fix (a PROVEN real
quality bug: oracle 23/40→9/60 [[tbr-rooted-vs-unrooted]]; Wagner +30%
[[wagner-insertion-cost-bug]]) was ported to the 3 MAIN kernels — `ts_tbr.cpp`(EW/IW
`_cached`), `ts_wagner.cpp`, `ts_sector.cpp`(#27) — but **NOT** to the secondary engines:
- `ts_drift.cpp` :492/499/547/555 (EW+IW) — drift candidate ranking + AFD-gate.
- `ts_prune_reinsert.cpp` :439/448 (EW).
- `ts_search.cpp` :358/365 (EW/IW).
- `ts_temper.cpp` :291 (NA+IW — NA out of scope).
- (tbr_search :1589/1596 NA path also `_bounded`, but NA is another agent's scope.)
Undercounting mis-ranks candidates ⇒ these engines perturb toward mis-scored targets ⇒
DEGRADED escape efficacy. Mission relevance hinges on (a) which are LIVE on the default
MaximizeParsimony path + their share, (b) magnitude of the discrepancy, (c) oversight vs
deliberate cheap-approx tradeoff (fixing = drift must pay the directional precompute the
way tbr_search does). Launched workflow `fuse-drift-scoring-audit` to settle these. If
real+live+impactful: scope+implement the directional fix per path on a worktree (mirror
tbr_search / #27), local-validate, queue a Hamilton QUALITY A/B; do NOT land unvalidated
(trajectory change). This is the anti-satisficing lead ([[tnt-outperformance-is-diagnostic]]).

### 2026-06-20 — SCOREAPPROX audit DONE — REAL-BUT-OFF-DEFAULT-PATH (not a recorded-quality bug)
Workflow `perturbation-scoring-audit` (3 lenses + synth). VERDICT: the discrepancy is
REAL (the `_bounded` union-of-finals = the proven Wagner-+30% / TBR-oracle-23→9
undercount; it mis-ranks candidates AND shifts the drift AFD/RFD accept band), BUT it
is **NOT a recorded-quality bug on any path**. DECISIVE RECONCILIATION (code-read):
every secondary engine uses `_bounded` ONLY to RANK a perturbation/start, then
EXACT-reconverges via `tbr_search`/`nni` (`_cached`) and KEEPS only on STRICT
improvement (prune ts_prune_reinsert.cpp:558-580; anneal ts_driven.cpp:482-486; drift
ts_drift.cpp:760-777) ⇒ a mis-rank only WASTES A CYCLE; the recorded optimum is always
produced by the exact kernels. So gapB=0/efficiency≈1 (measured ~70t) is NOT
contradicted, and this is **NOT the throughput gap's cause**.
- **Live-path map:** drift=opt-in (0% default); **prune_reinsert=preset-only, AUTO in
  `large` (≥120t), cycles=5L** (MP.R:224) — the ONLY auto exposure; temper=preset-only
  (large, annealCycles=1L, but a DEFENSIBLE tradeoff — 1 edge/step, non-amortizable —
  DO NOT convert); spr_search=opt-in (sprFirst=FALSE everywhere). gapB=0 was never
  established for the large preset (≥120t) ⇒ prune_reinsert there is "live but untested".
- **EXCEPTION (honest):** spr_search (ts_search.cpp:365) is a HILL-CLIMBER w/
  single-best verify ⇒ a mis-ranked-away improver is SILENTLY MISSED (the real harm
  mechanism) — but OFF every preset; ownership uncertain (confirm vs legacy before edit).
- **DOC BUG:** ts_fitch.cpp:385-391 comment "union exact; intersection overcounts" is
  BACKWARDS (authoritative header ts_fitch.h:118-126 = union UNDERCOUNTS). Byte-identical
  one-line fix; worth landing to stop the misconception re-spawning.
- **FOLLOW-UP (task #53, gated, throughput-NEGATIVE so NOT assumed-good):** (1) fix the
  backwards comment; (2) cheap local FLIP-PROBE — enable pruneReinsertCycles on a small
  dataset, compute `_cached(edge_set[D])` alongside `_bounded` at ts_prune_reinsert.cpp:
  439/448, log value-disagreements + argmin/best-edge flips; if flips≈0 → CLOSE with no
  port; (3) only if material flips → port to `_cached` (mirror #27 build_ras_sector:
  one compute_insertion_edge_sets per dropped tip, reused over DFS edges) on a worktree +
  time-matched Hamilton A/B (≥120t, ≥10 seeds; ship ONLY if neutral-to-better, since it
  adds the ~30%-EW precompute). drift/spr_search ports = defer to opt-in/human.
- **SCOPE GUARDS:** NA scoring paths OUT OF SCOPE (other agents). No measured quality
  loss on ANY path ⇒ do NOT over-claim; not the throughput gap. anneal stays `_bounded`.

### 2026-06-20 — #52 DISPOSITION — fuse/drift gate-1 CLOSED; gate-2 = Hamilton-confirmatory
Gate-1 (AT-LIMIT): CLOSED for fuse + drift (at-limit-by-inheritance) — completes gate-1
across ALL FOUR components. Crash fix landed (worktree da21f5dc, flagged for merge).
Byte-identical micro-levers found but all sub-0.1% mission (bank as hygiene, optional).
Gate-2 (TNT race): for fuse it is intrinsically awkward (fuse needs a diverse POOL, not
a single shared start) and throughput is inherited from the closed TBR kernel ⇒ a race is
confirmatory; the meaningful fuse/drift question is RECIPE value = composition #40 (gated).
Recommend: gate-2 race for fuse/drift is LOW-priority Hamilton-confirmatory, NOT a blocker.

### 2026-06-20 — SCOREAPPROX ELEVATION (read ts_prune_reinsert.cpp:412-468) + RATCHET gate-1
**prune_reinsert is STRONGER than the synth's "wasted cycles":** `expand_and_reinsert`
does INCREMENTAL GREEDY WAGNER reconstruction (wagner_incremental_rescore per tip, :467),
scoring each candidate edge with `_bounded` (:439/448) on a PARTIAL tree = the
CONSTRUCTION regime where the union undercount was measured at +30% (the original Wagner
bug), and it is the VERBATIM greedy-placement pattern Wagner/build_ras_sector (#27) were
fixed for — the sibling was missed. Strict gate still protects the RECORDED score, but
hampered reconstruction ⇒ prune_reinsert escapes LESS effectively ⇒ the large preset
(≥120t, the only auto path) may reach worse optima ⇒ a likely real LARGE-TREE EFFICACY
loss, not just wasted cycles. Flip-probe subtlety flagged in #53 (exact `_cached` needs a
current `prelim` downpass; incremental Wagner maintains `final_`). Recorded in #53 (do NOT
rush a probe that could give false flips). **RATCHET gate-1:** recon (ts_ratchet.cpp) =
NO getenv; ratchet = `perturb_upweight` (cheap O(chars) reweight) + `tbr_search`
(:153/203/209, the closed kernel) ⇒ AT-LIMIT-BY-INHERITANCE, recipe banked (12→6 T-P5d).
⇒ **gate-1 (AT-LIMIT) now COMPLETE across ALL components** (scoring/TBR/sectorial/fuse/
drift/ratchet). Remaining isolation work = gate-2 TNT races (Hamilton-confirmatory:
sectorial=other agent, ratchet+fuse/drift low-priority) + #53 + composition #40 (gated).

### 2026-06-20 — #53 RESOLUTION — backwards comment fixed; prune_reinsert Δ-probe DONE; port PREPARED; A/B composition-gated
**(1) DOC BUG FIXED + LANDED on cpp-search** (8671fdaa): ts_fitch.cpp:385-391 backwards
comment ("union exact; intersect overcounts") corrected to match header :118-126 (union
UNDER-counts; the directional edge_set is exact). Doc-only, mission-safe.

**(2) Δ-PROBE (not flip-count — advisor: tally exact-suboptimality, not edge-identity):**
gated `-DTS_SCOREAPPROX_PROBE` in expand_and_reinsert, non-perturbing (production still
inserted at the `_bounded` choice). Per placement tallied Δ = exact_cost(E_bounded) −
min_E exact_cost(E), exact scorer = compute_insertion_edge_sets + fitch_indirect_length_cached,
no cutoff. `prelim` confirmed current (wagner_incremental_rescore maintains it) + in-tree
fully binary from root ⇒ precompute safe. Result on Zanol (forced pruneReinsertCycles):
**~62% of placements strictly worse, mean ~6 steps, max 37, ~48% greedy-regret SHARE**
(bounded_exact_sum 6406 vs min_exact_sum 4315). Corroborates the validated +30% Wagner bug.

**(3) PORT PREPARED + VALIDATED (worktree claude/scoreapprox-probe, 41b0d237; NOT cpp-search):**
swapped the two `_bounded` calls for the exact `_cached`+edge_set (mirror ts_wagner.cpp:487).
After port the probe reports **Δ=0 at every placement** (production == exact argmin). Tests:
prune-reinsert 44/0, drift 22/0, ratchet 17/0, tbr 28/0.

**(4) PATH-RELEVANCE KILL for the heavy A/B (the decisive gate):** prune_reinsert auto-enables
ONLY at nTip≥120 (`large` preset, MP.R:249). **NO mission dataset reaches 120t** — full
inapplicable.phyData roster max = Dikow2009 88t; Zhu2013 75t, Zanol 74t, Giles 78t (all
`thorough` or smaller). So this path runs on ZERO default mission searches ⇒ the 48% is
greedy-regret SHARE on a config the mission suite never triggers, NOT a wall-clock
opportunity. blame: `_bounded` = afbf531f (2026-03-27, original T-266) PREDATES the June
directional fix ⇒ a genuine MISS, not a deliberate large-N tradeoff. `large` polish is NNI
(weaker than TBR) ⇒ regret survives more ⇒ fix WOULD matter at ≥120t.

**(5) DISPOSITION:** land + time-matched A/B (needs a ≥120t dataset + the `large`-preset
budget tradeoff, since the exact scorer ADDS an O(N·blocks·9-states) precompute the bounded
path skips) = **COMPOSITION #40** (user: composition waits until all pieces finished). Port
is ready + cost-characterizable for that phase. Component made best-known-correct in
isolation; the enable/wall-clock decision is recipe-level. #53 investigation CLOSED.

**(6) spr_search loose-end RESOLVED (the T-F1 "could silently miss" exception):**
ts_search.cpp `spr_search` (the Fitch SPR, :197) uses the bounded scorer (:365) BUT (a)
fires ONLY when sprFirst=TRUE — FALSE in every preset (off the default path); (b) accepts
ONLY on EXACT `full_rescore` improvement (:388-402) ⇒ can never false-accept, recorded
score always exact (gapB=0 preserved); (c) is a one-shot SPR WARMUP immediately followed by
exact `tbr_search` (ts_driven.cpp `if(!nni_wagner && spr_first){spr_search;} ... tbr_search`)
which re-explores and catches any improver a bounded mis-rank missed. ⇒ the silent-miss is
real-in-principle but MOOTED; porting adds the per-clip precompute for ~zero benefit. NO
ACTION. Remaining bounded sites all benign: drift (opt-in, rank-then-reconverge, T-F1),
temper (preset-only defensible tradeoff, T-F1), ts_rcpp.cpp:2339 (standalone export, not
the recipe). **Scoring-approximation sweep now COMPLETE across the whole search.**

### 2026-06-20 — PHASE-0 CONNECTIVE TISSUE — CLOSED, no addressable production fat
Read of the driven-search orchestration loop (ts_driven.cpp). Full `score_tree`
(O(N·chars)) call inventory at the DEFAULT verbosity (`verbosity=1L`, MP.R:505):
- **All per-phase score prints are `verbosity>=2`-gated** (XSS/RSS/CSS/ratchet/post-sect/
  NNI/drift/SA/PruneRI/TBR/fuse, ts_driven.cpp:249-588) ⇒ DO NOT fire at default v=1.
- **Interrupt/timeout exit branches** (257/269/299/307/431/505/539/563) ⇒ run once on exit.
- **Per-outer-cycle, un-gated:** `score_before_cycle` (:224) + `score_after_cycle` (:594)
  for the convergence/reset check = 2 full rescores/cycle; `score_before_cycle`(N+1) ≡
  `score_after_cycle`(N) (tree unmodified between :594 and next :224) ⇒ one is REDUNDANT.
- **Final:** `result.score = score_tree` (:617) once per replicate.
A full score_tree on Zanol ≈ O(74·210·9) ≈ 140K ops ≈ µs; ~1–few outer cycles/replicate
(outerCycles=1 in `large`) ⇒ total ≈ **0.001% of wall** (seconds of phase work dominate;
score_tree was NOT in the T-P5o hotspot list — consistent). **Step-switching:** each phase
owns/maintains its own prelim/final_ incrementally; the only orchestrator-level state
rebuild is intra-fuse `build_postorder()+reset_states()` (:581-582, preset-only, 1/cycle).
R/C marshalling already T-P5o'd (R.dll 12% = amortizable GC/glue + one-time LoadLibraryA,
startup-inflated by the tiny profiling workload). **VERDICT: Phase-0 AT-LIMIT** — the one
redundant `score_before_cycle` is a sub-floor (~0.001%) bit-identical micro-bank, NOT worth
the convergence-logic risk. This closes the last undone NON-GATED, non-other-agent aspect of
the component-isolation plan. Remaining: gate-2 races (Hamilton-confirmatory; sectorial=other
agent) + composition #40 (gated, where the addressable wall now lives: orchestration / T-S6e).

### 2026-06-20 — bit-packing reopen CLOSED + cherry-pick build-check + Hamilton-KPI BLOCKED
- **ns=9 representation/bit-packing reopen CLOSED analytically (T-P5r, advisor-gated, no build):**
  transposed bitset already bit-dense (9 state-words × 64 patterns = 0.14 op/pattern, 4 states/
  AVX2 instr); states-per-word packing SERIALIZES patterns → strictly worse at 210 patterns;
  the scalar/representation reopen is **ns≤4 only**, deader at ns=9 ⇒ residual ~2× heavy-
  multistate is a genuine ACCEPTED CONSTANT FACTOR, no representation lever. (The T-P5p
  "UNPINNED" tag was the tell — a 21-agent audit had found no concrete scheme.)
- **Cherry-pick build-check PASSED:** clean detached-worktree build of cpp-search HEAD
  (ac8e808a fuse fix + 8671fdaa comment) = INSTALL exit 0; fuse 22/0, tbr 28/0, prune-reinsert
  44/0. No stale-object ABI issue ([[stale-object-abi-gotcha]] cleared). Shared branch safe.
- **Hamilton mission-KPI re-measurement (advisor's highest-value non-gated item) — BLOCKED,
  FLAGGED FOR USER:** the stale TS-vs-TNT wall gap is worth refreshing (predates getenv ~20-26%
  + ratchet 12→6 ~20-38%, which shifted the phase mix). BUT a clean dispatch is blocked: the
  **ratchet 12→6 flip is UNCOMMITTED in the shared working tree** (R/SearchControl.R wt=`6L`;
  origin/cpp-search AND local HEAD both =`12L`), alongside `M` R/MaximizeParsimony.R +
  R/RcppExports.R — another session's in-flight work I must not touch/commit (concurrent-git-
  hazard) and not authorized to push. Cloning origin → measures stale ratchet=12; transferring
  the wt → bundles unowned multi-session WIP. ⇒ cannot define a clean reproducible code-state
  unattended. NEEDS USER: commit the ratchet flip (it's a major banked lever sitting only in the
  working tree — at risk of loss on any `git checkout -- .`) + authorize the Hamilton run.
