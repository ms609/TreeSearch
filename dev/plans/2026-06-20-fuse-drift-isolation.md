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
