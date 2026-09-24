# TBR per-candidate micro-lever banking sweep (task #48)

**Branch:** `claude/tbr-microlevers` off cpp-search `da0f203f`.
**Mandate:** unsupervised, aggressive — bank ANY exact bit-identical non-regressing
win in the TBR per-candidate path. Supervisor decides final merge to cpp-search.
**Autonomy:** commit to THIS branch only; no push; no merge to cpp-search/main;
stage named files only; per-agent installs; single-thread measurement.

## Method
- **Gate (correctness):** `verify_l1.R` — score + candidates_evaluated must be
  BYTE-IDENTICAL to clean da0f203f on Wortley2006/Zhu2013/Zanol2014 x seed{1,2}.
  Documented values confirmed reproduced by a clean-da0f203f BASE build.
- **Wall (magnitude):** `ab_wall.sh` + `analyze_ab.R` — paired SAME-SEED base vs
  mod; identical (byte-identical) work per seed, so the per-seed wall difference
  isolates the per-candidate cost change. Interleaved, order-alternated, sign-test.
- Libs built per-agent (`.agent-mod` / `.agent-base`), RELEASE (no load_all).

## Levers

### L-cutoff-hoist — COMMITTED 6295c401, GATE PASS, wall measuring
Maintain the bail cutoff `(best_candidate - divided_length + 1)` across the clip,
recompute only inside accept blocks instead of per candidate. divided_length is
clip-constant; best_candidate changes only on improvement -> cutoff VALUE at every
scorer call is unchanged. Covers SPR (1578) + scalar-reroot (1798/1816, the
RATCHET hot path) + flat cutoff_b (1697); IW untouched (uses best_candidate
directly). **Gate: 6/6 bit-identical (base==mod==documented).**
**Wall A/B (Zanol, 16 same-seed pairs, reps6): INCONCLUSIVE — noise-swamped.**
Per-pair rel delta ranged -16.5%..+11.9% (machine loaded by the concurrent
discovery workflow); paired median -0.15%, mean -0.47%, sign-test p=1.0. The
~0.15-0.5% effect is ~30x below the wall noise floor. `cand` identical within
every seed-pair (re-confirms byte-identical work; seeds 6/11 score 1262 but
base==mod). **VERDICT: BANK on merits** — gate-proven byte-identical + strictly
fewer per-candidate instructions (cannot regress); wall win real but
sub-noise-floor. NOT claimed as a measured speedup.

## Measurement-strategy pivot (IMPORTANT)
End-to-end wall is the WRONG instrument for sub-1% levers (noise +/-2% quiet,
+/-16% loaded). For the rest of the sweep:
- Bank each EXACT lever individually on gate (byte-identical) merits.
- Measure the CUMULATIVE bundle (base vs all-levers) with a LOW-NOISE IN-PROCESS
  instrument: env-gated std::chrono around the per-clip candidate-scoring region
  (SPR loop + reroot loops), summed over ~50k clips, file-scope static (NOT
  thread_local; single-thread). Same-seed => identical clip/candidate counts =>
  the loop-ns delta isolates the cumulative per-candidate saving, excluding
  ratchet-reweight / R-glue / alloc variance that dominates wall. A compounded
  ~1-3% is resolvable there even if each lever alone is not.

### L-trivial-hoists — IMPLEMENTED, gating
Hoist per-call-invariant checks out of hot loops: `use_collapsed = !collapsed.empty()`
(replaces per-candidate `!collapsed.empty()` at 4 sites) + `revert_check`/`iw_scanchk`
getenv (per-clip / per-accept -> once per call). Byte-identical (invariants). Gate
running (MOD-only; BASE==documented already confirmed).

## Discovery workflow inventory (51 levers, 40 byte-identity-confirmed)
Adversarial verifiers consistently DEFLATE magnitudes -> the per-candidate path is a
long tail of sub-0.3% exact micro-levers (near/below noise) + a few uncertain
structural levers. 2 lenses FAILED (idle timeout), re-running as agents:
per-pass-precompute + structural-incremental-length (the big lever-b question).

**BANK (9, all trivial <0.1%):** structured-binding copy, scores-init dead-store,
collapsed-empty hoist, getenv hoists, do-reroot-nz/ns specialize, saved_postorder hoist.
**MEASURE-FIRST (13), meatiest:**
- `batch-scalar-x4-for-ratchet` — x4-batch the RATCHET scalar scorer (60% wall path);
  needs a NEW weighted-x4 kernel; verifier: latency-hiding only (reduce at-limit),
  may ~wash at 37-78t (L1-resident). MED risk. HIGHEST mission headroom.
- `skip-vroot-directional-copy-index-directly` — drop the per-clip vroot memcpy, index
  edge_set_buf directly. Real, small-moderate. MED (index alignment).
- `hoist-skip-checks-scalar` — hoist sub_edge-invariant skip predicates out of the
  scalar ei-loop. ~1-3% of scalar-loop self-time.
- int-accumulator / int-cutoff (EW double->int) — near-noise, MED risk.
**DESIGN-ONLY (5):** `ew-flat-x4-in-spr-loop` (batch EW SPR via existing T-245 kernel;
largest EW-SPR win, HIGH risk), `incremental-postorder-maintenance` (build_postorder
5.2% CPU, rebuilt per clip; HIGH risk).
**REJECT (24).**

PLAN: bank trivial bundle (on gate); then attempt ONE mission-impact structural lever
(batch-scalar-x4-for-ratchet OR ew-flat-x4-in-spr-loop) with the in-process
candidate-region chrono instrument for low-noise measurement; gate byte-identical.
Await structural-incremental-length verdict before committing to a major rewrite.

## Structural-incremental-length (lever-b) — DEAD (opus feasibility agent, decisive)
The standing "only substantial route" for the 2.5x per-candidate gap is NOT a
buildable prize:
- **O(1)-sliding (TNT quick-TBR): INFEASIBLE** — same non-invertibility / no-locality
  wall as L3b (footprint 41-68%, Euler delta 1.14-1.24x); AND cannot touch the
  irreducible up-pass (`up[D]=combine(up[parent],prelim[sib])`, O(N)/pass, which TNT
  also pays; ts_rate flat in N = already O(1)/cand amortized).
- **Skip-combine / 3-way-direct = sub-lever (d): FEASIBLE-BUT-WASHES**, worse than its
  own estimate. **LOAD-BEARING FACT:** on internal-node clips (dominate candidate
  work) each `edge_set_buf[below]` is REUSED `n_sub_edges`x by the reroot loop
  (vroot_cache, ts_tbr.cpp:1611 + 1651-1759) => materialization is amortized+beneficial;
  fusing/skipping the combine RECOMPUTES it per re-read = NET LOSS on the dominant path.
  Only tip-clip SPR (small fraction) benefits; combine is ~10% of EW, L2-bandwidth-bound
  (T-P5f), cache-resident at 37-78t (T-P5i fused=0). Realizable: low-single-digit% at
  best, likely washes, regresses reroot unless tip-clip-gated.
=> **No major TBR rewrite.** The "56% precompute/2.5x" framing conflated the
irreducible up-pass with the avoidable (~10%) combine. Reopen only at n_blocks>>4 /
molecular scale (same condition as L3b/M46). Converges with T-P5j + M46: per-candidate
path AT-LIMIT at this dataset class.

## Batching levers (flat-x4 -> ratchet/SPR) — REGRESS at mission scale (resolved by reading the kernel, no build)
`batch-scalar-x4-for-ratchet` and `ew-flat-x4-in-spr-loop` were the last levers
touching the 60% ratchet path / EW SPR with claimed headroom. Reading
`fitch_indirect_cached_flat_x4` (ts_fitch.cpp:709-739) resolves it: the x4 batch
issues **4 SEPARATE any_hit_reduce calls per block** (723-726) — it does NOT reduce
the reduce COUNT (the dominant, at-limit cost T-P5l), only interleaves them for
ILP/latency-hiding, and uses a COMBINED bail (734) that runs until ALL 4 exceed
cutoff. So vs scalar it: (a) saves only tiny block-loop/call overhead (n_blocks=4),
(b) LOSES per-candidate early-bail (~2.85/4 -> ~4/4 blocks = ~40% more at-limit
reduce work), (c) the latency-hiding benefit (its OWN stated rationale) WASHES when
vroot_cache is cache-resident at 37-78t (T-P5i). Net: REGRESSES at mission scale.
The existing flat-x4 reroot batch is a LARGE-TREE (180t+, vroot_cache > L1)
optimization; extending it to the ratchet/SPR path nets negative here. NOT built.

## vroot memcpy elimination — likely REGRESSES (per-pass agent)
Exact (byte-identical, 3 legs verified) but trades a one-time per-clip pack for an
extra `main_edges[ei].second` load on EVERY candidate access; each row reused
n_sub_edges x => per-access overhead > saved copy. The memcpy EXISTS to pack rows
for cache-friendly scan. NOT built (measure-first, but analysis says regress).

## !!!!! MAJOR CORRECTION (quiet-machine measurement overturns "sub-noise") !!!!!
The banked bundle (cutoff hoist + invariant hoists) is NOT sub-noise. On a QUIET
machine, base(da0f203f) vs mod(cutoff+hoists), 20 same-seed pairs/dataset, REPS=6:
- **Zanol2014: -13.2% median wall, 20/20 pairs faster, sign-test p=0.000**
- **Zhu2013:   -19.2% median wall, 20/20 pairs faster, p=0.000**
Byte-identical work (gate-proven) => pure per-op cost. The first A/B's "-0.15%
sub-noise" was an ARTIFACT of running under the 30-agent discovery workflow (±16%
load noise swamped everything). LESSON: do NOT trust verifier "sub-nanosecond"
deflations or loaded-machine A/Bs; quiet-machine measurement is authoritative.
**ROOT CAUSE — CONFIRMED by 3-way attribution A/B (base / nogetenv / mod, Zhu, 12 same-seed
triples, scores byte-identical): the per-clip `std::getenv("TS_REVERT_CHECK")`** (a DIAGNOSTIC
left in the per-clip teardown, ts_tbr.cpp:1852, ~100k+ calls/search). On Windows/ucrt getenv
is µs-scale (locked linear environment scan), NOT sub-ns.
- **getenv hoist ALONE: -19.1% (Zhu, 12/12 faster).** (Zanol cumulative -13.2%, 20/20.)
- **cutoff + collapsed hoists: +0.00%** (nogetenv == base) -> genuinely negligible (the verifiers
  were RIGHT about those; harmless cleanups, keep but they carry ~no win).
=> The ENTIRE 13-19% is the getenv hoist. REAL, banked, byte-identical.
**CAVEAT (honesty): getenv cost scales with ENVIRONMENT SIZE (it scans the env block) and is
platform-dependent (Windows/ucrt has a lock + UTF conversion; Linux cheaper). The 13-19% is
THIS test env (Rscript via Git Bash); a smaller env or Linux/Hamilton => smaller but still
strictly-positive, byte-identical win. Cross-platform magnitude => Hamilton confirmation
(queued for supervisor). Unambiguously beneficial regardless of magnitude.**
**WHY IT WAS MISSED: in the T-P5a VTune the getenv cost likely hid inside ucrtbase self-time
(mislabeled "memory traffic"), so no prior round flagged a per-clip getenv. Lesson: re-examine
ucrtbase/CRT self-time attribution for other hot-path stdlib calls.**
**PATTERN: diagnostic/env checks in per-clip/per-candidate hot loops are a hidden cost
class.** Hunted src/: cancel-file getenvs (ts_driven:654/ts_parallel:300) ALREADY
hoisted (fine); TS_EV_AUDIT (ts_tbr:895) convergence-frequency (fine); **ts_sector.cpp
TS_FREE_HTU_PROBE (802/848/895) + TS_SECT_DEBUG (1147) are PER-SECTOR — moderate, FLAG
for the sectorial agent (smaller than per-clip but free to hoist).**

## ===== (SUPERSEDED) earlier conclusion: per-candidate SCORING path at-limit =====
NB the SCORING-kernel verdict below still stands (lever-b dead, batching regresses,
reduce at-limit). What was WRONG was calling the banked control-flow hoists "sub-noise":
they are 13-19% via the getenv. The 2.5x gap framing is also revisited — a chunk of it
was this Windows getenv overhead, NOT a fundamental TNT per-candidate advantage.
## ===== FINAL CONCLUSION: TBR per-candidate/per-clip path is AT-LIMIT at mission scale =====
Exhaustive sweep (51 levers enumerated + adversarially verified, 2 deep agents,
direct kernel reads). Every route closed:
- **Per-candidate scorer reduce:** AT-LIMIT (T-P5l, AVX2 optimal n_states=9).
- **Structural incremental-length (lever-b):** DEAD — O(1)-slide infeasible (L3b wall +
  irreducible up-pass TNT also pays); skip-combine = sub-lever(d) regresses the
  dominant reroot path (view reused n_sub_edges x).
- **Lazy/incremental precompute (L3b, M46):** DEAD (no locality / union saturates).
- **Batching (flat-x4 -> ratchet/SPR):** regresses at mission scale (lost early-bail,
  cache-resident).
- **vroot memcpy elim:** likely regresses.
- **build_postorder incremental (5.2% CPU):** DESIGN-ONLY, order-dependent, HIGH risk —
  the ONE remaining item with real magnitude, NOT attempted unsupervised (correctness-
  critical postorder splice; needs supervised build + oracle). Flagged for review.
- **BANKED (exact, byte-identical, gate-proven):** the per-clip getenv hoist (in 3a50537e) =
  **13-19% of EW wall** (the headline; see MAJOR CORRECTION above). cutoff (6295c401) + collapsed
  hoists = ~0 (harmless). kept_ei (8291bbec) = marginal.
**=> The per-candidate SCORING gap vs TNT does NOT live in recoverable TBR work** (constant-factor
+ the NON-TBR half: ratchet reweight / sectorial / R-glue = Phase-0 + sectorial #39). **BUT a
13-19% per-clip getenv WALL overhead WAS recoverable and is now removed** — so part of the measured
2.5x TS-vs-TNT wall gap was Windows getenv overhead, not a per-candidate disadvantage. Reopen the
per-candidate STRUCTURAL levers only at much larger N / molecular scale (n_blocks>>4).

## ===== FINAL SUPERVISOR SUMMARY (branch claude/tbr-microlevers, off cpp-search da0f203f; NOT pushed/merged) =====
All commits gate-proven BYTE-IDENTICAL (score + candidates_evaluated identical, Wortley/Zhu/Zanol x seed{1,2}).

**THE WIN — MERGE THIS:**
- **3a50537e** perf(tbr): hoist call-invariant collapsed.empty()/getenv checks. Contains the
  per-clip `getenv("TS_REVERT_CHECK")` hoist = **13-19% of EW wall** (Zanol -13.2% 20/20 p=0,
  Zhu -19.1% 12/12 p=0; 3-way-attributed as the ENTIRE win). A diagnostic getenv left in the
  per-clip teardown, µs-scale on Windows/ucrt. Magnitude env-size/platform dependent (Hamilton
  confirmation owed) but byte-identical + strictly removes ~100k getenv/search => unambiguous.

**OPTIONAL / ~0 (exact, harmless, keep-or-drop your call):**
- **6295c401** perf(tbr): cutoff hoist. +0.00% (attribution proven). Exact cleanup.
- **8291bbec** perf(tbr): kept_ei (hoist sub_edge-invariant skip predicates out of reroot loops).
  MARGINAL: Zanol -0.1% (wash), Zhu -2.3% median (p=0.18). Byte-identical; scales favorably with
  tree size; reopen at larger N. Optional merge.
- **18d70dde / a5e434cf / 13c57946** chore/docs: sweep notes + gate/wall drivers + attribution.

**FLAGGED, NOT DONE (need supervised build):**
- `build_postorder_prealloc` incremental maintenance — 5.2% CPU, per-clip O(n) DFS rebuild. The
  ONE remaining real-magnitude TBR item. DESIGN-ONLY: order-dependent (postorder splice on
  clip/unclip), HIGH risk; a few-seed gate is insufficient. Supervised build + oracle if pursued.
- **ts_sector.cpp per-sector getenvs** (TS_FREE_HTU_PROBE 802/848/895, TS_SECT_DEBUG 1147) —
  SAME pattern as the TBR win, moderate frequency. -> sectorial agent (possible 30%-phase win).
- **Re-profile** (VTune): the 13-19% getenv removal invalidates the prior hotspot attribution
  (the getenv hid in ucrtbase self-time). A fresh quiet symboled VTune should reveal the new
  distribution + any other mislabeled hot-path stdlib.
- **Hamilton/Linux** cross-platform magnitude confirmation for the getenv win.

**NET:** the SCORING-kernel verdict stands (per-candidate path at-limit: lever-b dead, batching
regresses, reduce at-limit T-P5l) — BUT a 13-19% per-clip getenv WALL overhead was found+removed,
which the prior "at-limit/sub-noise" framing missed. Part of the measured 2.5x TS-vs-TNT wall gap
was this Windows getenv overhead, not a fundamental per-candidate disadvantage.
