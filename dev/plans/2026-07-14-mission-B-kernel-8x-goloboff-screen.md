# Mission B — 8× faster per-move kernel (objective (ii): be as fast as possible)

**Framing (2026-07-14).** The SPEED mission, split off the per-move investigation. A faster per-candidate
TBR scoring kernel serves the WHOLE corpus (hit the optimum sooner everywhere), independent of 5432-reach
(that is Mission A). "As fast as possible" is a standing objective, so kernel speed always pays — subject
to the gating discipline below.

## The gap (measured, not inferred)
TS's per-candidate TBR reconnection scorer is **~30–33× slower than TNT's** (corrected from an earlier
confounded 43×): near-optimal regime TS **46.8 ns/candidate** vs TNT **1.4 ns**; fingerprint job 17866795
confirmed ~28–42× roughly CONSTANT across n_chars 47/94/141/189. This is a raw KERNEL gap — a per-move
**screen and/or constant-factor packing** — NOT a missing operator (see [[tnt-per-move-kernel-gap]]).

## The live lever — the Goloboff 1996 UP-AWARE approximate screen
"Reject several locations as suboptimal by checking just ONE node" (Goloboff 1996, *Cladistics* 12:199;
abstract: approximate method **3–8× faster than exact**). Per the red-team note
`dev/red-team/proofs/tnt-quick-tbr-views-literature.md`, this is **"the one genuinely live thread."**
It is DISTINCT from:
- the DEAD lever-c (which killed only up-IGNORING admissible bounds), and
- TS's exact directional edge_set `E[D]=combine(prelim,up)` (up-aware but EXACT/expensive).
An up-aware *approximate* screen sits between them: cheaper than the exact edge_set, sound enough to reject
most non-improving candidates by touching one node, with exact re-score of survivors.

## Existence-before-build — MANDATORY gate (the campaign's hard rule; the advisor's explicit ask)
Do these THREE before writing any C++:
1. **Read Goloboff 1996 full text** (user has Durham–Wiley access). Get the EXACT screen spec — do NOT
   reverse-engineer from the abstract.
2. **Run the vary-n_STATES discriminator (NOT yet run).** The fingerprint could not separate "one-node
   screen" (buildable, ~flat in n_states) from "constant-factor packing" (TS cost ∝ n_states: 32-state →
   ~33 words/64-char block, vs TNT's 2/4-bit packed fields). Hold n_chars fixed, vary n_states: flat →
   it's a screen; ∝ n_states → the gap is packing. This says WHICH gap you're closing.
3. **Name the dataset(s) where TS is kernel-bound AND losing the anytime race.** The corpus wall gap is
   "default-budget mismatch, NOT algorithmic" ([[kpi-2026-06-21]]) and quality is CLOSED (TS≥TNT except
   5432). So a kernel win may only make already-won curves faster. State concretely where it PAYS
   (e.g. a large matrix where TS's anytime curve trails TNT's and the cause is per-move throughput), or
   the lever is speculative. Mission A's 5432 is NOT it — its curve is basin-capture-limited, not
   throughput-limited (the score-vs-candidate trajectory flattens).

## Exactness discipline (a screen is a soundness obligation, NOT a head start)
An approximate screen must be a **SOUND bound**: a lower bound on insertion cost, used to reject a
candidate ONLY when the bound already ≥ cutoff, with exact re-score of survivors. NOTE: the union-of-finals
work PROVED the up-IGNORING screen is UNSOUND (it over-counts — `dev/red-team/union-of-finals-bound-proof.md`).
That is a warning, not scaffolding: an up-aware approximate screen needs its OWN soundness proof. Use the
`dev/profiling/exactness-gate.R` harness as the empirical check and the math-prover lane for the proof
before any deployment. Byte-identity is NOT required (it's a screen), but soundness IS.

## Alternative / complementary lever — constant-factor bit-packing
If the vary-n_states discriminator shows the gap is dominated by PACKING (TS ∝ n_states): a state-count-aware
bit-packing (à la TNT's 2/4-bit fields for low-state data) is a bigger rewrite but a real corpus-wide
constant-factor win, especially on binary / few-state matrices. Gate on the discriminator; don't build both.

## Context — what is already LANDED (don't repeat)
- #7 scorer monomorphization (merged, ~5% whole-search, exact): removed dead per-candidate flavour
  dispatch. `dev/profiling/s7-land.md`.
- #6 incremental edge-set (~1.5× raw TBR, ~1× production wash — insurance, default-OFF): `dev/plans/
  2026-07-14-lever6-incremental-edgeset-land.md`, [[l3b-lever6-landed]].
- These addressed ORCHESTRATION + edge-set recompute (~5% + wash), NOT the 30× core. The 30× core is the
  SCREEN and/or PACKING above.

## Discipline
No local heavy compute (Hamilton). Existence-before-build (all three gates above). Sound-bound-or-exact
only. Verify anchors against the current tip.

## Durable artifacts (Hamilton `/nobackup/pjjg18/`)
`tnt_disasm_64.txt` (+`_intel`); `reeval/` micro_tbr_k, fp_gen.R, fp_run.sh (fingerprint harness — extend
it to vary n_states). Red-team note `dev/red-team/proofs/tnt-quick-tbr-views-literature.md`. Profiling:
`dev/profiling/s7-fastpath-sizing.md`, `s7-land.md`. Memory [[tnt-per-move-kernel-gap]],
[[sectorial-throughput-at-limit]], [[getenv-ucrt-cost]].
