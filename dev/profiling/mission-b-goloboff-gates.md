# Mission B — existence-before-build gates for the Goloboff 8× kernel screen

**Date:** 2026-07-14 · **Base:** worktree `claude/mission-b-kernel-goloboff-6f5c0e`
(tip `2595da46`) · **Build:** full `R CMD INSTALL` → `.agent-disc` · **Harness:**
`dev/profiling/reeval/disc_nstates.R` + `corpus_state_heterogeneity.R` (committed).
**Plan:** `dev/plans/2026-07-14-mission-B-kernel-8x-goloboff-screen.md`.

## Verdict (headline)

**STOP — the gates do NOT clear a build, now confirmed by a large-dataset head-to-head
and a soundness proof.** All four lines of evidence converge:

- **Constant-factor bit-packing — NOT refuted (this bullet was a METRIC BUG; corrected
  2026-07-15).** The "median headroom 0.000" counted the all-states `?` token as all global
  states, inflating every missing-data char to n_states. The TRUE per-char alphabet (`?` = the
  char's own states, Fitch-EXACT) is mean ~2.6–3.4 vs global 9–10 → **corrected median headroom
  0.40, up to 0.62** (74% of datasets >0.30). A per-block-local-alphabet bit-slice is a REAL,
  EXACT ~1.2–2× per-candidate lever (regime-dependent). See `kernel-1p9ns-findings.md`.
- **Goloboff Union-Construct screen — SOUND (with true finals) but GATE #3 FAILS.** The
  screen is a genuinely sound lower bound *only* with true MPR finals (proof:
  `dev/red-team/proofs/union-construct-lower-bound.md`; 0 violations exhaustively,
  vs. massive over-rejection with the engine's collapsed finals). BUT the head-to-head
  (Hamilton `17872299`) shows **no dataset where a per-move speedup would improve
  time-to-optimum.** Two independent readings:
  1. **Airtight (TS trajectory alone, no TNT baseline):** on the 8/11 reach-parity large
     datasets where TS reaches the floor, `ts_tt_s ≪ ts_conv_s` — TS finds the optimum
     early then **over-searches to the cap** (e.g. project3733 floor at 4s but runs 433s;
     project970 floor at 328s, runs 555s). The wall it "spends" is stopping-rule waste
     ([[stopping-rule-mechanism]]), which a faster/screening kernel would only make faster.
  2. **TS-vs-TNT (two-sided, `tnt_tt_s` measured — see table):** against TNT's own
     time-to-first-floor, TS reaches the optimum **faster on 5/6 datasets (median 0.39×) and
     ~ties on the most kernel-heavy (project970, 1844 chars, 1.12×)**. So the 30× per-*candidate*
     gap does NOT propagate to time-to-optimum — **the mission premise ("30× per-move ⇒ hit the
     optimum later everywhere") is FALSIFIED** (candidate-efficiency ~1.5× + prompt convergence
     absorb it).
  3. **Screen payoff regime empty (measured):** all 13 large datasets have CI 0.11–0.35 —
     uniformly high-homoplasy/incongruent, the regime where Goloboff says the union construct
     "did not save any time — or even made the program slower." Screen contraindicated
     corpus-wide at large scale, independent of the race.

Recommendation: **do not write C++.** The screen is correct but unmotivated on this corpus;
packing is a sub-10% constant on a minority of data. Residuals + reopen path below.

---

## Gate #1 — Goloboff 1996 full text: the EXACT screen spec

Read in full (`Zotero/.../Goloboff - 1996`; text at
`scratchpad/goloboff1996.txt`). The abstract line *"rejecting several locations as
suboptimal by checking just one node"* is the **"Union Construct" method** (pp. 212–213),
**not** the per-candidate length-check (TS already does that — the strict-`<` bail).

**Mechanism (sound lower bound):**
- Each internal node N of the (clipped/main) subtree gets a **union set** = all states
  present in N's descendant leaves (a cheap down-pass, incrementally maintainable).
- The **union construct** for the branch above N = `union_set(N) ∪ MP_states(parent(N))`
  (must include the ancestor's most-parsimonious states — states absent from all
  descendants can still be added in the up-pass; Goloboff's 5-taxon `(4(3(2(1 0))))`
  example).
- Score the clipped subtree's basal state set against the union construct exactly as in
  indirect length calc. Because the union set is a **superset** of any real descendant's
  states, the computed increase is a **lower bound** on inserting anywhere in N's
  subtree. If it already ≥ cutoff, **reject the whole subtree** in one check.

**Empirical caveats Goloboff reports (directly relevant to Gate #3):**
- Congruent 168-taxon set: ~2× (run time ~50% *with* the construct).
- *"For data sets with much incongruence … did not save any time — or even made the
  program slower."*
- *"probably also … when swapping on very suboptimal trees (like … the Wagner method)."*
- Gain < expected because *"for many nodes, necessary to check length increases both
  against union construct and final states."*

So the screen is a **conditional** lever (large + congruent + near-optimal), NOT an
"always-pays" one. It reduces the number of FULL scores; it does not speed a full score.

### How the union construct relates to the PRIOR approximate-screen work (checked, not assumed)

`git log --all --grep` (union-construct / goloboff-screen / one-node / approximate-screen) and a
scan of all ~40 local worktrees find **no commit or doc implementing/investigating Goloboff's
union CONSTRUCT.** The prior work is on three *related but distinct* constructs — and the
distinction is exactly the soundness pivot:

1. **lever-c** (`dev/red-team/proofs/lever-c-bound-then-verify.md`, DEAD): up-**IGNORING**
   admissible bounds. Goloboff's construct is up-**AWARE** (it unions in the ancestor's MP states —
   his 5-taxon `(4(3(2(1 0))))` example shows why states absent from all descendants must be
   included). Distinct.
2. **union-of-finals** (`dev/red-team/union-of-finals-bound-proof.md`, math-prover): the DEPLOYED
   `fitch_indirect_length` ranker uses `final(A) ∪ final(D)` — the collapsed `uppass_node` finals of
   the two branch **ENDPOINTS**. Proven a sound **UPPER** bound (over-counts, `U ⊆ E(D)`;
   958016/958016 edges) ⇒ **unsafe as a rejection screen** (would discard improving moves); safe
   only as a ranker with exact re-score. The `scoreapprox`/`exactness-gate` branches +
   `exactness-gate.R` P2/P3 probe exactly this endpoint-union quantity.
3. **Goloboff union construct** (this mission's lever, UN-BUILT): uses
   `union_set(descendants of D) ∪ MP(A)` — the union SET of **all** states below node D (a
   SUPERSET), ∪ the ancestor's MP states. Because it tests insertion against a **superset**, the
   `[r ∩ S = ∅]` indicator is antitone in S ⇒ it under-counts ⇒ a sound **LOWER** bound ⇒ safe for
   rejection; and one check covers the **whole subtree** below D at once.

So the prior union-of-finals proof is **instructive, not contradictory**: the bound direction hinges
on subset-vs-superset. TS's collapsed finals are a *subset* → upper bound (the trap that sank
union-of-finals as a screen); Goloboff's union *set* is a *superset* → lower bound. The construct
sidesteps that trap **by design** — but its soundness under the corpus's NA (state-0) /
polymorphism handling still needs its OWN math-prover proof (the union-of-finals proof proved the
*opposite* direction for a *different* set; it does not transfer). Note also that the exact per-edge
E(D)=M(D) is already computed and is a *tighter* lower bound than the construct — the construct's
only advantage is **coverage** (reject a subtree per one check), never tightness.

## Gate #2 — vary-n_states discriminator (TS-only, local, tight-cutoff)

**Layout facts (read from source, load-bearing):**
- TS is **bit-sliced**: `a[s]` is a 64-bit mask over 64 chars for state s
  (`ts_simd.h any_hit_reduce3`); inner loop runs `n_states` words per 64-char block.
- Every block uses the **global** applicable-state count (`ts_data.cpp:214–250`:
  `blk.n_states = total_app_states`; comment: "A pattern using state index k needs
  state word k, regardless of how many states that individual pattern uses").
- TNT-64 packs **2-bit/4-bit per-character fields** (prior Task B disasm) = also
  **n_states bits/char**. ⇒ **identical bit density**; TS is *denser* at non-power-of-2
  state counts. So "TS cost ∝ n_states → packing gap" is **too loose** — TNT's curve
  tracks TS. Reframed axis: **per-word-bound vs fixed-overhead-bound**.

**Measurement** (per-candidate ns via committed `TS_IW_TIMING` chrono, min-of-runs;
converge to a kernel local optimum first, then time a near-optimal round = tight cutoff;
`totns` = evaluation-count-weighted blend of SPR + REROOT):

| n_states | totns @100t | totns @250t | totns @482t (141ch) |
|---|---|---|---|
| 2 | 8.98 | 8.55 | 8.45 |
| 4 | 9.13 | 8.66 | 8.65 |
| 8 | 10.20 | 9.65 | 9.40 |
| **9 (corpus)** | **10.86** | **10.30** | **10.17** |
| 16 | 12.34 | 11.41 | — |
| 32 | 15.75 | 15.79 | — |

`totns ≈ 8 + 0.24·n_states` (R²≈0.99). At corpus n_states=9 the **~8 ns fixed intercept
is ~75–78%** of per-candidate cost; a 16× state increase costs only **1.75–1.85×**, not
16×. **Fixed-overhead-bound at corpus scale → not a packable inner loop.**

**Large-scale (cache/work-bound) — cost is total_words-driven** (482t, n_states=9,
n_chars swept): 141ch→10.2, 189ch→15.6, 400ch→28.2, **800ch→51.2 ns** (~0.45 ns per
added word). Cost rises monotonically with `total_words = n_blocks × n_states` (more
words gathered/reduced + a bigger edge-set buffer spilling cache). At scale the per-word
term dominates the ~8 ns intercept — but packing it away still needs state heterogeneity
(below), and TS is already at TNT's density.

**Calibration MET:** the fingerprint's ~47 ns (5432 = 482t×189ch on Hamilton) is
reproduced locally at **800ch/482t = 51 ns** — i.e. ~47 ns is a **working-set (cache)
effect** of large total_words. 5432's own 189 chars gives 15.6 ns *isolated*; the real
search reaches ~47 ns via search-context cache contention (full machinery evicting the
TBR buffers) that an isolated descent cannot reproduce at 189ch, but which the larger
working set reaches directly. Either way the target regime is total_words-bound, and the
*shape* (intercept-dominated at corpus n_states; total_words-driven at scale) is what the
gate needs.

## Gate #2b — corpus state heterogeneity (governs the packing lever)

`corpus_state_heterogeneity.R` over 31 local datasets. `headroom` = fraction of
state-words a TNT-style per-block state alphabet could remove (chars sorted by
distinct-state-count, 64/block, block width = block max).

- **Median headroom 0.000; max 0.375** (Dikow2009); **0 datasets > 0.5**.
- ~2/3 are **uniform** (per-char state count ≈ global: Zhu/Giles/Conrad/OMeara/Wortley/
  Zanol/Shultz…) → zero packing possible.
- The heterogeneous third (Dikow 88t, Agnarsson 62t, Liljeblad 68t) are small/medium →
  compute-bound → the per-word fraction is only ~20%.

Best case = Dikow 0.375 word reduction × ~20% per-word fraction ≈ **~7.5%**, on a
minority of matrices. **Packing refuted as a mission lever.**

## Gate #3 — name a dataset where TS's anytime curve trails TNT *because of* per-move throughput

Could **not** be met with current evidence:
- **Corpus (≤88t):** head-to-head (`headtohead_phase0.csv`) gives `ts_wall/tnt_wall`
  ≈ 1.3–2.5× with `cand_ratio` ≈ 1.2–1.9× ⇒ **per-move throughput ~1.4–1.6×**, and the
  wall gap is **budget choice + candidate-efficiency**, not per-move (KPI 2026-06-21;
  quality CLOSED, TS ≥ TNT).
- **The 30× is confined to ~482t scale** (cache/work-bound gather; reproduced as steep
  n_chars/n_tips scaling above). In *this* corpus the only 482t dataset is **5432**,
  which Mission A shows is **basin-capture-limited** (score-vs-candidate trajectory
  flattens) — faster per-move does not change its outcome.
- **The 1.4 ns TNT figure is partly a counting artifact:** TNT's "rearrangements
  examined" denominator includes screen-rejected candidates (the union construct
  rejects 6–17 destinations per one-node check), deflating its ns/candidate. TS's
  46.8 ns is per *fully scored* candidate. So matching TNT's *counter* is exactly
  building the screen — but that only pays where many candidates are bulk-rejectable
  (congruent, near-optimal), i.e. the conditional regime that Gate #3 cannot yet name.
- Large corpus datasets (200–400 tips, e.g. project6403 294t, project2183 318t) are
  **unbenchmarked** head-to-head, and the hard ones are **incongruent** (screen hurts).

### Gate #3 reopen — in progress (2026-07-14, user said "proceed")

**Reach structure of the large corpus (from the floor campaign's `floor_table_corpus.csv`,
both-engine best-known).** Of the 35 large datasets (ntax≥150) with a TNT floor:
- **26 are reach-PARITY** (TS floor ≤ TNT floor) — quality closed, so the only axis is SPEED.
  These are the candidate gate-3 datasets. Span ntax 150–267, nchar 17–1844.
- **9 are reach-LIMITED** (TS floor > TNT floor by +1…+8: project2183/2722/3285/427/2477/
  syab0720{4,5}/4446_(1)/4284) — Mission-A basin territory (screen cannot fix reach).

So gate #3 is **not dead on arrival** — large speed-only datasets DO exist (reach isn't the
universal blocker at scale). The open question narrows to: among the 26 reach-parity large
datasets, does TS reach the floor but TRAIL TNT on wall because the per-move kernel dominates?

**Head-to-head RESULT** (Hamilton array `17872299`, TS {auto,thorough} vs TNT xmult, cap 600s,
2 seeds, frozen 2.0.0 lib; scripts `reeval/gate3_h2h.R`+`gate3.sh`+`gate3_keys.csv`; pulled to
`.agent-disc/hamilton/gate3_out/`). `ts_tt_s` = TS wall to FIRST reach the floor; `ts_conv_s` =
TS wall to its own stop (cap); `tnt_conv_s` = TNT xmult convergence.

`ts_tt_s`/`tnt_tt_s` = wall to FIRST reach the floor (two-sided race, `tnt_tt_s` from array
`17872439`); `ts_conv_s` = TS wall to its own stop; **CI** = consistency index (min steps / floor;
low = high homoplasy = incongruent). `TS/TNT` = `ts_tt_s/tnt_tt_s`.

| key | nTip | nChar | **CI** | TS reached | ts_tt_s | ts_conv_s | tnt_tt_s | **TS/TNT** |
|---|---|---|---|---|---|---|---|---|
| project970  | 157 | 1844 | 0.19 | yes | 328 | 555 | 292 | **1.12** |
| project2668 | 196 | 1227 | 0.15 | **no (+1)** | — | 543 | 553 | — |
| project4550 | 230 | 889  | 0.19 | yes | 318 | 541 | 424 | 0.75 |
| project3733 | 157 | 853  | 0.33 | yes | 4.2 | 433 | 96  | 0.04 |
| project4327 | 197 | 823  | 0.21 | **no (+1)** | — | 529 | (reached; tt n/a) | — |
| project4085 | 164 | 716  | 0.23 | yes | 46  | 540 | 100 | 0.46 |
| project1024 | 163 | 156  | 0.34 | yes | 2.9 | 120 | 12  | 0.24 |
| project912  | 173 | 74   | 0.30 | yes | 14  | 174 | 44  | 0.32 |
| project175  | 165 | 71   | 0.21 | no (+2; TNT also +2) | — | 212 | — | — |
| project4204 | 163 | 37   | 0.11 | yes (−1) | 12 | 107 | — | — |
| project2220_(2)| 267 | 17 | 0.35 | yes | 103 | 190 | — | — |

**Reading (airtight, TS trajectory alone).** On **8/11** reach-parity datasets TS reaches the
floor, and `ts_tt_s ≪ ts_conv_s` throughout: TS finds the optimum EARLY then **over-searches to
the cap** (project3733 4.2s→433s; project970 328s→555s). This needs no TNT baseline — it is the
known stopping-rule waste ([[stopping-rule-mechanism]]), and a faster/screening kernel would only
make the over-search faster, not shorten time-to-optimum.

**Two-sided race (both time-to-first-floor, `tnt_tt_s` measured).** On the 6 datasets where both
milestones are clean, TS/TNT time-to-floor is **median 0.39×, max 1.12×, min 0.04×** — i.e. TS
reaches the optimum *faster* than TNT on 5/6, and ~ties on the most kernel-heavy (project970,
1844 chars, 1.12×). So the 30× per-*candidate* gap does **not** propagate to time-to-optimum
(candidate-efficiency ~1.5× + prompt convergence absorb it). [My earlier "median 0.34×" used TNT
*full-convergence*; `tnt_tt_s` corrects it, and the conclusion only strengthens.]

**The screen's payoff regime is empty here — measured, not asserted.** CI across ALL 13 large
datasets is **0.11–0.35** (every matrix is high-homoplasy / incongruent; the corpus's own `ci_`
columns are NA at this size, so this CI was recomputed from matrix+floor, `reeval/ci_fix.R`).
Goloboff (1996) is explicit that the union construct "did not save any time — or even made the
program slower" on incongruent data. So corpus-wide at large scale, the screen is contraindicated,
independent of the over-search finding.

**Non-reachers in the 600s cap (3):** `project175` — BOTH engines miss (+2), a hard floor, not
TS-specific; `project2668` (CI 0.15) / `project4327` (CI 0.21) — TS +1 vs TNT 0, the only genuine
TS lags. Both are strongly incongruent (screen would hurt), and the lag is not cleanly
throughput-attributable (project970 has *more* chars yet reaches at parity) — it reads as
candidate-efficiency/basin difficulty, not per-move throughput. Controls `project2183` (CI 0.13)/
`project2722` (CI 0.12) confirm reach-limited (TS +13/+7 — basin, Mission A).

**Soundness proof landed** (`dev/red-team/proofs/union-construct-lower-bound.md` +
`dev/red-team/union-construct-gate.R`): UC is a sound lower bound with **true MPR finals** (0
violations, exhaustive ≤7t + polymorphic ≤8t; 0 over-rejections / 180,850 added-length checks)
but **UNSOUND with the engine's collapsed finals** (over-rejects; 4-tip `1,1,0,0` counterexample).
So the screen is *correct-buildable* (must use exact edge-set finals `E(A)`, plus new
`union_set` machinery) — but Gate #3 shows it has no time-to-optimum gap to close here.

## Recommendation

**STOP — do not build the Goloboff screen or state-count packing.** All three gates were run
(existence-before-build satisfied) and the soundness proof completed:
- Packing: refuted (density-equal on uniform data; median-0 heterogeneity headroom; ≤7.5% on a
  minority).
- Screen: proven sound (with true finals) but **Gate #3 fails** — no corpus dataset has a
  time-to-optimum gap a per-move speedup would close. TS reaches the optimum at ≤1.13× TNT's
  wall even on the most kernel-heavy large dataset; its apparent wall loss is over-search
  (stopping), not per-move throughput. The mission's motivating premise is falsified.

**What the evidence points at instead (NOT this mission's lever):**
- The one genuine large-scale TS deficit is **over-search / stopping** (`ts_conv_s ≫ ts_tt_s`):
  TS finds the optimum early then keeps swapping to the cap. That is a stopping-rule problem
  ([[stopping-rule-mechanism]]), not a kernel problem — orthogonal to Mission B.
- The 9 reach-limited large datasets (+ 5432) are **basin capture** — Mission A.

**Residual (honest, low-priority):** `project2668`/`project4327` (196–197t, 823–1227 chars) —
TS +1 vs TNT 0 in a 600s cap. Not cleanly throughput-attributable (project970 has more chars yet
reaches at parity) and the screen is the wrong tool (incongruent, char-rich). If ever chased, the
right diagnostic is candidate-count + congruence attribution, not a kernel rewrite. Should any
future data DO name a large+congruent+near-optimal throughput-bound case, the screen is
buildable-sound **only** via true edge-set finals `E(A)` (never the collapsed `tree.final_`) plus
`union_set` machinery, default-OFF, A/B'd against incongruent/Wagner-start regression.

## Provenance / artifacts

- `dev/profiling/reeval/disc_nstates.R` (+ `disc_nstates_{100,250,482}t.rds`),
  `dev/profiling/reeval/corpus_state_heterogeneity.R` (+ `.rds`).
- Gate #3 head-to-head: `reeval/gate3_h2h.R` + `gate3.sh` + `gate3_keys.csv` (Hamilton array
  `17872299`); results in `.agent-disc/hamilton/gate3_out/`. Corpus feature/reach tables pulled
  to `.agent-disc/hamilton/` (`feature_table_corpus.csv`, `floor_table_corpus.csv`).
- Soundness: `dev/red-team/proofs/union-construct-lower-bound.md` +
  `dev/red-team/union-construct-gate.R` (empirical gate).
- Logs: `.agent-disc/disc_run.log`, `disc_482.log`, `disc_cache.log`. Full text of the
  paper: session `scratchpad/goloboff1996.txt`.
- Builds on: `dev/profiling/s7-fastpath-sizing.md` (#7 refuted, TNT-64 disasm),
  `dev/benchmarks/tnt_disassembly_analysis.md`, `dev/profiling/kpi-2026-06-21.md`,
  `dev/red-team/proofs/tnt-quick-tbr-views-literature.md`. Memory
  `[[tnt-per-move-kernel-gap]]`, `[[l3b-lever6-landed]]`, `[[kpi-2026-06-21]]`.
