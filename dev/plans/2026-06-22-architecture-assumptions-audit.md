# Architectural-assumption audit (2026-06-22)

**Question asked:** What structural decisions did we commit to early in development —
data representation, search shape — that the optimization effort has taken for
granted? The kind of thing an agent skips because it's "out of scope," "too much
work to change this late," or simply never named as a decision.

**Method:** Read the core C++ representation (`ts_tree.h`, `ts_data.h`, `ts_tbr.h`,
`ts_data.cpp`) and dispatched three mapping passes — R/C++ boundary, search-population
structure, binary/rooted commitments. Cross-checked against the profiling/red-team
memory so "discovery" is separated from "known-and-parked."

The findings are split along **two axes the original question conflates**:
- **Reachability** — what the search can *find at all* (bounds the QUALITY mission).
- **Throughput** — where time is *spent* (bounds speed within fixed reachability).

---

## A. Premise correction (the most useful single finding)

**The worked example — "we represent trees as edge matrices throughout all search
operations" — is empirically false.** The hot path never touches an edge matrix.

- Internal representation (`ts_tree.h:29-66`): flat `parent[]`, `left[]`, `right[]`
  index arrays plus node-major bit-packed state sets (`prelim`, `final_`, `down2`,
  `subtree_actives`). Rearrangements mutate these in place with a snapshot/undo stack
  (`clip_undo_stack`, `PreallocUndo`, `ts_tree.h:68-128`).
- Edge matrices exist **only at the R↔C++ I/O boundary**: `init_from_edge` on warm-start
  (replicate 0 only, `ts_driven.cpp:748-754`) and `tree_to_edge` on output
  (`ts_rcpp.cpp:116-131`). Nothing in TBR/SPR/NNI/sector/ratchet/fuse reads an edge matrix.

**Corollary that kills the user's hypothesis:** a pointer-based ("nodes-as-pointers")
tree would buy ≈nothing here. The dominant per-move cost is **state-set recomputation**
(`compute_insertion_edge_sets` ≈ 30% of EW CPU, memory T-P5o), not topology bookkeeping.
We already moved off the slow part of topology editing (per-clip heap allocation →
`PreallocUndo`). Topology mutation is a few index swaps; pointers would not make it cheaper,
and would *lose* the contiguous node-major memory layout the SIMD scorer depends on. So the
specific structural worry in the prompt is not where the leverage is — but the *instinct*
(question the representation, not the inner loop) points at real items below.

---

## B. Reachability decisions — bound what the search can find

These are the candidates most likely to bound the actual mission (best-possible trees,
not parity). None of them are throughput levers; they change the *set of reachable optima*.

### B1. Single-tree multistart + passive pool, not a maintained diverse population  ★ HEADLINE

**What we committed to.** The search is a serial sequence of *independent single-tree
hill-climbs*. Each replicate starts from a fresh Wagner tree, runs the whole pipeline
(Wagner → NNI/SPR → TBR → sectorial → ratchet → drift → fuse) on **one** in-place
`TreeState` (`ts_driven.cpp:54-620`), and dumps the result into a pool. The `TreePool`
(`ts_pool.*`) is a **write-once deduplicated archive** with capacity ~100 and
diversity-aware eviction — not a population that is actively toured.

The pool *does* feed back into the search, so this is not "purely independent restarts":
- intra-replicate fusing pulls pool donors (`ts_driven.cpp:568-590`),
- periodic inter-replicate fusing every `fuse_interval` (`ts_driven.cpp:954-1009`),
- consensus-tightening constraint from pool strict consensus (`ts_driven.cpp:866-876`),
- conflict-guided sector selection from pool split frequencies (`ts_sector.cpp:99-188`),
- Thompson-sampling over Wagner start strategies (`ts_driven.cpp:757-813`).

**What is structurally absent.** Sectorial search — our highest-value escape mechanism —
runs on the **single current tree** and reduces space *within that one tree*
(`ts_driven.cpp:246`, `ts_sector.cpp`). TNT's documented strength (memory:
[tnt-sectsch-superpower]) is sectorial search *over a retained diverse set* — it holds many
basins simultaneously and recombines across them mid-search. We approximate that only
through the indirect channels above; we never hold and tour a diverse elite set.

**Why it's been taken for granted.** Every profiling round has measured throughput
*within a single hill-climb* and concluded "AT-LIMIT." The shape of the search above that —
one tree at a time — was never on the table; it's the substrate the profiler runs inside.
The memory *diagnoses* TNT's diverse-set advantage but there's no evidence the search was
ever restructured around it.

**Status: genuinely under-examined as an architecture.** This is the headline because it
is both consequential (it's a *quality* mechanism, and quality is the mission) and the one
big decision no profiling round could see.

**But it is a claim that must be tied to a LIVE gap before acting.** KPI memory says quality
is "largely CLOSED" (TS ≥ TNT except the reliable-1261 config on Zanol). So the consequence
is an **open question, not a verdict**:

> *Does the residual gap appear in multimodal cases where independent restarts keep
> re-finding the same basins — exactly what a retained-and-toured diverse set would escape?*

**Discriminating check (cheap, do before any rebuild):** instrument basin diversity across
restarts — how many *distinct* topologies the pool holds at convergence, and how often a
fresh restart rediscovers an existing pool member rather than a new basin. Run it on the
configs where TNT still wins per-iteration (memory: any per-iteration TNT win is diagnostic,
[tnt-outperformance-is-diagnostic]). If restarts collapse onto a few basins there → B1 is
bounding quality and a population-based restructure is the highest-leverage move available.
If basins are already diverse → demote B1 to "structural but not currently bounding."

### B2. Binary-only resolution; no native collapsed-space search

**What we committed to.** `n_internal = n_tip-1`, `n_node = 2·n_tip-1`, `left[]/right[]`
(`ts_tree.h:31-39`) hardcode strictly bifurcating trees. Input polytomies are auto-resolved
(`MakeTreeBinary`, `MaximizeParsimony.R:743-752`); the search only ever produces binary trees.

**Consequence.** Parsimony is naturally defined over *collapsed* trees: zero-length
(zero-support) branches should not be distinguished. We instead search the space of binary
*resolutions*, so distinct resolutions of the same true polytomy are explored as separate
topologies — wasted candidate evaluations and a divergence from how TNT searches (collapsed
space). `ts_collapsed.*` marks zero-length edges to *skip clips* (`ts_tbr.cpp:61,2150`) but
does **not** contract them or merge regraft positions; the topology stays binary.

**Status: known as a feature, under-framed as a search decision.** There is a planned branch
(`dev/plans/2026-03-22-...-full-polytomy-search.md`, "Approach B" = keep binary internals,
prune harder). The reframe worth recording: this is not only an output-format nicety — it's a
*search-efficiency and reachability* decision (collapsed-space search changes the neighbourhood
graph), and it interacts with B1 (pool dedup already uses collapsed-form hashing,
`ts_driven.cpp:913-921`).

**CLOSED for the morphological roster 2026-06-22 (efficiency probe under the wall-clock
mission `[[mission-wallclock-to-optimum]]`).** Two findings invert the plan's framing:
1. **The plan's "main win" (Phase 3 regraft-region merging) is ALREADY BUILT AND ACTIVE** in
   production hill-climbing / sectorial / ratchet TBR (`ts_tbr.cpp:1581` SPR-loop skip +
   `:1692` reroot `kept_ei` filter; `use_collapsed` is on whenever `collect_pool==nullptr`,
   i.e. every optimisation pass — only the final MPT-enumeration pass at :1168 disables it).
   So B2-Approach-B is not a 16–24 day greenfield build; its core is running.
2. **It captures ~nothing on this data class, two independent ways.** The conservative,
   score-exact merge (the only *free-lunch* lever — skips provably-identical-score positions)
   fires on **0%** of work: `n_zero_skipped = 0` over 11.8M candidate evals (5 random-start
   descents), AND an env-gated Condition-1 counter (`TS_B2_CEILING`, `ts_tbr.cpp`) reads
   **0.000%** over 620k×6 regraft evals. `compute_collapsed_flags`' Condition 1 (zero *downpass*
   cost at parent) + Condition 3 (`prelim[sib]==prelim[parent]`) is precisely calibrated for
   identical-score safety and is NOT safely relaxable.
The remaining lever is the *aggressive* TNT `collapse 3` (collapse minimum-length-zero edges;
density 5.6% at the optimum / 27.8% at a random start, `b2_collapsed_density.R`). It is NOT a
free lunch: it skips *different-score* resolutions (Goloboff's asymmetric reachability), so on
**closed-quality data the speedup and the quality risk are the same mechanism** — the magnitude
gate doesn't even apply. Three independent reasons it's a wash here: (a) B1 located the residual
in *per-restart-reach* (restart count), not within-restart wasted evaluation, which a smaller
neighbourhood doesn't touch; (b) the density sits at 5.6% exactly at the optimum where
wall-clock concentrates; (c) the plan's own literature: union-construct gave 50% on *congruent*
168t data but **"no gain on incongruent data"** — morphological matrices are incongruent.
**Verdict: no free efficiency win; the only lever is a risky heuristic on the axis already won.
Closed for the morphological roster.** REOPEN = large / molecular / congruent data (where
Goloboff's 50% lives) or ns≤4. The env-gated `TS_B2_CEILING` counter is left in `ts_tbr.cpp`
(byte-identical when unset) as the reopen measurement tool.

**Aggressive-collapse PROTOTYPE built 2026-06-22 (user-directed: "build it so it's available
when we propose recipes for large datasets" — i.e. pre-build the reopen lever).** Not a roster
win (B2 stays closed here) but a SHELF knob for the large/molecular/congruent regime. New
`compute_collapsed_flags_aggressive` (`src/ts_collapsed.cpp`) implements TNT `collapse 3`
(min-length-0): an internal edge is collapsible iff `final[p] & final[c] != 0` for every
character. Criterion validated bit-for-bit vs a brute-force MPR oracle — R formula 0/3192,
kernel port 0/206 internal edges, **fp=0 (never over-collapses)** —
`dev/benchmarks/b2_minlength_oracle.R` + `b2_collapsed_kernel_validate.R`. Env-gated
`TS_COLLAPSE_AGGRESSIVE` (default OFF, byte-identical), scoped to `tbr_search` neighbourhood
reduction (pool dedup keeps exact flags); NA falls back to conservative. Gates: exercised
(18–22 internal collapses on random Zanol vs conservative 0), score==full_rescore (exact),
no quality regression (Zanol 1261 / Dikow 1606). Handoff for #40:
`dev/plans/2026-06-22-collapse-aggressive-strategy-briefing.md` (score on large/congruent
corpus, anytime-curve speed metric, MANDATORY quality-not-regressed gate — the heuristic can
skip improving moves). Uncommitted in worktree.

### B3. Rooted representation of an inherently unrooted problem (+ a coupling to B1)

**What we committed to.** A designated root pseudo-node at `n_tip`, `parent[root]=root`
(`ts_tree.h:9`), rooting fixed on tip 0 (`build_postorder.h:37-39`). Fitch parsimony is
root-invariant, so an unrooted problem is solved on a rooted structure. This forces a whole
sub-architecture: directional up-messages / virtual-root state sets
(`ts_fitch.cpp`, `fitch_indirect_length_cached`), a two-phase TBR (rooted inner loop over
2n-4 edges + explicit root-edge enumeration `try_root_edge_moves`, `ts_tbr.cpp:2135-2151`),
and historically a real soundness bug — the directional-vroot under-count, fixed at
cpp-search `2b299e4b` (memory: [tbr-rooted-vs-unrooted]).

**Status: bug fixed, structural cost remains.** Don't re-open the bug; the reframe is that the
representation *permanently carries* the vroot bookkeeping + root-edge enumeration tax, and it
was the root cause of a class of completeness bugs. A natively-unrooted formulation is a
plausible-but-large alternative — flagged, not recommended.

**The coupling worth elevating (B3 × B1).** The unrooted-TBR *completeness certification* is
**gated OFF whenever sector / constraint / tabu / pool is active** (`ts_tbr.h:39-44` /
`ts_tbr.cpp:1384`: `do_reroot` false under those modes). The gate is a genuine *correctness*
incompatibility, not conservatism — the completeness mechanism's outer loop physically
re-roots, permuting node ids and invalidating the **node-id-keyed `sector_mask`**. So our
highest-value mode, sectorial search, never gets the unrooted-TBR-optimum guarantee, and
forcing it on is real plumbing (mask remap), not a flag flip.

**REFUTED as a mission lever 2026-06-22** (falsification probe
`dev/benchmarks/b3b1_endstate_probe.R`, result `b3b1_endstate_result.rds`). The decisive
end-state test: of the 28/30 INDEP single-replicate Zanol2014 end-state trees that MISS the
1261 optimum, how many have an improving unrooted-TBR neighbour (run the kernel's *validated*
unrooted-complete path — the one that scored 0/60 on the completeness oracle — on each)?
**Answer: 0/28** (and 0/28 even SPR-improvable; regime check 0.000 drift, so these are valid
trees not phantom witnesses). The misses are already TRUE unrooted-TBR optima sitting in the
wrong basin — completeness cannot rescue a locally-complete tree. So the sector gate is **not
costing final quality** on the one open quality config; the residual reliable-1261 gap is
basin/budget (corroborates B1 + KPI), not a completeness ceiling. Confirms all three priors
(main-path completeness "does not close the level gap"; good-Wagner-starts "already reach
0-improving"; B1 budget attribution). Scope caveat: tests Zanol end-state only — sufficient
to refute the *mission* lever (Zanol-1261 is the only open quality item), not a claim the
gate is globally costless. Do NOT build the sector_mask-remap plumbing.

---

## C. Throughput decisions — bound speed within fixed reachability

### C1. Eager edge-set precompute vs incremental reinsertion length

**What we committed to.** Per clip we eagerly build full directional state views for the
insertion edges (`compute_insertion_edge_sets`, ≈30% EW CPU). TNT's quick-TBR derives
reinsertion length incrementally and never materialises the views.

**Status: heavily known and parked.** Memory T-P5n/o/p settled this: TS already implements
quick-TBR's *level-1* amortization; the only unbuilt form (level-2 incremental view
derivation = "lever-b") is deferred and judged large-N/molecular-only (the within-clip up-pass
`up[D]=combine(up[parent],prelim[sib])` is already at the one-combine-per-node floor). **No new
discovery — recorded here only so the audit names it as the architectural commitment behind the
single biggest hotspot,** with the reopen condition being a data-class shift (large N, denser,
ns≤4), not a code lever on the current roster.

### C2. Node-major transposed bitset, sized to applicable-state count

**What we committed to.** State sets are node-major, transposed (bit *i* of state-word *j* =
"character *i* can be in state *j*"), `n_states` words per block, ≤64 chars/block grouped by
shared weight (`ts_data.h`, `ts_data.cpp:133-168`). At Zanol's ns=9 this is 4 blocks × 9 words.

**Status: known, analytically at-limit** (memory T-P5r): the transposed layout is already
bit-dense; states-per-word packing would serialise patterns and is strictly worse at high
pattern counts. Residual ~2× on heavy multistate is an **accepted constant factor with an
unpinned mechanism**. Reopen only at ns≤4 (binary/DNA), where a scalar path might beat AVX2.

### C3. Inapplicable (NA) characters as a bolted-on scalar three-pass

**What we committed to.** NA characters use a separate De-Laet-style three-pass
(`down2`, `subtree_actives`) that is irreducibly per-character scalar and cannot share the
SIMD scorer; it also forces a full `exact_verify_sweep` at TBR convergence instead of the
cheap root-edge check (`ts_tbr.cpp:2135-2147`).

**Status: known (na-directional thread).** An exact O(1)/char directional NA scorer was built
and validated bit-for-bit but **refuted on perf** (24-89× slower; memory [na-directional-regimeC]).
The architectural point that remains: NA is a second-class scoring path, so any NA-heavy dataset
pays a structural penalty the EW/IW path avoids. Not currently actionable; recorded for completeness.

---

## D. Checked and ruled out (so they don't resurface as "discoveries")

- **R/C++ marshalling overhead** — clean. Edge matrices convert once in / once out; the entire
  search stays in C++. Not a leak into the hot path (A).
- **Pattern "expansion by weight"** — the header comment is stale: in the live build path each
  unique pattern is one block character and `weight` is a per-block *multiplier*
  (`blk.weight`, score = weight × popcount, `ts_data.cpp:153-166`). Representation size is
  **independent of weights** for both EW and IW — no physical row duplication. Ruled out.
- **Pointer vs array topology** — moot; see A (cost is state recompute, not topology editing).

---

## E. Recommendation — what to actually question

Ranked by (consequence to the mission) × (how under-examined it is):

1. **B1 (population vs single-tree)** — ~~run the basin-diversity instrumentation first~~
   **DONE + SETTLED 2026-06-22.** Pairwise-fuse A/B + INDEP basin diversity on the 4 hardest
   roster datasets: recombination adds no quality, multistart saturates the TNT optimum
   everywhere; residual is budget/per-restart-reach. Handed to #40 as a default-OFF
   efficiency/reliability strategy variant for the full corpus.
2. **B3×B1 coupling** — ~~concrete testable quality risk~~ **REFUTED 2026-06-22** (0/28 missed
   Zanol end-state trees have an improving unrooted-TBR move; see B3 above). The sector
   completeness gate is not costing final quality; do not build the remap plumbing.
3. **B2 (collapsed-space search)** — ~~the only open architectural item~~ **CLOSED for the
   morphological roster 2026-06-22 (see B2 above).** Its "main win" is already built and active
   and captures 0% (free-lunch lever dead two ways); the only remaining lever is the aggressive
   TNT collapse, a heuristic whose speedup on closed-quality data *is* the reachability risk, on
   data the literature says sees no gain, against a bottleneck (per-restart-reach) it doesn't
   touch. Reopen = large/molecular/congruent or ns≤4.
4. **C1/C2/C3** — leave parked; reopen only on the recorded data-class triggers (large N,
   molecular, ns≤4, NA-heavy).

**Frontier status after this audit:** every architectural item is now resolved — A corrected,
B1 settled, B3×B1 refuted, **B2 closed**, C1/C2/C3 parked-with-measurement. No reachability or
element-efficiency probe remains open on the morphological roster. The live frontier is the
composition/orchestration layer (#40), where the addressable wall-clock actually lives. (Caveat:
B2 closure is data-class-scoped — a large/molecular/congruent corpus would reopen it as a real
efficiency lever, since that is where Goloboff's 50% collapse win is documented.)
