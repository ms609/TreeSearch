# MaximizeParsimony search switches — reference for composition (#40)

**Audience:** the agent composing dataset-tailored recipes (#40).
**Scope:** equal-weights Fitch parsimony (the mission objective; `m[m=="-"]<-"?"`).
NA/inapplicable and IW/XPIWE/profile knobs are listed but flagged out-of-mission.
**Authoritative source:** `R/SearchControl.R` (params + defaults), `R/MaximizeParsimony.R`
(presets + `.AutoStrategy`). Defaults below are the *formal* `SearchControl()` defaults
as of cpp-search `1284bdf2` (2026-06-21).

**The one framing that should drive #40** (from the 2026-06-21 KPI, `dev/profiling/kpi-2026-06-21.md`):
quality is CLOSED (TS reaches the optimum ≥ TNT on every mission dataset; on Zanol TS is
the *only* reliably-1261 config). Every component is measured **at-limit**. So the wall gap
is **not per-component throughput — it is budget/composition**: the eye-catching 8–110×
KPI ratio is a *default-budget mismatch* (TS `default` runs a heavy search; TNT `xmult`
default runs a light one). #40's job is to spend the right amount of the right effort per
dataset class — **not** to make any single component faster. Recompose from scratch if any
step's cost changes.

---

## 0. Mission dataset roster (what "per class" means here)

| Dataset | tips | n_states | landscape | notes |
|---|---|---|---|---|
| Wortley2006 | 37 | mixed | flat-ish | small; `sprint`/`default` reach 480 in seconds |
| Giles2015 | 78 | **mixed** | structured | col-reduce helps (17%); reaches 670 |
| Zhu2013 | 75 | **mixed** | structured | col-reduce helps (9%); reaches 624 |
| Zanol2014 | 74 | **uniform ns=9** | hard/structured | the load-bearing case; reaches 1261; col-reduce ~0% |
| Dikow2009 | 88 | mixed | structured | roster max tips (still < 120 → never `large`) |

**No mission dataset is ≥120 tips** — so the `large` preset and `pruneReinsert` auto-enable
(nTip≥120) never fire on the mission roster. The mixed-vs-uniform-`n_states` split is the
key discriminator for `TS_SECT_COLREDUCE` (below).

---

## 1. Top-level `MaximizeParsimony()` arguments (budget + objective)

| Arg | Default | When relevant to #40 |
|---|---|---|
| `strategy` | `"auto"` | The starting point. `auto` → `.AutoStrategy(nTip,nChar)` (§2). #40 will likely **override per class** rather than trust auto. |
| `control` | `SearchControl()` | The expert knob bag (§3). Pass a tuned `SearchControl(...)` here. `...` args to `MaximizeParsimony` also forward into the control. |
| `maxReplicates` | `96L` | **The dominant budget lever.** Wall ≈ replicates × per-rep cost. Most of the "8–110×" is here: TS keeps searching long after the optimum is hit. Pair with a stop criterion (`targetHits`, `consensusStableReps`, `perturbStopFactor`) so it *stops* once converged. |
| `targetHits` | `NULL` | Stop after the best score is independently re-found this many times. **The cleanest convergence stop** — set it (e.g. 3–10) to avoid burning budget post-optimum. Interacts with `perturbStopFactor`. |
| `maxSeconds` | `0` (off) | Wall cap. `0`=use replicate budget. For race-style/time-matched composition, set this; reserves `enumTimeFraction` (10%) for MPT enumeration. |
| `concavity` | `Inf` (EW) | `Inf` = equal weights = **the mission**. Finite = IW; `extended_iw`/`xpiwe_*` only matter then. Leave `Inf`. |
| `inapplicable` | `"bgs"` | NA handling — **out of mission** (EW converts `-`→`?`). Ignore. |
| `nThreads` | `1L` | Parallel replicates. >1 speeds wall but (a) RNG/repro differs, (b) a pre-existing NA parallel crash exists (nThreads≥2 on the NA path). For EW mission timing keep `1L` unless deliberately testing throughput. |
| `verbosity` | `1L` | `0` silent; `1` per-phase (production); `≥2` adds per-phase `score_tree` prints (measurable overhead — Phase-0 finding). Keep `0/1` for timing. |
| `tree` | — | Seed a start tree (shared-start races). |
| `constraint` | — | Topological constraints; clears `consensusConstrain`. Out of mission unless asked. |

---

## 2. Strategy presets (the starting recipes) + auto-selection

`.AutoStrategy(nTip, nChar)`:
- `nTip ≤ 30` → **sprint**
- `nChar < 100` → **default** (flat landscape; thorough is pointless — 0/6 benefited)
- `nTip ≥ 120` → **large** (scaled big-tree preset; **never fires on the mission roster**)
- `nTip ≥ 65` (and nChar ≥ 100) → **thorough**
- else → **default**

| preset | ratchetCycles | xss/rss/css rounds | sectorMax | wagnerStarts | outerCycles/resets | fuseInterval | extras |
|---|---|---|---|---|---|---|---|
| **sprint** | 3 | 1/0/0 | 50 | 1 | 1/0 | 5 | tabu off; light — `nTip≤30` |
| **default** | **6** (was 12, T-P5d) | 3/1/0 | 50 | 3 | 1/**2** | 3 | `adaptiveLevel=TRUE` |
| **thorough** | 20 | 5/3/2 | 80 | 3 | 2/3 | 2 (acceptEqual) | `ratchetAdaptive`, `adaptiveStart`, ratchetMode=2 |
| **intensive** | 20 | 5/3/2 | 80 | **5** | 2/3 | 2 | opt-in only; +Wagner starts for hardest datasets (±1 tradeoff) |
| **large** | **12** (kept, T-179) | 3/2/1 | 100 | 1 (biased) | 1/0 | 3 | `annealCycles=1`, `pruneReinsertCycles=5`+NNI, biased Wagner; **≥120t only** |

**Composition note:** the mission roster (37–88t, ≥100 patterns for the hard ones) auto-selects
**`default`** (Wortley, few chars) or **`thorough`** (Giles/Zhu/Zanol/Dikow). The proven headroom
is `thorough → default` on Zhu/Giles (same score, ~2× less wall) — but on Zanol the thoroughness
is **load-bearing** for the reliable 1261. #40's core question: *how far below `default` can each
class go without losing the reliable optimum?*

---

## 3. `SearchControl()` switches, by component — with #40 relevance

### 3a. TBR core
| switch | default | relevance |
|---|---|---|
| `tbrMaxHits` | `1L` | Equal-score trees held per TBR pass. `1`=fastest descent; thorough uses 3 (more plateau capture, slower). Raise only when MPT diversity matters. |
| `clipOrder` | `0L` (random) | **MEASURED (jobs `17533071`@20-rep + `17541277`@40-rep, 3 seeds) = a per-class TRADEOFF, NOT a safe global win.** `2L` (tips-first) is ~1.25× faster overall / ~26% fewer candidates, but it biases the *trajectory* (not byte-identical): **clean win on Zanol-class** (uniform ns=9 — 3/3 reach 1261, consistently ~1.5× faster at 40-rep); **quality tradeoff on Zhu** (loses +1 on 1 seed even at 2× budget — doubling reps did NOT recover it ⇒ not a budget artifact); **wall unstable on Giles** (one seed 60% *more* candidates). So enable `2L` **per-class on Zanol-type data**; do NOT apply blindly. Complements `TS_SECT_COLREDUCE` — clipOrder helps the uniform-ns case col-reduce can't, and hurts the mixed-state case col-reduce helps. N=3 ⇒ directional. |
| `tabuSize` | `100L` | TBR plateau tabu list. `0`=off (sprint). Larger = more plateau exploration, more memory. Marginal for EW; leave at preset. |

### 3b. Starting trees (Wagner)
| switch | default | relevance |
|---|---|---|
| `wagnerStarts` | `1L` | Independent random-addition starts per replicate. `default`/`thorough`=3. `intensive`=5 helped the *hardest* datasets (Wortley −3, Zhu −2) but +1 on Zanol/Giles → **per-class**, not global. More starts = more basin diversity = more wall. |
| `wagnerBias` | `0L` (random) | `1`=Goloboff non-ambiguous priority, `2`=entropy. `large` uses `1` (near-optimal Wagner at 180t, saves restarts). For mission sizes random is fine; bias mainly pays at large t. |
| `wagnerBiasTemp` | `0.3` | Softmax selectivity for biased addition. Only matters if `wagnerBias>0`. |
| `nniFirst` | `TRUE` | NNI pass before SPR/TBR. Negligible ≤88t; **accelerates the Wagner descent at ≥100t**. Keep TRUE. |
| `sprFirst` | `FALSE` | SPR before TBR. Off-default; washed by TBR; benign. Leave FALSE. |
| `adaptiveStart` | `FALSE` | Thompson-sampling over start strategies. `thorough`/`intensive` use it. Needs several replicates to learn → **regresses at large-t/low-replicate**; helps multi-rep mid-size. |

### 3c. Ratchet (the load-bearing perturbation; ~60% of full-EW phase wall)
| switch | default | relevance |
|---|---|---|
| `ratchetCycles` | **`6L`** | **The single biggest banked recipe lever** (12→6 = 20–38% wall, 0 quality loss, T-P5d). `large` keeps 12 (big-tree tradeoff). Ratchet is **load-bearing — do NOT drop to 0** except <~30t (truly-off ≠ TNT; gap is structural). #40 may tune per class (6 is provisional; a size grid will refine). |
| `ratchetPerturbProb` | `0.25` | Per-character perturbation prob. The perturbation *space* was NOT swept (the isolated race used production params). A scheme reaching the optimum in fewer cycles is an **open #40 question** (audit #55-adjacent). |
| `ratchetPerturbMode` | `0L` (zero-weight) | `1`=up-weight, `2`=mixed. `thorough`/`large` use `2`. |
| `ratchetPerturbMaxMoves` | `5L` | TBR moves per perturbation (`0`=auto). Short perturbation + many cycles (ratchet design). |
| `ratchetAdaptive` | `FALSE` | Adjust prob by escape rate. `thorough`/`large` ON. |
| `ratchetTaper` | `FALSE` | Taper prob as pool stabilizes (finer late exploration). Untested mission-wide. |
| `stallEscalateFactor` | `1.0` (off) | >1 escalates perturbation on cross-replicate stall (auto-discovers needed strength). A **runtime-adaptive alternative to hand-tuning** per class — worth a #40 trial on Zanol. |
| `adaptiveLevel` | `FALSE` (TRUE in `default`) | Scale ratchet+drift effort by hit rate. |

### 3d. Sectorial (TNT's workhorse; ~30% of full-EW phase wall; 96% of *its* wall is `tbr_search`)
| switch | default | relevance |
|---|---|---|
| `xssRounds` / `rssRounds` / `cssRounds` | `3` / `1` / `0` | Exclusive / random / constrained sectorial rounds. The 3 run in **sequence**, each with its own trailing full-tree TBR → **a consolidation candidate** (T-S6e: fusing the 3 sequential trailing TBRs into one is a recipe redesign, needs broad e2e). Tune counts per class; sectorial is where TNT escapes via a diverse retained set. |
| `xssPartitions` / `cssPartitions` | `4` / `4` | Partitions (must be ≥1 — SIGFPE guard). thorough/large use 6. |
| `sectorMinSize` / `sectorMaxSize` | `6` / `50` | Clade-size window. thorough=80, large=100. TNT uses ~min(n/2,45). Bigger sectors = coarser moves, more per-sector cost. |
| `rasStarts` | `1L` | Per-sector RAS+TBR restarts. `3` (TNT-faithful) **closes the rss-ONLY gap (+7/+8→+1, wins time-matched)** but is **REDUNDANT in the full thorough pipeline** at mission sizes (Zanol/Zhu reach optimum at `1`, 60s). **Revisit for larger datasets / shorter budgets** where the full search can't converge. |
| `sectorAcceptEqual` | `FALSE` | Accept equal-score sector resolutions (plateau walking). For flat/NA landscapes; gated out at plateau for EW mission (#24). |
| `sectorMaxHits` | `1L` | Equal-length trees the inner sector TBR holds. Pairs with `sectorAcceptEqual`. |
| `sectorCollapseTarget` | `0L` (off) | Collapse a big sector into ~this many composite terminals (coarse skeleton; Goloboff-1999 reduced dataset). The **row-axis** reduction (cf. `TS_SECT_COLREDUCE` = column-axis, §4). Worth pairing for large sectors. |
| `postRatchetSectorial` | `FALSE` | Re-run XSS+RSS+CSS after ratchet (TNT-interleaved). Adds a full sectorial pass per ratchet — expensive; only if it earns score. |

### 3e. Drift / annealing / NNI-perturb (alternative escapes)
| switch | default | relevance |
|---|---|---|
| `driftCycles` | `0L` (off) | Tree drifting (Goloboff). #25 added TNT-faithful drift for the +1 datasets. Off in all EW mission presets — a per-class escape to trial. `driftAfdLimit`/`driftRfdLimit` bound accepted suboptimal moves. |
| `nniPerturbCycles` | `0L` (off) | NNI-topology perturbation (complements weight-ratchet). **Measured 69% overhead, zero time-adjusted benefit (T-274)** — leave off unless a class proves otherwise. `nniPerturbFraction=0.5`. |
| `annealCycles` | `0L` (off) | PCSA simulated-annealing perturbation. Effective ≥100t (the `large` preset uses 1 cycle to replace drift). Below mission sizes, unproven. `annealPhases/TStart/TEnd/MovesPerPhase` shape the schedule. |

### 3f. Prune-reinsert (T-266; a strong large-tree perturbation)
| switch | default | relevance |
|---|---|---|
| `pruneReinsertCycles` | `0L` (off) | Drop tips → restructure backbone → reinsert. **Auto-on only via `large` (nTip≥120) — never on the mission roster.** An **exact-scorer port is prepared** (worktree `claude/scoreapprox-probe 41b0d237`): `expand_and_reinsert` used a union-of-finals approximation (the #27 miss); the exact `fitch_indirect_length_cached` port is validated (Δ=0) and parked for #40 to land + A/B *if a class ≥120t enters scope*. |
| `pruneReinsertDrop` | `0.10` | Fraction dropped/cycle (≥3 tips, keep ≥4). |
| `pruneReinsertSelection` | `0L` (random) | `1`=instability, `2`=missing-data, `3`=combined tip selection. |
| `pruneReinsertTbrMoves` / `FullMoves` | `5` / `0` | Backbone / full-polish TBR budgets. |
| `pruneReinsertNni` | `FALSE` | NNI full-polish instead of TBR — **~5× faster at ≥120t**; `large` uses it (TBR polish was catastrophic at 206t/60s). |

### 3g. Pool / fusing
| switch | default | relevance |
|---|---|---|
| `fuseInterval` | `3L` | Fuse pool trees every n reps. **AUDIT #55: fuse is DEAD WEIGHT on the >64t mission class** — fires (7×/run, ~70 exchanges) but **0 improvements** on Zanol/Zhu/Giles. Already effectively off (`poolSuboptimal=0`→pool size 1→fuse SKIPPED). **#40: drop fuse / raise `fuseInterval` on this class to reclaim the inter-replicate `score_tree`+`tree_fuse` overhead.** (Fuse may still matter at >88t / with a diverse pool — re-measure if the class changes.) |
| `intraFuse` | `FALSE` | Within-replicate fuse vs pool donors. Also 0 improvements in the #55 probe. |
| `fuseAcceptEqual` | `FALSE` | Accept equal-score fused trees (plateau). thorough/large ON. |
| `poolMaxSize` | `100L` | Max pool trees (≥1; segfault guard). |
| `poolSuboptimal` | `0` | Retain trees within N steps of best. **`0` ⇒ pool collapses to optimal-only ⇒ fuse rarely fires.** Raise (e.g. 5) ONLY if you want fuse/diversity to actually run — but #55 says fuse doesn't pay here, and a diverse pool is TNT's fuse fuel (untested as a TS lever). |

### 3h. Stopping / outer loop (the budget governors — high-leverage for #40)
| switch | default | relevance |
|---|---|---|
| `consensusStableReps` | `0L` (off) | Stop when the strict consensus is unchanged for N reps. **A convergence stop** (3–5 typical). Pair with `targetHits`; stops at whichever fires first. |
| `perturbStopFactor` | `2L` | Patience: stop after consecutive non-improving reps exceed `(targetHits/hits)·nTip·factor`. Scales patience with progress. `0`=off. **The main "don't over-search" governor** — tune per class. |
| `outerCycles` | `1L` | Repeats of [XSS/RSS/CSS→ratchet→NNI→drift→TBR] per replicate. thorough/intensive=2. More = TNT-style interleaving, more wall. |
| `maxOuterResets` | `0L` | Improvement-triggered resets of the outer counter (`-1`=unlimited). `default`=2, thorough=3. Lets a productive replicate keep going. |
| `enumTimeFraction` | `0.1` | Fraction of `maxSeconds` reserved for MPT enumeration. `0`=disable reserve. Only matters with `maxSeconds>0`. |
| `consensusConstrain` | `FALSE` | Lock pool-consensus splits as constraints after ≥5 reps (focus on uncertain regions). Off-default; only when no user constraint. |

---

## 4. Opt-in env-var levers (not in `SearchControl`)

| env | default | relevance |
|---|---|---|
| `TS_SECT_COLREDUCE` | unset (off) | **AUDIT #56, shipped `830b8cc3`.** Per-sector column-axis reduction: drops chars constant-within-{sector tips+HTU} + repacks → smaller inner-sector block scan. **Saving: Giles 17%, Zhu 9%, Zanol ~0%** (uniform ns=9 = least reduction = the load-bearing case). **Enable per-class on MIXED-`n_states` data (Giles/Zhu/Dikow); skip on Zanol.** Validated bit-exact (dScore=0 9/9, valgrind clean) but **changes the search trajectory on mixed-state data** (`dCand≠0`, equally-optimal different path) ⇒ **opt-in, NOT a default flip.** **Before enabling by default for any class, run a sector-score ORACLE** (reduced vs full score, same topology, mixed-state) — an accept-gated search can't discriminate a masked packing bug. Read once at static init ⇒ set in the env **before** the R process starts; one process per arm. |

(`TS_AUDIT_PROBE` is a *compile* flag for measurement counters — not a runtime recipe lever.)

---

## 5. Composition cheat-sheet — session-derived starting recommendations

These are **hypotheses for #40 to validate**, not settled recipes:

- **Small / few-char (Wortley-class, <100 patterns or ≤30t):** `default` (or `sprint` ≤30t). Thorough is pointless on flat landscapes (0/6 benefited). Low `maxReplicates` + `targetHits` stop early.
- **Mixed-`n_states` structured (Giles/Zhu/Dikow-class, 65–88t):** `thorough`-ish but trim toward `default`; **enable `TS_SECT_COLREDUCE`** (9–17% sectorial); ratchet 6; **drop fuse**; add a `targetHits`/`perturbStopFactor` stop. `rasStarts=1` (full pipeline converges).
- **Uniform ns=9 hard (Zanol-class, ~74t):** thoroughness is **load-bearing** for the reliable 1261 — do NOT strip ratchet/sectorial aggressively. `TS_SECT_COLREDUCE` gives ~0% here (skip). **`clipOrder=2L` IS a clean ~1.5× throughput win here (measured, 3/3 optima)** — the one place it's safe. The win is otherwise *stopping at the right time*, not running lighter. Consider `stallEscalateFactor>1` to auto-find the perturbation strength.
- **Large (≥120t, not in mission roster):** `large` preset; `pruneReinsert`+NNI; biased Wagner; anneal; **land the prepared exact-scorer port** (41b0d237) + A/B; `rasStarts=3` may re-enter (short-budget/large).
- **Cross-cutting cheap trials:** `clipOrder=2L` is **measured = per-class, Zanol-only safe** (see §3a — ~1.5× on Zanol, but +1 quality cost on Zhu at 2× budget), NOT a global flip; consolidate the 3× sequential trailing sectorial TBRs (T-S6e); the redundant trailing TBRs generally.

## 6. Hard constraints for #40

- **Quality first:** any recipe must still reach the class's known optimum (Wortley 480 / Giles 670 / Zhu 624 / Zanol 1261 / Dikow's best) at full budget. Speed that loses the reliable optimum is a regression, not a win.
- **Recompose from scratch if any step cost changes** (the reason composition waits until all pieces are final).
- **Validate with separate processes per arm** for any env-flag (static-init read) and prefer **candidates_evaluated / score** equality over wall for correctness; wall is the *saving* axis.
- Ratchet is load-bearing (don't zero it >30t); fuse is dead weight on the mission class (do zero it); `TS_SECT_COLREDUCE` is mixed-state-only and never a default.

---
*Maintained by the component-isolation/audit workstream. Companion: `dev/plans/2026-06-19-component-isolation-profiling.md` (component verdicts), `dev/profiling/kpi-2026-06-21.md` (the gap reframe).*
