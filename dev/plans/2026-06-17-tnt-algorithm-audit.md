# TNT 1.6 "new technology" search vs TreeSearch cpp-search — line-level audit

Date: 2026-06-17. Branch: `cpp-search`. READ-ONLY audit (no `src/` edits).
Goal: map TNT's xmult / sectsch / ratchet / drift / fuse to our exact code,
flag divergences, rank them by likelihood of explaining TNT's *escape*
advantage on EW-Fitch, and confirm each phase fires.

Primary sources consulted (all load-bearing):
- TNT 1.6 console help: `help xmult/sectsch/ratchet/drift/tfuse/rseed/mult/bbreak`
  -> `dev/benchmarks/tnt_help.txt`
- TNT 1.6 **default settings dump** (`sectsch:; xmult:; ratchet:; drift:; tfuse:;`)
  -> `dev/benchmarks/tnt_defaults.txt`  (THE authoritative default values)
- Goloboff 1999, *Analyzing Large Data Sets in Reasonable Times* (Cladistics
  15:415-428) — SS/RSS/CSS/MSS, tree-fusing, tree-drifting. Pages 417-422 + App 1.
  Local PDF: `C:/Users/pjjg18/Zotero/storage/TETHI9A5/`.
- Nixon 1999 (ratchet) — referenced; algorithm is weight-perturbation + TBR.
- Prior reverse-engineering: memory `tnt-sectorial-recipe.md`, `sector-resolve-status.md`,
  `untried-search-ideas.md`; `dev/expertise/tnt.md`.

The EW-Fitch comparison is triggered by replacing `-` with `?` (no inapplicables),
so **all parsimony scoring is plain Fitch**. This fact is decisive (see §3).

---

## TNT defaults (from `tnt_defaults.txt`) — the reference behaviour

```
Sectorial (sectsch):
  * Separate matrix-buffer for sectors (xbuf ON)
  * Random sector selections; min/max size 0 (=> size law min(n/2,45)/max(45,n/10))
  * Sectors of size BELOW 75 analyzed with 3 RAS+TBR
       (+ EXTRA 3 starts if the first 3 produce score differences)   <-- combstarts
  * Global TBR every 10 substitutions (small AND large sectors)
  * NOT accepting equally good subtrees                              <-- equals OFF by default

xmult (Extra search settings):
  * 4 replications as starting point for each hit
  * Each replication autoconstrained (previous + wagner)
  * Each replication: random sectorial searches, NO ratchet,
       WITH drifting (5 iters), NO hybridization, fusing (1 round)   <-- drift+fuse, NOT ratchet
  * 1 hit; multiply trees by fusing after hitting best score

ratchet:   50 iters, 40 subs, equal cycle yes, up/down weight 4    (NOT used by default xmult)
drift:     30 iters, 60 subs, AFD 1, RFD 0.20, reject factor 3.00
tfuse:     5 rounds, start from best, accept-equal OFF, TBR-swap after fusion
```

Goloboff 1999 RSS algorithm (p.418-419), verbatim structure:
> (1) Select a sector at random so the reduced data set has S terminals.
> (2) Do **R replications of RAS + TBR** (saving a single tree) for the reduced
>     data set. If the R replications produce trees of the **same length** (S was
>     in a non-conflictive region) go to (3); otherwise do **r additional** reps.
> (3) **Choose the best among the R + r replications and the present resolution**
>     for the sector and place it in the whole tree.
> (4) Do a round of **global swapping, but only if replacements at (3) have been
>     made more than X times.** Go to (1), N times.
>
> "The best size seems to be 35 to 55 nodes. For that value of S, **R = 3 and r = 3**.
>  A round of global swapping of the entire tree is made **every 5 to 10 replacements**,
>  as that number makes it likely that (through clade substitution) the tree will
>  have become **globally suboptimal under TBR**."

Reduced dataset (p.418 + App 1): internal nodes represented by their **first-pass
state sets**; the basal node (HTU) by the first-pass set calculated upward. The
**HTU is an ordinary terminal** in the reduced RAS+TBR. (App 1 `marknodes`
accumulates a connected sector by clade-size until `<= sector_sz`.)

---

## §1. MAPPING TABLE  (TNT step | TNT behaviour | our code | match | divergence)

| # | TNT step | TNT documented behaviour | Our function : line | Match | Divergence detail |
|---|----------|--------------------------|---------------------|-------|-------------------|
| 1 | Starting trees | `xmult`: 4 RAS replications per hit; addition seq = `ras` (random) | `driven_search` loop `ts_driven.cpp:728`; per-rep Wagner in `run_single_replicate:99-144` (`wagnerStarts` random-order Wagner) | partial | We use 1 replicate stream with `wagnerStarts` (1 sprint / 3 default+thorough) Wagner starts kept-best, then iterate `maxReplicates`. TNT does 4 *independent* RAS per hit then fuses them. We fuse across the replicate pool periodically (`fuseInterval`). Functionally similar; not identical. |
| 2 | Wagner addition | RAS = random addition sequence + greedy placement | `random_wagner_tree` (ts_wagner.cpp) | yes | Matches. (Biased Wagner is an extra TreeSearch option, off by default for EW < 120t.) |
| 3 | Initial hill-climb | RAS+TBR (or SPR then TBR) to local optimum | `run_single_replicate:165-173`: optional NNI/SPR warmup then `tbr_search` to convergence | yes | Matches; `nniFirst=TRUE` adds an NNI warmup TNT lacks but that is only a speed lever. |
| 4 | Sector ORDER in xmult | "css first, then rss, then xss last"; ratchet/drift always FOLLOW sectorial | `run_single_replicate:231-303` runs **XSS, then RSS, then CSS**; ratchet/drift after | **partial** | Order is reversed (we do XSS->RSS->CSS; TNT does CSS->RSS->XSS). Low impact, but note default xmult uses **RSS only** (css/xss off by default). |
| 5 | Sector SELECTION (RSS) | Random connected sector, S terminals, **size law min(n/2,45)**; reduced data set built from first-pass state sets | `rss_search:993-999` selects **eligible internal nodes whose clade size ∈ [min,max]** (a CLADE, not a constructed connected region); `build_reduced_dataset:273` | **partial** | We SELECT an existing clade in a size band ([6,50] default / [6,80] thorough). TNT CONSTRUCTS a connected sector of ~n/2 by accumulating clades (App 1 `marknodes`). Memory's band-shape probe found this lever **real but marginal-to-noise** end-to-end. |
| 6 | Reduced-dataset HTU | Internal nodes = first-pass state sets; basal HTU = first-pass set calculated upward; **HTU is an ordinary terminal** | `build_reduced_dataset:462-475` sets HTU tip-state = `compute_from_above_for_sector` (exact rest-of-tree first-pass set). For EW-Fitch this is **exact** (reduced length = full length − const). | **partial** | Scoring fidelity matches and is EXACT for EW (verified §3). BUT we anchor the HTU at the synthetic root and **freeze it there** — TNT lets it float as a terminal. This is divergence D1 (§2). |
| 7 | Within-sector search | **3 RAS+TBR** (+3 more if scores differ); accept best of R+r + present resolution | `search_sector:773-847`; start 0 = TBR on existing subtree, starts 1.. = `build_ras_sector` RAS rebuild | **partial/no** | Default `rasStarts=1` => **only the TBR-polish of the existing subtree, NO RAS rebuild**. ALL TreeSearch presets leave `rasStarts=1` (never overridden). TNT's default is 3+3. Divergence D3a. Even when `rasStarts=3`, the RAS rebuild keeps the HTU frozen (D1). |
| 8 | Sector ACCEPT | "if a better configuration is found, it is replaced on the whole tree" (best reduced score); full tree may become "globally suboptimal", cleaned by periodic global TBR | `search_sector:829-831` (keep best REDUCED score) THEN `rss_search:1069-1148`: reinsert, full-rescore, **revert unless full tree STRICTLY improves** | **partial** | Double gate: (a) reduced strict-best, (b) full-tree strict-improve revert. Goloboff accepts on the reduced score and tolerates transient full-tree suboptimality. **HOWEVER for EW-Fitch the from-above HTU makes reduced-improve ⟺ full-improve, so gate (b) is a NULL divergence here** (verified §3, 0 gate-bites). Bites only under inexact NA scoring. |
| 9 | Global TBR cadence | Round of global TBR **every 5-10 sector replacements** | `rss_search:1154-1165` runs **one** global TBR at the END of all `rss_picks` | partial | We TBR once per RSS round (and `rssRounds` loops the whole thing). Coarser cadence than TNT's every-5-10. Minor; the per-round outer loop approximates it. |
| 10 | accept-equal subtrees | `equals` OFF by default | `sectorAcceptEqual=FALSE` default; `search_sector` `accept_equal` param | yes | Matches TNT default. (Prior probes: accept_equal "wanders", undirected; not the lever.) |
| 11 | Drift (DFT) | xmult default: **5 drift iters per replication**; suboptimal accepted by RFD = (F−C)/F rejected when > Z = X/(F+J−C); AFD limit; alternate perturbed/unperturbed | `drift_search` (ts_drift.cpp:721); `drift_phase` AFD/RFD logic at :578-695 | yes | Algorithm faithfully implemented (AFD `driftAfdLimit`, RFD `driftRfdLimit`). BUT **all TreeSearch presets set `driftCycles=0`** — drift never runs by default. TNT's xmult default leans on drift. Divergence D3b. |
| 12 | Ratchet | NOT in default xmult; standalone: 50 iters, perturb (up/down weight 4%), short TBR, restore weights, full TBR | `ratchet_search` (ts_ratchet.cpp:136); perturb modes zero/upweight/mixed; initial TBR + cycles | yes | Algorithm matches Nixon 1999 / TNT ratchet. BUT TreeSearch makes ratchet the **primary** escape (12-20 cycles in every preset), the inverse of TNT's xmult default (ratchet off). This is *how TreeSearch approaches TNT* but is a different engine balance. |
| 13 | Tree-fusing | Exchange shared clades (>=5 taxa) between trees; accept improving (equal off by default); TBR after; 5 rounds | `tree_fuse` (ts_fuse.cpp:325): shared-split detection, `replace_subtree`, accept `< score`, TBR cleanup, `max_rounds` | yes | Faithful. `fuseAcceptEqual` off by default (matches). Pool-level fusing every `fuseInterval` reps + optional intra-rep fuse. Reroot-segfault fix on >64t noted in memory NOT yet in cpp-search (`fuse-reroot-segfault.md`) — orthogonal to the gap. |
| 14 | CSS | Consensus-based sector selection; sector-restricted TBR, exact scoring | `css_search` (ts_sector.cpp:1330): `xss_partition` + `tbr_search` with `sector_mask` | partial | We partition (XSS-style) rather than select from a consensus of conflict; the TBR is exact (no HTU). BUT the mask **freezes the sector's attachment** (`ts_tbr.cpp:807,910,1130`): clips and regrafts confined to the clade, so CSS-TBR cannot relocate the whole sector. Divergence D1 (CSS variant). |
| 15 | XSS | N exclusive non-overlapping sectors tiling the whole tree, R rounds, global TBR after each round | `xss_search` (ts_sector.cpp:1172) + `xss_partition:905` | yes | Partitioning + per-sector reduced-dataset rebuild matches. Same HTU-freeze (D1) as RSS since it uses `build_reduced_dataset`+`search_sector`. |
| 16 | Stop rule | `hits N` (default 1); consense stabilization options | `driven_search`: `targetHits`, `consensusStableReps`, `perturbStopFactor` | yes | Matches in spirit; TreeSearch adds time-budget + MPT enumeration tail. |

---

## §2. THE KEY QUESTION — what lets TNT's sectorial cross barriers ours cannot

Setup the analysis must respect: the memory's **shared-start probe is the clean
signal**. From an *identical* TNT `mult` T0, TNT's RSS sectorial improves +3..+11
on the gap datasets (Zanol +11, Wortley +7, Zhu +4, Giles +3); **ours improves 0
on all 6**, having examined up to 1.26M candidates. T0 is already OUR global-TBR
optimum. So the gap is a genuine SECTORIAL quality/escape gap, independent of
wall-clock, and independent of the starting tree. Every prior refutation
(`sector-resolve-status.md`, `tnt-sectorial-recipe.md`) varied the sector REBUILD
or SCORING while holding two things fixed: (i) the strict full-tree accept gate,
and (ii) the **frozen sector↔rest-of-tree attachment**.

### Candidate (a) — RAS multi-start banking best topology. REFUTED (prior).
`build_ras_sector` + `search_sector(ras_starts)` implemented; shared-start gap
UNCHANGED, 10/12 cases `ts_sect == start` (`sectorial_shared_greedy.csv`). Rebuild
*alone* is null.

### Candidate (b) — accept equal-length sector solutions. REFUTED (prior).
`accept_equal` walks the plateau but undirected; drifts AWAY from the target basin
(`sector-resolve-status.md` §2). Matches TNT default anyway (`equals` OFF).

### Candidate (c) — reduced-dataset scoring is approximate where TNT's is exact.
REFUTED for EW, and now **verified empirically** (this audit, §3): under EW-Fitch
the from-above HTU gives reduced-improve ⟺ full-improve (0 gate-bites). Exact-CSS
probe also null (`sectorial_shared_css.csv`). So scoring fidelity is NOT the lever
for the EW case. (It *does* bite under native NA/Brazeau — see §3 — but that is not
the audited comparison.)

### Candidate (d) — sector SELECTION (clade-band vs constructed region). REAL but SMALL.
Band-shape probe: real (band wins 30/40 rss-only, median −6/−7.5) but ratchet
SUBSUMES it; end-to-end collapses to noise, harms Giles/Zanol on some seeds. Not
the +0→+11 lever on its own. Secondary.

### >>> THE LEADING, UNTESTED MECHANISM: frozen sector↔rest-of-tree attachment (D1) <<<

Goloboff's reduced dataset treats the **HTU as an ordinary floating terminal**. His
R RAS+TBR replications therefore choose, *jointly and simultaneously*:
  (i) the internal topology of the sector, AND
  (ii) **which node of the sector is basal/adjacent to the rest of the tree** —
       i.e. where the rest-of-tree (HTU) reattaches.
TBR on the reduced data set can put the HTU terminal anywhere, which corresponds
to **re-rooting the sector relative to the rest of the tree** and simultaneously
re-resolving it. This is a move that crosses an uphill barrier in the FULL tree
in a single accepted step, because the reduced score (exact for EW) drops even as
the global arrangement reorganises.

TreeSearch forecloses this in **three independent places** — the HTU/attachment is
frozen so the search only ever explores *rebuild-with-fixed-attachment* or
*reroot-alone*, never their product:

1. **`search_sector` root-structure revert** — `ts_sector.cpp:808-819`.
   After the internal TBR, if the HTU and `sr_mapped` are no longer the two direct
   children of the synthetic root, the move is **discarded** and the topology reverts
   to the pre-TBR snapshot. Any TBR move that floats the HTU (the very move that
   re-roots the sector against the rest of the tree) is thrown away.

2. **`build_ras_sector` insertion restriction** — `ts_sector.cpp:716-748`.
   The RAS rebuild seeds with HTU at `new_root`, content rooted at `sr_mapped`, and
   restricts candidate edges to the subtree **below `sr_mapped`** ("never a root
   edge ... so the HTU stays anchored at new_root"). The HTU is never an addable/
   movable terminal; the rebuild can only re-resolve the clade with the *same* node
   kept basal. This is exactly why `build_ras_sector` reproduces the start 10/12.

3. **CSS sector mask** — `ts_tbr.cpp:807` (`if (sector_mask && !mask[clip_node]) continue`),
   `:910` and `:1130` (`if (sector_mask && !mask[below]) continue`).
   The mask is the clade only. CSS-TBR can neither clip the `sector_root` edge nor
   regraft outside the clade, so it cannot relocate the sector as a whole — the
   exact attachment is frozen.

Global TBR (`rss_search:1154-1165`, the outer-cycle TBR) DOES cover re-rooting the
sector relative to the rest of the tree — **but alone**, on a *fully-resolved* tree
that has already converged (=0 from T0). The sector rebuild covers re-resolution
**but alone**, with attachment frozen (=0). TNT's RAS+TBR-on-reduced-data does
**both at once** (rebuild × free-attachment). That joint move is the one neither
TreeSearch operator can reach, and it is the only structural lever the prior
sessions never isolated (`sector-resolve-status.md` itself flags "its RAS rebuild
holding/accepting equal trees during the rebuild; next probe = log sectors picked"
— the right instinct, wrong mechanism: it is attachment-freedom, not equal-holding).

"Low root_revert frequency rules this out" is a non-argument: revert only counts
HTU-displacing moves the *polish-TBR proposed*; the RAS rebuild never proposes one
(it is forbidden by construction at :716), so the freeze is baked into construction,
not observable as revert frequency.

---

## §3. EXACTNESS VERIFICATION (settles the accept-gate ranking)

Discriminating trace (`dev/benchmarks/diag_accept_gate_trace.R`, using the
`TS_SECT_DEBUG=1` REprintf at `ts_sector.cpp:1081`): does a sector ever improve on
the reduced score (`red_best < red_cur`) while the full tree does NOT improve
(`full_new >= full_best`)? That is the only condition under which the strict
full-tree accept gate (item 8b) can foreclose a real sectorial improvement.

EW-Fitch, Zanol2014, seed 1, rasStarts=3 (`dev/benchmarks/trace_fitch.txt`):
```
sect[117] red_cur=820 red_best=816 full_new=1331 full_best=1335 STRICT
sect[ 95] red_cur=197 red_best=196 full_new=1330 full_best=1331 STRICT
sect[ 75] red_cur=445 red_best=441 full_new=1326 full_best=1330 STRICT
gate-bites (red improved, full did NOT): 0
```
**Every reduced improvement translated 1:1 to a full-tree improvement.** Confirms
the from-above HTU is EXACT for EW-Fitch, so the strict full-tree accept gate
(item 8b) is a **NULL divergence for the audited case**. (The trace only fires on
the accept branch, so it reports the improving sectors; the point is none of them
needed transient full-tree worsening.) Native NA (`-` kept) is where exactness
breaks and the gate can bite (the `WORSE-revert` branch at `rss_search:1143-1148`
exists precisely for that) — but EW-Fitch is the comparison in scope.

=> The acceptance double-gate is **demoted**: it is the tidy-looking smoking gun
that does not fire for EW. The real lever is the attachment freeze (D1).

---

## §4. PER-PHASE EXECUTION CONFIRMATION (does each phase fire?)

Install: `.agent-audit` per-agent lib (NOT the shared default), built clean.
`bench_phase_yield.R` with `TS_LIB=.agent-audit`, Wortley2006 + Zanol2014, seeds
1-3, 20s, nThreads=1, `strategy="auto"` (Wortley 37t -> "default"; Zanol 74t ->
"thorough"). Phase columns = % of wall-clock.

`dev/benchmarks/phase_yield_audit.csv`; medians over 3 seeds:

| dataset | tips | preset | score_med | reps | late_frac | wagner | init_tbr | sector | ratchet | final_tbr | fuse |
|---------|------|--------|-----------|------|-----------|--------|----------|--------|---------|-----------|------|
| Wortley2006 | 37 | default | 481 | 216 | 0.81 | 6% | 2% | **7%** | **83%** | 2% | 0% |
| Zanol2014 | 74 | thorough | 1264 | 27 | 0.41 | 5% | 5% | **22%** | **66%** | 2% | 0% |

ALL phases fire (wagner, init_tbr, sector, ratchet, final_tbr all > 0). The
distribution is the **mirror image of TNT**:
- TreeSearch: **ratchet-dominated** (83% Wortley / 66% Zanol), sectorial a thin
  slice (7% / 22%), drift absent (presets set `driftCycles=0`), periodic fuse ~0%.
- TNT xmult default: **~67% sectorial**, drift(5)+fuse(1), ratchet **OFF**.

So the phase that consumes TreeSearch's budget (ratchet) is the phase TNT does NOT
run by default, and the phase TNT leans on (sectorial, with floated-HTU RAS+TBR and
drift) is the under-weighted, structurally-limited slice in our presets (D1+D3).
`late_frac` 0.81 on Wortley means 81% of replicates ran AFTER the last improvement
— effort spent re-finding the same optimum, consistent with the 2x wall-clock gap.
Scores (Wortley 481 vs TNT best 479; Zanol 1264 vs TNT best 1262) reproduce the
documented +1..+3 EW gap, so this install behaves as the gap reports describe.

---

## §5. TOP-3 DIVERGENCES (ranked by likelihood of explaining TNT's escape) + experiments

### D1 (RANK 1) — Frozen sector↔rest-of-tree attachment; HTU never floats
- Files/lines: `ts_sector.cpp:808-819` (root-structure revert), `ts_sector.cpp:716-748`
  (RAS insertion restricted below `sr_mapped`), `ts_tbr.cpp:807/910/1130` (CSS mask).
- Hypothesis (falsifiable): TNT's per-sector RAS+TBR treats the HTU as an ordinary
  terminal and so jointly re-resolves the sector AND re-roots it against the rest
  of the tree in one accepted, reduced-score-improving step. TreeSearch can do
  rebuild-alone (=0 from shared T0) or reroot-alone-via-global-TBR (=0), never the
  product. Freeing the HTU will turn the shared-start result from +0 to a
  substantial fraction of TNT's +3..+11.
- Experiment: In `build_ras_sector`, make the HTU a normal addable terminal (allow
  insertions on ALL reduced-tree edges including the synthetic-root edge), and in
  `search_sector` drop the `root_ok` revert (or, equivalently, define reinsertion
  by *whichever sector node ends up adjacent to the HTU terminal* and reattach the
  rest-of-tree there). Re-run `bench_sectorial_shared.R` from the identical TNT T0.
  PREDICT shared-start 0 -> +N on Zanol/Wortley. If still 0, demote and proceed to D2.
  (Cheaper precursor: log, per accepted sector, whether the node adjacent to the HTU
  in the rebuilt sector differs from the original basal node — if it is ALWAYS the
  same, attachment is provably frozen.)

### D2 (RANK 2) — `rasStarts=1` in every preset (vs TNT 3 + 3-on-disagreement)
- Files/lines: `R/MaximizeParsimony.R:106-216` (no preset sets `rasStarts`; stays
  `1L` from `R/SearchControl.R:290`). Engages `search_sector` start-0-only path.
- Hypothesis: even rebuild-alone is null from a converged T0, but rebuild is a
  *precondition* for D1 — you cannot exploit a floated HTU without re-resolving the
  sector. With the HTU frozen, `rasStarts=3` is null (already shown). With the HTU
  floated (D1), `rasStarts>=3` becomes necessary to realise the joint move, and the
  "+3 extra starts on score disagreement" (`combstarts`) matters. So D2 is *coupled*
  to D1: the only meaningful test of D2 is rasStarts>=3 *with* a floating HTU.
- Experiment: factorial on `bench_sectorial_shared.R` — {HTU frozen, HTU floated} ×
  {rasStarts 1, 3, 6}. Expect improvement only in the (floated, >=3) cell.

### D3 (RANK 3) — Engine balance: ratchet-primary vs TNT's sectorial+drift+fuse
- Files/lines: every preset sets `driftCycles=0L` + `ratchetCycles` 12-20
  (`R/MaximizeParsimony.R:107-216`); TNT xmult default = RSS sectorial + drift(5) +
  fuse(1), ratchet OFF (`tnt_defaults.txt`). Plus global-TBR cadence: `rss_search:1154`
  runs one TBR at end vs TNT every-5-10 replacements.
- Hypothesis: TreeSearch reaches near-TNT quality by substituting an aggressive
  ratchet for the sectorial escape it cannot perform (D1). This explains why
  ratchet-off TreeSearch trails by +4..+8 (memory) and why ratchet is "necessary
  but never erases the gap" — ratchet is a *workaround* for the missing sectorial
  move, not an equal. Closing D1 should let drift+fuse (cheaper) replace some
  ratchet load, attacking the 2x wall-clock gap simultaneously.
- Experiment: after D1 lands, A/B a TNT-faithful preset (rss sectorial w/ floated
  HTU + rasStarts=3, drift=5, ratchet=0, fuse every 2) vs current default at matched
  wall-clock on Wortley/Zanol/Zhu/Giles. Predict equal-or-better quality at lower
  time if D1 is the true lever.

---

## §6. Caveats / notes for the orchestrator
- `src/` is READ-ONLY here; all line numbers are against the working tree at audit time.
- The `TS_SECT_DEBUG` trace only prints on the sector *accept* branch
  (`rss_search:1081` is inside `if (accept && sector_best <= sector_current)`), so it
  reports improving sectors, not rejected ones. For a full gate-bite census, the
  orchestrator may want a trace on the reject/`WORSE-revert` path too — but the EW
  exactness argument (reduced=full−const) makes a bite impossible for EW regardless.
- `candidates_evaluated` is NOT a clean cumulative counter (goes negative across
  sector rounds, per memory) — use *score*, not candidate deltas, as the signal.
- Reduced-dataset alloc churn (~19% VTune, memory) and `xbuf` reuse are a wall-clock
  lever (TNT reuses a buffer; we rebuild per sector) — orthogonal to the escape gap.
- TNT run-script filenames must be multi-char alphabetic (`helpdump.run`), else TNT
  parses the basename as a command.
