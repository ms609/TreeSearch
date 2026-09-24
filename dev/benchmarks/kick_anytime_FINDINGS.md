# ratchetPerturbMaxMoves 5 → auto/scaled — corpus ANYTIME test (VERDICT: not a default)

**Date 2026-07-14. Mission A precursor** (`dev/plans/2026-07-14-mission-A-cold-start-basin-capture.md`,
first task). Hamilton array **job 17872154**, current cpp-search engine built fresh into `kick-lib`
(job 17872123). Harness `kick_anytime.R`; analysis `kick_analyze.R`.

## Question
The reweighting kick `ratchetPerturbMaxMoves` moved project5432's floor 1949→1945, but the fixed default
`5L` is well below the engine auto `max(20, min(200, n_tip/8))` for every tree ≳40 tips, and the *reach*
benefit was **5432-only (n=1)**. The open question (memory `project5432-basin-structure`: "wall benefit on
solved matrices UNTESTED"): does 5→auto help — or **regress** — wall-clock **time-to-optimum** broadly?
The suspected regression was on **small/easy** matrices (auto=20 is a near-restart on a 25-tip tree).

## Design (buildless; advisor-reviewed)
- Drive the real production `strategy="auto"` presets; override ONLY `ratchetPerturbMaxMoves`. Three arms:
  **fixed5** = 5 (current default) · **auto** = 0 (engine `max(20,min(200,n/8))`) · **scaled** =
  `max(5, round(n/8))` (the "scale with n_tip" candidate; no 20-floor → gentle on small trees).
- **Metric = per-replicate best-vs-wall ANYTIME trace** via `progressCallback` (improvement events).
  `nThreads=1` (deterministic + serial `candidates_evaluated`). Per (matrix,seed) the target is the
  **union-best final across arms**; report reach, time-to-first-hit, total wall. Production stopping
  (maxReplicates=300, targetHits=auto, tier `maxSeconds` cap).
- Corpus = the fixed 25-matrix **training** sample (7 small / 7 medium / 7 large / 4 xlarge), **validation
  split sequestered** (harness hard-refuses non-training keys). 5 seeds → 125 cells × 3 arms.
- Regime = EW Fitch, gaps→missing (the 5432 regime).

## Result — no broad wall win; mild wall regression where it doesn't help reach

| tier (n matrices) | reach fixed5 / auto / scaled | paired median tt_hit ratio vs fixed5 (auto / scaled) | read |
|---|---|---|---|
| small (7)  | 1.00 / 1.00 / 1.00 | — (all runs <0.05 s) | **TIE** (differences are proc.time noise; `scaled`==`fixed5` by construction) |
| medium (7) | 1.00 / 1.00 / 1.00 | — (all runs <1 s) | **TIE** (sub-second; tt noise) |
| large (7)  | 0.91 / 0.94 / 0.94 | ~1.07 / ~1.01 | ~tie (reach spread all from project2771) |
| xlarge (4) | 0.95 / 0.80 / 0.80 | ~1.17 / ~1.11 | reach regression all from project4284 |
| **clean large+xlarge** (9, excl. the 2 variance matrices; 45 cells) | tie | wall2hit **1.122 / 1.073**; **rep2hit 1.000 / 1.000** | see granularity note |

Net reach vs fixed5 (of 125 cells): **auto 3 better / 5 worse**; **scaled 4 better / 5 worse**.

### Granularity check (advisor) — the wall2hit gap is NOT a slowdown-to-optimum
`tt_hit`/`wall2hit` is replicate-granular (the callback fires at replicate boundaries). On the clean-tie
matrices the optimum is first hit at **replicate 1 in 33–35/45 cells (median rep2hit = 1, all arms)**, and
the **paired replicates-to-first-hit ratio is exactly 1.000** for both auto and scaled. So the bigger kick
reaches the optimum in the *same replicate* — **no loss of per-replicate search effectiveness**. The
~7–12% `wall2hit` gap is the bigger kick's **costlier ratchet phase within each replicate** (costlier per
replicate, largely post-optimum on these easy matrices), i.e. an accounting cost of the replicate-boundary
timestamp — **not** "slower to the optimum." For the mission metric (time-to-first-reach), auto/scaled are
**neutral** here; the added per-replicate cost buys nothing on the easy-tie bulk (and adds to total
wall-to-stop under production stopping).

### The two matrices that move (everything else ties)
- **project2771 (94 t) — a marginal reach benefit reappears** (auto/scaled reach 0.6 vs fixed5 0.4;
  per-seed 3 better / 1 worse / 1 tie, n=5). Suggestive, not established: the one training matrix where a
  bigger kick helps reach — a genuinely hard **multi-basin** mid-size landscape, echoing the 5432 (482 t)
  benefit and the `intensity-thorough-hardtail` regime, but too noisy at n=5 to call an analogue.
- **project4284 (4062 t × 27 chars) — bigger kick HURTS reach** (auto/scaled 0.2 vs fixed5 0.8; finals
  1–2 steps worse in 4/5 seeds). A near-degenerate, **rep-starved** giant: a bigger per-replicate kick
  buys fewer replicates in the budget → worse. One seed went the other way.

## Verdict
- **Do NOT change the production default `ratchetPerturbMaxMoves` 5→auto (or scaled).** It buys **no broad
  wall-clock benefit** (14/25 matrices tie exactly; on the clean large/xlarge tier it reaches the optimum
  in the **same replicate**, rep2hit ratio 1.000, so no gain and no genuine slowdown), adds per-replicate
  ratchet cost that buys nothing on the bulk, and **regresses reach** on very large low-signal matrices at
  production budgets. A default change needs a broad win; there is none.
- **The over-perturbation-on-small-trees fear is REFUTED** — small/medium tie exactly (auto=20 on a 25-tip
  tree costs nothing measurable), and even on large trees the kick is not *less effective* per replicate
  (rep2hit ratio 1.000) — it is only *costlier* per replicate.
- **The kick's value is landscape-specific** (hard multi-basin: project5432 482 t, project2771 94 t) →
  it belongs on the **opt-in thorough / hard-tail intensity path** (`intensity-thorough-hardtail`,
  `auto-vs-thorough-objective`), **not** the default. `scaled` is marginally less harmful than `auto`
  (gentler) but still net-neutral-to-negative broadly.
- Vindicates the advisor's earlier "dodged a bad default" flag — now on the **wall** axis, not only reach.

## Reproduce
`dev/benchmarks/kick_anytime.R` (harness) + `kick_anytime_hamilton.sh` (SLURM array 1-125) +
`kick_analyze.R` (per-cell union-best, time-to-first-hit, reach, size-stratified). Hamilton:
engine `/nobackup/pjjg18/TreeSearch/kick-lib` (current), source `/nobackup/pjjg18/TreeSearch-kick`,
outputs `/nobackup/pjjg18/TreeSearch/kick_anytime_out/cell_*.csv`. Per-arm summary
`kick_per_arm.csv`.
