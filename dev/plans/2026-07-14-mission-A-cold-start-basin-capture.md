# Mission A — Cold-start basin capture (objective (i): always find the best score)

**Framing (2026-07-14).** One of two distinct missions split off the per-move investigation. This is the
REACH mission. Canonical hard case: **project5432** (482t, 189 chars, ~44% missing, EW Fitch). TS's
production search reliably lands at ~1947 (best-ever across every config); TNT reaches **1943** ~1/3 of
runs and **1944** reliably. Gap = 3 steps / ~0.15% on this one matrix.

## The diagnosis is DONE and it is actionable (read [[project5432-basin-structure]] in memory first)
- **1943 is a VALID, STABLE TS optimum.** Seeded at 1943, TS's own TBR HOLDS it (0 improving moves);
  MaximizeParsimony retains it. NOT a completeness / representation / rooting problem.
- **TS's TBR REBUILDS the whole 1943 tree from any start within ~8–12 TBR moves of it** (restore238 job
  17866684: break the deepest 238-tip backbone split, 22/24 truly-adjacent starts recover 1943 with all
  64 skeleton splits; capture radius ~8–12 generic).
- **The wall is PURELY cold-start basin capture:** cold RAS+Wagner lands ~109 TBR moves from 1943 and
  never routes within the ~8–12 capture radius. ~2150+ cold reps → 0 hits. The mechanism gap is
  CONFIRMED BY OBSERVATION (head-to-head job 17863597/8: TS 0/20 vs TNT 4/12 → 1943).
- **=> The mission: get a cold start (or a search trajectory) to land within ~8–12 TBR of the best basin
  on large trees.** You do NOT need an operator that hits the optimum exactly — normal TS descent
  finishes from within capture radius.

## What is REFUTED — do NOT re-derive (each has a Hamilton job + memory entry)
- Distant-basin / diversity-generation pivot — DEAD (1943 intermingled, CID 0.33 ≈ TS-internal spread).
- **BSS** (backbone-skeleton / freeze-a-seed-and-rearrange) — STRUCTURALLY refuted (ceiling job 17863572:
  K_CEILING=323; 1945's block internals are wrong pervasively → needs near-per-tip re-opt = the full
  search, not a coarse freeze-and-rearrange operator).
- **Fuse / recombination** — DEAD: TS `tree_fuse` is pure clade-exchange; 85/479 of 1943's clades are
  absent from ALL near-optimal TS trees (both a 10008-tree diverse pool and the genuine best-34) → cannot
  be assembled, only GENERATED. Iterative/transient autoconstraint bootstraps the WRONG (1947-basin)
  backbone. Rigid consensus autoconstraint compat-killed.
- **Random-TBR perturbation strength** — universal landscape property: reach-from-near-miss-basin = 0 at
  every kick size on BOTH reachable (4138) and unreachable (5432) floors. Not the discriminator.
- Clip-ordering (largefirst/antitip/invweight) — clean null. Suboptimal buffer, addseq variety,
  sector-size growth, tbrMaxHits 3→20, rasStarts>3, big-pool fuse, deep-narrow (outerCycles+unlimited
  resets) — all refuted or null (see [[intensity-thorough-hardtail]], [[tnt-feature-gap-audit]]).

## The ONE lever that moved the floor — the live starting point
**`ratchetPerturbMaxMoves` (character-reweighting kick), 5 → ~40–auto, moved 1949 → 1945.** This is a
DIFFERENT mechanism from random-TBR (reweighting distorts the landscape → re-optimizes the backbone
skeleton sectorial can't touch); it coherently assembles ~49/64 skeleton splits, then plateaus at 1945
(2 steps from 1943). The fixed default 5 is an undersized kick on large trees (engine auto = n_tip/8 ≈ 60);
this is a CONFIRMED recipe defect.
- **First task (cheap, owed, no build): close its corpus-wide anytime/wall-clock test** (memory: "wall
  benefit on solved matrices UNTESTED"). If it helps broadly → a deployable default fix (needs small-tree
  non-regression + user OK to change a default). Reach benefit was 5432-only (n=1); wall benefit is open.

## Genuinely-open angles (the mission's search space — none yet refuted)
1. **A cold-start CONSTRUCTOR that generates starts near the target basin.** The fan-out's original goal
   stands even though BSS (one instantiation) is dead. Candidates: constraint-guided / char-informed
   addition sequence; progressive assembly; a start-generator that lands within capture radius. The bar:
   does it land TS within ~8–12 TBR of the best basin?
2. **Push the reweighting kick as a basin-HOPPING schedule** (the one thing that worked) — strength,
   schedule, multi-kick — beyond the single-knob ratchet default. Plateaus at 1945; can a schedule cross
   the last ~2 steps?
3. **Smarter breadth**, not naive: 2150+ uniform cold reps = 0 hits, so uniform restarts are refuted;
   the open question is whether STRUCTURED / diversified starts raise the per-rep basin-capture rate.

## Discipline (the campaign's hard-won rules — obey them)
- **Existence-before-build.** BSS was saved from a 250-line dead build by a structural precondition test.
  Buildless probes first (`ts_tbr_diagnostics` runs the real main-loop TBR from ANY start tree). Build
  only with a confirmed, ground-truthed target.
- **No local heavy compute** — Hamilton SLURM for anything CPU-greedy; builds + ~30s tests local.
- Set a clear terminal bar: if no mechanism gets TS within capture radius after a bounded effort,
  characterize honestly and stop — 0.15% / 1 matrix / TNT-1/3-reliable is a legitimate wall.

## Precursor housekeeping (do this first)
Rewrite `dev/plans/2026-07-11-project5432-hard-tail-characterization.md` — it still asserts the FALSIFIED
"clade-generation / TS can never generate the split" story (restore238 killed it). Correct it to the
cold-start-basin-capture framing. (Currently untracked in the working tree; do NOT commit as-is.)

## Durable artifacts (Hamilton `/nobackup/pjjg18/`)
`floors/project5432_tnt_floor.tre` (TS-1943); `tntwatch/` (TNT watch logs, watchlog_*.txt); `reeval/`
(restore238.R, basin_5432.R, scaffold.R, ceiling.R, geom_full.R, clade_check.R, hold_check.R). Working
TNT: `TERM=dumb /nobackup/pjjg18/TreeSearch/tnt/TNT-bin/tnt < script`; floor recipe
`sect:slack 20; xmult=replic 120 hits 20; bbreak=tbr`. Matrix `.../lgsweep/matrices/project5432.nex` = m01.
Full diagnostic chain: memory [[project5432-basin-structure]] (long node — read it).
