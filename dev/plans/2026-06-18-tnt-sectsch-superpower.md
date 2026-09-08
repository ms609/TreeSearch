# What gives TNT 1.6 `sectsch` its escape — bare-bones reverse-engineering

Date: 2026-06-18. TNT-ONLY investigation (no TreeSearch `src/`/`R/` edits). Primary
dataset Zanol2014 (74 tips, equal-weights Fitch), confirmed on Wortley2006 & Giles2015.
TNT exe `C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe` (32-bit; `mxram` ≤ 1024).
All scratch + scripts under `dev/benchmarks/tnt_bare/`.

## CONCLUSION (one paragraph)

TNT `sectsch` does **not** escape by any clever per-sector move engine, sector geometry,
global-TBR cadence, drift, or extra search effort. It escapes because it runs sectorial
search over a **retained set of several equally-optimal trees** — `mult` keeps ~10 trees
(default `hold`), and a strict (`noequals`) sectorial sweep over that *topologically diverse*
set reaches the target on every dataset, whereas the identical sweep on a **single** tree
plateaus well above target and *stays there no matter how long it runs* (10× the rounds gives
the identical plateau). The mechanism is **equal-length topological variety**: a strict sector
re-solve only yields a strictly shorter tree when it has access to a *neighbouring* optimal
topology with a different sector arrangement, and a single frozen tree never exposes one. TNT
supplies that variety two ways — (1) the retained diverse optimum set (what it does by
default), and (2) the `equals` option, which lets a *single* tree accept equal-length lateral
sector re-solves and plateau-walk into an improvable configuration. **What TreeSearch should
replicate: drive/maintain the sectorial over a diverse set of equally-optimal trees rather
than polishing one tree** (effort is not a substitute); secondarily, turn ON the already-wired
`sector_accept_equal` lever for the single-tree case.

Two calibrations (detail in "Is the set… SHARING?"): (i) the cross-lane benefit is real but
**modest and tree-level, not sector-level** — there is NO recombination of sectors across trees
(Goloboff: each sector is re-solved against its own tree); the shared *tree buffer* reuses whole
best-found trees as starting points, which beats "10 independent lanes, pick best" at equal
compute (1261 in 7/15 vs 1/15) though medians tie at 1262. (ii) 1261 is reachable by a lone tree
too, just rarely (~0.5%/restart) — the set raises the hit-rate ~7×, it does not unlock a
forbidden score.

This **corrects** `2026-06-17-sectsch-escape-mechanism.md`, whose RANK-1 ("sector geometry +
replacement count") is refuted as the escape source and whose RANK-2 demotion of `equals` was
based on a single, atypical seed and a worse (1275) starting basin. See "Reconciliation".

---

## Method & fixtures

- Fitch matrix: `-`→`?`, equal weights. `WriteTntCharacters` (TreeSearch 2.0.0, lib `.agent-aband`).
- **The `hold` (tree-buffer) value moves the starting basin**: `mult=replic 1` (rseed 1) gives
  1271 under `hold 1000` (10 trees retained) but only **1275** under `hold 1` (1 tree). The TBR
  buffer size is itself a lever. The canonical gap is defined at **T0 = 1271**, so the
  single-tree start is tree #1 of the 1271 set (scores 1271 in both TNT `score;` and `TreeLength`).
- Bare runs read the start tree **fresh** (`proc <tree>`); **no mult/bb/tbr/xmult before sectsch**.
  TNT score authoritative; final tree re-scored with `TreeLength` (min over all saved trees,
  with `tipLabels` resolved from the data path TNT writes) — matches TNT every time (mapping OK).
- Scripts: `harness.R` (reusable runner+parsers), `driver1..5.R`, `confirm.R`,
  `make_single.R`/`setup.R` (fixtures). Live TNT default reports **global TBR every 2
  substitutions in small sectors** — the `tnt_defaults.txt` dump saying "10" is stale.

### The exact BARE-BONES script (deliverable requirement)
```
mxram 1024;
report+;
proc data.tnt;                 ' Fitch matrix
rseed 1;
hold 1000;                     ' working buffer
proc tee.tre;                  ' the FIXED single T0=1271 (no search run yet)
sectsch: noglobal noequals nofuse godrift 9999 ;   ' strip global-TBR, equals, fuse, drift
sectsch = rss ;                ' ...repeated up to 12-300x
score ;
```
**Result: 1271 → 1271, every round, forever (zero escape).** The simplest possible strict
sectsch on the real T0 does NOT reach target. (Default sectsch — global TBR every 2, 3
RAS+TBR/sector — is also 1271 from this tree.) So the vanilla strict move-engine is *not* the
source of the escape.

---

## Knob sweep — three start conditions × acceptance rule (Zanol, 12 rounds; `driver2.R`)

Per-round running-best; final = min `TreeLength` over saved trees (verified equal).

```
A: SINGLE 1271 tree (hold 1000)
  default            1271 ...                                 -> 1271
  noglobal noequals  1271 ...                                 -> 1271
  global 1           1271 ...                                 -> 1271   (max global TBR: no help)
  equals             1265 1263 1263 1263 1262 ... 1261        -> 1261   *** reaches target
  equals global 1    1265 1263 ... 1263                       -> 1263   (global TBR HURTS equals)

B: 10-tree 1271 SET (hold 1000)
  default (noequals) 1263 1261 ...                            -> 1261   *** TNT's actual behaviour
  noglobal noequals  1265 1264 1263 1263 1261 ...             -> 1261   (strict reaches it)
  equals             1261 ...                                 -> 1261

C: in-memory hold-1 mult = 1275 (the prior doc's start)
  noglobal noequals  1272 1269 1264 1262 ...                  -> 1262   (strict stalls 1 ABOVE target)
  equals             1265 1262 1261 ...                       -> 1261
```

## Controls (`driver3.R`, `driver4.R`)

- **Acceptance, not cadence/geometry/fuse.** `nofuse` on the set-strict run is byte-identical
  to default → the set route is **not** tree-fusing. `global 1` never helps and *hurts* when
  combined with `equals`.
- **`equals` accumulates NO buffer diversity** — with `equals` the tree count stays **1** all
  the way down to 1261. So `equals` = single-tree lateral *plateau-walking* (temporal variety),
  mechanistically distinct from the set route (stored variety).
- **Tree count vs diversity.** 10 *identical* copies + strict → 1263 (partial escape: the
  sectorial buffer self-diversifies with tied-length alternatives); 10 *different* trees +
  strict → 1261 (full). A single fixed topology + strict is the *only* fully frozen cell.

## Seed robustness + DIVERSITY-vs-EFFORT (`driver5.R`, Zanol, 30 rounds, seeds 1–6)

```
SINGLE-T0 strict          : min 1264  median 1267  max 1271      {1271,1264,1267,1267,1267,1266}
SINGLE-T0 equals          : min 1261  median 1263  max 1267
SET(10 diverse) strict    : min 1261  median 1261.5 max 1264
SINGLE-T0 strict @300 rnds : min 1264  median 1267  max 1271     (10x effort = IDENTICAL plateau)
```
**Single-tree strict plateaus at median 1267 and 10× more rounds changes nothing** → the set's
advantage is **diversity, not compute**. Note single-strict median **1267 == the task's
"TreeSearch reaches ~1267"**, and set/default == TNT's 1261 — the gap is exactly this lever.

## Cross-dataset confirmation (`confirm.R`, seeds 1–4, 30 rounds, median [min–max])

| dataset (target)   | single strict | single `equals` | SET strict | **SET default = TNT** |
|--------------------|---------------|-----------------|------------|------------------------|
| Zanol2014 (1261)   | 1267 [1264–1271] | 1262.5 [1261–1263] | 1261 [1261–1262] | **1261 [1261–1262]** |
| Wortley2006 (479)  | 485 [482–485]    | 480.5 [479–485]    | 480 [479–482]    | **479.5 [479–480]**  |
| Giles2015 (670)    | 671.5 [671–672]  | 670 [670–671]      | 670 [670–670]    | **670 [670–670]**    |

Universal ordering: `single-strict` (worst, above target) ≫ `single-equals` ≈ `set-strict` ≈
`set-default` (= target). The retained-set route reaches target on all three; the single-tree
strict route never does.

## Is the set "10 independent lanes, pick best", or genuine SHARING? (`driver6–C.R`)

The N trees are processed with `tree` = "all trees" (default): **10 different trees → 10
incomparable sector sets**, NOT one tree's sectors re-solved 10× (that is the flat 300-round
effort control). The honest answer took several careful tests and overturned earlier wording:

- **NO sector-level recombination.** Goloboff 1999's RSS re-solves each sector against *its own
  tree's* scaffold ("best among the R+r replications AND the present resolution… place it in the
  whole tree") and `nofuse` is a no-op. Structure from tree A is **never** spliced into tree B.
- **1261 IS reachable by a lone lane — just rare.** 200 independent single-tree restarts: 1/200
  reached ≤1261 (~0.5%/restart); the set is *not* reaching an otherwise-impossible score
  (`driverA.R` TEST1). (An earlier "0/50 → only the set can" was undersampling.)
- **Identical copies in a shared buffer ≈ separate restarts.** 10 copies of one tree, shared
  buffer, vs 10 separate restarts of it: ties (1262–1263) — no coupling detectable when the
  starting tree is fixed (`driverA.R` TEST2; uninformative for trees that *can* reach 1261).
- **But "together" beats "apart" at EQUAL compute.** Paired, 15 reps, both = 10 trees × 30
  rounds: SET reached 1261 in **7/15**; 10 independent lanes (same trees, take best) **1/15**;
  paired set<indep 8, tie 6, set>indep 1 (sign test p≈0.04) — `driverB.R`. The 1/15 for
  independent matches the ~0.5%/lane rate; the set's 7/15 does not. So the set is **NOT**
  equivalent to "10 lanes, pick best."
- **The advantage is NOT sector size.** Forcing a single tree onto sectors up to size 70
  (near-global, n=74) still gives median 1267, 0/6 reaching 1261 — identical to default size 37
  (`driverC.R`). So the shared-size-counter idea is refuted too.
- **By elimination, the channel is the shared TREE buffer (tree-level, not sector-level).** Not
  recombination, not size, not effort (300-round flat), not starting-diversity-reaching-the-
  unreachable. The 10 lanes draw from and write to one common tree pool, so **whole improved
  trees discovered by any lane become starting points reused by the ongoing search** — a
  beam/population effect. Corroboration: the set has a 1261 in its buffer by round ~8
  (`driver8.R`), yet *no* isolated lane reaches 1261 even in 30 rounds (0/36 in `driverC`), so
  that 1261 cannot be a single independent-lane trajectory. The exact retention/reuse rule is a
  TNT-internal (closed source) not line-traced here; the buffer stays diverse (lengths 1261–1271
  coexist at the plateau), so it is not simple best-culling.

**Bottom line on "collective":** there is cross-lane information flow, but at the granularity of
**whole trees** (the shared buffer pools and reuses the best trees found by any lane), NOT at the
granularity of sectors. Magnitude is modest (medians tie at 1262; the effect is a ~7× higher
chance of hitting the 1261 optimum), but real and reproducible.

## Sector-size schedule — TNT vs TreeSearch (`driverC.R`, `driverD.R`)

TNT default sector size = **min(n/2, 45)** → 37 = n/2 for Zanol (n=74). Within ONE `sectsch=rss`
invocation it ramps: do M = (T·100)/((100−selfact)·S) ≈ 3 selections at size S, then S → S×1.75
(`increase 75`), toward ~n. BUT the size **resets to n/2 at the start of every invocation**
(settings dump is byte-identical after each round: "run 3 sectors of 37 nodes"). So across a
looped search the operative size is just **n/2**.

And the ramp is inert here: SET with `minsize 37 maxsize 37`, with `… increase 0` (escalation
OFF), and default(escalating) are **identical** (3/6 reach 1261, med 1261.5). Forcing a SINGLE
tree onto fixed sizes 37→70 is also flat (`driverC`: all 0/6, med 1267). **Sector size is not a
lever for this escape; the buffer is.**

TreeSearch (read-only, `src/ts_sector.cpp:1044-1051`, `R/MaximizeParsimony.R`) does NOT use n/2:
it collects EXISTING internal clades with size in a band [`sectorMinSize`,`sectorMaxSize`] =
[6,50] default (80 thorough, 100 large) and picks ~`2·n_tip/avg_size` ≈ 5 of them (random or
conflict-weighted). So its sectors are a wide size *distribution* of existing clades (max can
EXCEED n/2), vs TNT's single ~n/2 walked-up clade. (The "n/2 capped at 45 / for n<~88" rule is
TNT's `selectem`, not TreeSearch's current code.) Since size is not the lever, this difference
is not what drives the gap.

## MINIMAL sufficient configuration
- **From the fixed canonical T0=1271 (clean isolation):** the decisive factor is **the number
  of distinct equally-optimal trees the sectorial operates over**. One tree (strict) → median
  1267, never target, even at 10× rounds. The 10-tree diverse set (strict, `noequals` =
  TNT default) → 1261. Same start, same effort budget; only the retained-set diversity differs.
- **From a forced single tree:** **`sectsch: equals;`** is the single sufficient knob (reaches
  1261 on Zanol seed 1; median 1262–1263 over seeds). Necessary too — every strict single-tree
  config plateaus above target.
- *(Aside, not a clean isolation:* in the end-to-end pipeline `hold 1` → 1262 vs `hold 1000` →
  1261, but `hold 1` also shifts the `mult` output to the worse **1275** basin, so that
  comparison conflates buffer size with start quality — the clean evidence is the fixed-T0
  single-vs-set contrast above.)

## Reconciliation with `2026-06-17-sectsch-escape-mechanism.md`
- Its "**-13 strict bulk escape with `noglobal`**" was measured from the **1275** (`hold 1`)
  in-memory tree, which is *not* sector-optimal, so strict moves trivially exist — but even
  there strict stalls at **1262**, one step above target (block C). From the real T0=1271,
  strict does *nothing* on seed 1 and plateaus at 1267 on average. The "bulk" was a worse-basin
  artifact, not the T0 escape.
- It ranked sector **geometry/replacement-count #1** and **`equals` a minor #2 bridge**. The
  controlled fixed-start experiment inverts this: geometry/cadence/effort do not move the single
  T0; variety (set or `equals`) does. Its `noglobal`-barely-changes-it observation is consistent
  with mine (`global` is not the lever) — but it concluded the *strict sector replacements* were
  the escape, which the diversity-vs-effort control refutes.

## What the parent TreeSearch project should replicate (priority order)
1. **PRIMARY — sectorial over a retained DIVERSE set of equally-optimal trees**, not a single
   polished tree. (Prior doc notes TS picks a single existing clade on one tree.) Effort cannot
   substitute. Keep the set of optimal trees `mult`/RAS produces and let sector improvements
   propagate across it.
2. **SECONDARY — flip ON `sector_accept_equal`** (`src/ts_driven.h:94`, default `false`,
   already plumbed → `SectorParams::accept_equal`, comment "Goloboff 2014 plateau lever"). This
   is the single-tree substitute (median ~1262–1263; sometimes hits target). Cheap to test.
3. **NOT the lever:** global-TBR cadence (flat alone, harmful with equals), sector
   geometry/sub-clade collapse, recursion, more rounds/effort. Do not invest there for THIS gap.
