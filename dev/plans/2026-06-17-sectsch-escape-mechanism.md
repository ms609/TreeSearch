# How TNT's `sectsch=rss` escapes a single-sector-optimal tree — mechanism, from primary sources

Date: 2026-06-17. READ-ONLY analysis (no `src/` edits). Independent of, and CORRECTING,
`2026-06-17-tnt-algorithm-audit.md` (whose RANK-1 hypothesis D1 is refuted below).

Sources: Goloboff 1999 (Cladistics 15:415-428), full text at
`C:/Users/pjjg18/Zotero/storage/TETHI9A5/.zotero-ft-cache`; TNT defaults
`dev/benchmarks/tnt_defaults.txt`; TNT help `dev/benchmarks/tnt_help.txt`.
New empirical probes this session: `dev/benchmarks/diag_tnt_noglobal_probe.R`,
`diag_tnt_seq_accum.R`; existing oracle `dev/benchmarks/d1_confirm.out`.

## Headline

TNT's escape from an identical TNT-`mult` T0 is large and FAST: one RSS round drops
Zanol 1275->1264, a second 1264->1262, then plateaus (`diag_tnt_seq_accum.R`). From the
SAME T0, TreeSearch's TBR finds 0 AND its sectorial finds 0 (`sectorial_shared.csv`:
ts_tbr==ts_sect==start on every gap dataset).

By process of elimination against the DEFAULTS, with two new TNT-side probes:

- (d) recursion — OFF by default (`tnt_defaults.txt` "Recursion ... disabled"); `recurse2`
  == `default` on all 4 datasets. NULL.
- (b) global-TBR cadence — REFUTED as the primary lever. `sectsch: noglobal;` (kills ALL
  global TBR) BARELY changes the escape: Zanol -13->-13, Zhu -8->-8, Wortley -5->-4,
  Giles -4->-2. `global 1` (max cadence) does NOT help and slightly HURTS. The -8..-13
  bulk happens with NO global swapping. Goloboff's "globally suboptimal under TBR" framing
  is real but is the CLEAN-UP, not the barrier-crosser, for these n=37-76 EW cases.
- (a) accept-equal laterals — REAL but SMALL. `sectsch: equals;` adds ~1 step on EVERY
  dataset and REACHES the sectsch target on 2/4 (Zanol -13->-14=1261; Giles -4->-5=670;
  Zhu -8->-9; Wortley +1). Default is `noequals`, so laterals are NOT how the bulk escapes,
  but they are the final bridge to the endpoint.
- THE BULK (-8..-13) is sequential strict-on-the-reduced-score sector REPLACEMENTS over
  TNT's large, overlapping, sub-clade-collapsed sectors — with NO global TBR and NO equal
  moves needed.

## (a) Acceptance criterion — EXACT answer

Goloboff p.418-419 step 3: "Choose the best among the R+r replications AND the present
resolution for the sector and place it in the whole tree." TNT default (`tnt_defaults.txt`
line 20) = `noequals` = "Not accepting equally good subtrees". So "best ... and the present
resolution" is STRICT: a re-solve replaces the present arrangement only if its REDUCED score
is strictly lower; an equal-length-but-different re-solve is NOT taken by default. With
`equals` ON, equal re-solves are taken (the help: "[no]equals accept equally good subtrees").

Crucial invariant (user-verified, re-confirmed in the audit's trace): for EW-Fitch the
from-above HTU makes `reduced = full - const` with const = rest-of-tree standalone Fitch
length, INVARIANT to how the sector re-roots. Therefore a strict reduced-score improvement
is identically a strict FULL-tree improvement (audit §3: 0 gate-bites). TNT's per-move accept
and TreeSearch's `new_score < result.best_score` (`ts_sector.cpp:1140`) are THE SAME GATE for
EW. The gate is NOT the gap.

## (b) Global-TBR cadence — role, and why it is NOT the escape

Goloboff step 4 + p.419: "A round of global swapping of the entire tree is made every 5 to
10 replacements, as that number makes it likely that (through clade substitution) the tree
will have become globally suboptimal under TBR." Mechanism as described: accepted sector
substitutions can leave the tree in a state where a cross-region TBR move now improves it;
periodic global TBR harvests that. TNT default (`tnt_defaults.txt` line 18) = global TBR
every 10 substitutions. TreeSearch runs ONE global TBR at the END of all picks
(`ts_sector.cpp:1199-1210`), looped by `rssRounds`.

EMPIRICAL REFUTATION as the primary lever (`diag_tnt_noglobal_probe.R`): with `noglobal`,
TNT still escapes -13/-8 on Zanol/Zhu. So for these EW cases the barrier is crossed by the
sector replacements themselves; global TBR is a secondary clean-up. (Goloboff's framing is
about Zilla, n=500 — a much larger, cleaner composite-optima case where the cadence matters
more.)

## (c) Large sector selection — the actual mechanism (Appendix 1 decoded)

`selectem()` (Appendix 1, OCR-decoded; `5`=`=`,`,`=`<`,`.5`=`>=`,`2`=`-`):
```c
min_sz = (sector_sz * 80) / 100;                       // sector is 80-100% of cap
for (nod = rand() % root; clad_sz[nod] < min_sz;)      // random node, walk UP
    nod = anc[nod];                                    //   until clade >= 80% of cap
items = marknodes(list, nod, 0, marker);               // mark clade
for (a=items; a--;) if (list[a] < ntax) marker[list[a]] = 2;   // all TIPS -> terminals
if (clad_sz[nod] >= sector_sz) {                       // clade too big: COLLAPSE sub-clades
    for (...) if (!marker[x] && ...>=min_sz) {
        marknodes(inlist, x, 1, marker);
        marker[x] = 2;                                 // sub-clade x becomes ONE composite terminal
        if ((cur_sz -= clad_sz[x]-1) <= sector_sz) break;
    }
}
marker[nod] = 2;                                       // basal node = HTU terminal
```
`marknodes` traverses left/right in RANDOM order (`side = 1 & rand()`).

Two consequences:
1. The sector is LARGE (cap = `min(n/2,45)`; here ~n/2 for n<=90) and obtained by walking UP
   from a random node, so it spans many small clades / much of the backbone.
2. When the clade exceeds the cap, whole sub-clades are COLLAPSED into single composite
   terminals (their first-pass state set). The reduced RAS+TBR then reshuffles ~n/2 UNITS
   that are themselves entire sub-clades, against each other and the rest-of-tree HTU. A move
   that relocates one composite unit = transplanting a whole multi-taxon clade across the
   backbone in the full tree — a large-radius move.

TreeSearch instead SELECTS a single EXISTING clade in a size band [6,50] (`ts_sector.h:21-22`,
`rss_search:1037-1044`) and TBR-polishes it (default `ras_starts=1`), which is REDUNDANT with
the global TBR that already produced T0. It never collapses sub-clades, never spans the
backbone, and samples only `n_picks = 2*n_tip/avg_size` ~ 4-5 sectors/round (`ts_sector.cpp:1025`)
vs TNT's ~20-25.

## (d) Recursion — not part of the default escape

`tnt_defaults.txt`: "Recursion (user-defined searches) disabled". `recurse2` == `default`
empirically. NULL.

## Why the audit's RANK-1 hypothesis (D1, frozen HTU attachment) is REFUTED

The audit argued TNT escapes by letting the HTU float (jointly re-resolve + re-root the
sector). The session left an oracle in code: `TS_FREE_HTU_PROBE` (`ts_sector.cpp:867-888`)
runs, per sector, an UNCONSTRAINED reduced search (HTU = ordinary floating (S+1)th leaf, 20x
RAS+TBR) and prints `<<D1-CONFIRM` iff `free < anchored`. By the `reduced=full-const`
invariance, free<anchored would prove a shorter FULL tree the anchored search cannot reach.

RESULT (`dev/benchmarks/d1_confirm.out`, 26 sectors, all 4 datasets): **ZERO `<<D1-CONFIRM`.**
Every `FREEHTU` line has `free >= anchored` — equal on small sectors (the float adds nothing),
strictly WORSE on large sectors (cold from-scratch RAS is weaker than the warm anchored polish:
e.g. Zhu sect76 free=663 vs anchored=631). The 392 `<<FLOAT-IMPROVES` lines are a RED HERRING:
they fire when a SCRAMBLED random restart (s>0, orig=1592 etc.) beats its own bad start by
floating the HTU; that floated score never beats the warm s=0 (=T0) arrangement. Freeing the
HTU forecloses NOTHING. D1 is dead. (This also matches the user's direct finding: 20 RAS+TBR
floating-HTU re-solves per sector of T0 beat T0 on no sector.)

## RANKING of candidates a-d (by evidence)

1. **(c) sector geometry + replacement count** — large, overlapping, sub-clade-collapsed,
   ~n/2-tip sectors, ~20-25 per round, walked UP from random nodes. This is what makes the
   sequence of strictly-full-improving moves AVAILABLE that TreeSearch's small-clade-band,
   ras1, ~5-pick sectorial never proposes. PRIMARY. (Consistent with the band-shape probe
   finding the lever "real"; the prior "marginal-to-noise" verdict was under END-TO-END with
   aggressive ratchet that SUBSUMES it — from a frozen T0 with ratchet OFF it is the lever.)
2. **(a) accept-equal laterals** — `equals` adds the final ~1 step and reaches target on 2/4.
   SECONDARY bridge. NOT default, NOT the bulk.
3. **(b) global-TBR cadence** — `noglobal` barely dents the escape on these n<=76 EW cases.
   TERTIARY clean-up (matters more at Zilla scale).
4. **(d) recursion** — null. Not in default.
5. **(D1) HTU float** — REFUTED by the in-tree oracle.

## MINIMAL modification to TreeSearch's strict-descent per-sector sectorial

To replicate TNT's escape from a frozen T0, in priority order:

1. **Change sector SELECTION to TNT's `selectem` (PRIMARY).** Replace the "existing-clade-in-
   [min,max]" pick (`rss_search:1037-1044, 1078-1088`) with: pick a random node, walk UP via
   `parent[]` until `subtree_size >= 0.8 * cap` (cap = `min(n/2,45)`); take that clade; if it
   exceeds cap, COLLAPSE descendant sub-clades into composite terminals until size <= cap
   (this is exactly `build_reduced_dataset_collapsed`, already in the tree — wire it to the
   walk-up selection, not just to oversized clades). This makes large, backbone-spanning,
   sub-clade-collapsed sectors whose RAS+TBR re-solve proposes large-radius full-tree moves.
2. **Raise sectors-per-round to TNT's count.** Set `rss_picks_per_round` ~ `T*100/(50*S)`
   (TNT's `selfact` law, `tnt_defaults.txt` line 13) ~ 20-25, not `2T/S` ~ 5
   (`ts_sector.cpp:1025`). More sequential replacements per round = more chances to chain.
3. **Set `rasStarts=3` (+3 on score disagreement).** No preset sets it; stays 1
   (`ts_sector.h:29`), which is pure TBR-polish (redundant with global TBR). TNT does 3+3
   (`tnt_defaults.txt` lines 15-16). The re-solve must be a RAS REBUILD to propose new
   topologies, not a re-polish.
4. **Turn `accept_equal` ON for the sector re-solve (the bridge step).** Already wired
   (`sectorAcceptEqual` -> `params.accept_equal` -> `search_sector`, `ts_rcpp.cpp:1399`).
   Worth ~1 step and reaches target on 2/4. NB: keep the strict FULL-tree gate (item 8b) for
   safety under NA; for EW it is equivalent anyway.
5. **Global-TBR cadence is LAST.** Optionally move the single end-of-round TBR
   (`ts_sector.cpp:1199`) to fire every ~10 accepted replacements INSIDE the pick loop.
   Lowest expected yield on these datasets (noglobal barely changed TNT), but cheap and
   matches Goloboff for larger data.

Expected: items 1-3 (geometry + count + RAS rebuild) recover the bulk (-8..-13 from T0); item
4 adds the last step to the sectsch target; item 5 is scale insurance. The HTU-float rewrite
the audit proposed is NOT needed and would not help (refuted).
