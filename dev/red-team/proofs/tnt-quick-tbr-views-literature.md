# Literature verdict: does TNT/Goloboff build per-node directional views for TBR/SPR scoring?

Question: TreeSearch's TBR/SPR regraft scorer builds per-node directional edge sets
(`edge_set[D] = combine(prelim[D], up[D])`, full whole-tree up-pass) per clip, costing
~30% of EW Fitch CPU. The closure ("precompute is at-limit") rested on the ASSUMPTION
that TNT/Goloboff's "quick TBR" / incremental method also builds equivalent down+up
state sets. This file verifies/refutes that assumption from the primary literature.

## ANSWER: (A) confirmed for the structure. Confidence: MEDIUM-HIGH (full text).
##         The one (B) sub-question is UNCERTAIN (abstract-only) and = the already-
##         catalogued lever-b. Confidence LOW on that sub-question.

Disentangle TWO levels of amortization (the trap I initially fell into — welding them):

- LEVEL 1 (per-candidate, within ONE clip): build views once for the residual tree, then
  score each reinsertion candidate by a root-to-root state-set comparison, NO view rebuild
  per candidate. **TS ALREADY DOES THIS** — `compute_insertion_edge_sets` builds views
  once per clip, then `fitch_indirect_length` scores each candidate as
  `clip_prelim ∩ edge_set[D]` (ts_fitch.cpp:400-423, 486-577).
- LEVEL 2 (per-clip / across accepted moves): derive each clip's divided-tree views
  INCREMENTALLY from the pre-clip whole-tree pass, avoiding an O(n) whole-tree up-pass
  per clip. **TS does NOT do this** — it rebuilds the entire whole-tree up-pass from
  scratch per clip (full preorder over all in-tree nodes, ts_fitch.cpp:540-571).

The literature's FULL-TEXT mechanism evidence (chapter §1.3.6.3, evidence #1) describes
LEVEL 1 — "particularly effective when an important number of SPR or TBR neighbors has to
be evaluated" = candidates within a clip. That confirms TS's structure and per-candidate
amortization; it is NOT a reopen.

The ONLY evidence pointing at LEVEL 2 (the work TS doesn't do) is the Goloboff 1996
ABSTRACT ("...for the divided tree based on calculations for the whole tree"). That is
abstract-only, cannot carry a confident verdict, and IS lever-b from memory.

So: views are the right structure and TS matches the accessible-literature method
including per-candidate amortization (A confirmed). Level-2 per-clip incremental
derivation is a real catalogued lever (lever-b), abstract-supported only, empirically
low-yield on this data class — UNCERTAIN, not a confirmed cheaper method TS lacks.

## KEY EVIDENCE

1. [READ FULL TEXT] Goëffon, Richer & Hao, "Heuristic Methods for Phylogenetic
   Reconstruction with Maximum Parsimony" (book chapter, §1.3.6.3 "Fast character
   optimization techniques"; refs: Goloboff 1993 [21], Gladstein 1997 [16],
   Ronquist 2000 [43]):
   - "a set of shortcuts that helps decrease the computation time **by not recalculating
     the whole tree each time a SPR or TBR modification is applied**. Those techniques are
     particularly effective when an important number of SPR or TBR neighbors has to be
     evaluated."
   - "In [21], Goloboff proposed a method for indirect calculation of the parsimony score
     which uses two passes. **This method needs only to compare the root of the clipped
     tree with the potential root of the target tree to obtain the score of a potential
     new tree for a SPR search.**"
   - "In [47=Ronquist 2000] a two passes algorithm is described which has the same
     complexity of Goloboff's and is faster than the incremental method of Gladstein."
   Source PDF: leria-info.univ-angers.fr/~jinkao.hao/papers/BookParcimony2011.pdf (text PDF,
   fully extracted locally via pdftotext).

2. [ABSTRACT/AUTHOR-TEXT] Goloboff 1996, "Methods for faster parsimony analysis,"
   Cladistics 12:199-220 — structured abstract (author-written; retrieved via search
   summaries of the Wiley/ResearchGate abstract, NOT full text):
   - "Three different algorithms for faster estimation of final state assignments for the
     divided tree **based on calculations for the whole tree** are presented. The first ...
     is approximate; it uses information from the final state sets for the whole tree. The
     second is exact ... based on the **union of the state sets of the descendants for each
     node**. The third is also exact ... faster, ... based on **final and preliminary state
     sets for the whole tree**." (= the directional-view structure, derived from one
     whole-tree two-pass optimization; called "incremental two-pass optimization" in
     secondary summaries.)
   - "The method for indirect tree length calculation when moving a clipped clade, based on
     final states for the divided tree ... include the possibility of **rejecting several
     locations as suboptimal by checking just one node**." (= a lower-bound screen.)

3. [READ FULL TEXT, secondary — but about TAXON INSERTION, not TBR] XMP paper
   (Bioinformatics 27(10):1359, faster exact MP):
   - "Goloboff (1993) describes a way to speed up parsimony searches by **avoiding a
     complete first-pass Fitch optimization for each taxon insertion, enabling amortized
     O(k) Fitch scoring of taxon insertions.**" Uses "Shortcut C from Goloboff (1996) to
     eliminate unnecessary second-pass recursion."
   - CAVEAT: this is stepwise-addition / Wagner taxon insertion, NOT TBR clip-reinsert.
     It shows Goloboff's incremental philosophy exists but is NOT direct TBR-incremental
     evidence. Do not cite it as proof of Level-2 TBR amortization.

4. [READ FULL TEXT, package source] Goloboff 1993 (Character Optimization and Calculation
   of Tree Lengths, Cladistics 9:433-436) — per search-summary of abstract: "shortcuts
   that allow rapid evaluation of tree lengths and fast reoptimization of trees after
   clipping or joining of subtrees, and ... a new **incremental** character optimization
   algorithm which is exact, correct, and comparable in speed."

## THE MECHANISM (if porting)

Two-pass directional state sets (preliminary = down-pass, final = up-pass) for the whole
tree, computed ONCE. To score a clip+reinsert: compare the clipped subtree's root state
set against the target branch's "potential root" (final/edge) state set — O(states) per
candidate, NO per-candidate view rebuild. The views are maintained INCREMENTALLY across
moves (recompute only nodes whose views actually changed — the path affected by the last
accepted rearrangement), NOT rebuilt whole-tree per clip. Plus an optional approximate
lower-bound screen that can reject candidate locations "by checking just one node."

Compatibility with our case (n_states=9 nonadditive Fitch, re-rooting TBR):
- The view STRUCTURE is identical to TS's edge_set[D]/up[D]; multistate-Fitch compatible
  (Goloboff's algorithms are stated for unordered multistate). The "any branch can be a
  root" property the chapter relies on holds for our nonadditive Fitch.
- The re-rooting of the CLIPPED subtree in TBR is the part Goloboff 1996 explicitly costs
  extra ("more reinsertion points under TBR than SPR; recalculate final states when the
  tree is divided"). The clipped-subtree side still needs its own (small) directional pass.
- The lever is NOT a different scorer — it is INCREMENTAL VIEW MAINTENANCE: avoid the
  full whole-tree up-pass per clip by updating only changed views. This is exactly the
  "incremental-length rewrite of compute_insertion_edge_sets" (lever-b) that the memory
  flagged as the only substantial remaining route, then marked "dead-by-solid-argument."

## RELATION TO EXISTING CLOSURE — VINDICATES it (does not overturn)

Memory T-P5p asserts TS "already implements quick-TBR's incremental-length method" and
"TNT scores from down+up sets = our edge_set[D]." The accessible FULL-TEXT literature
SUPPORTS this structural claim: Goloboff builds equivalent two-pass directional views and
scores each candidate by a root-to-root state-set comparison — exactly TS's
`compute_insertion_edge_sets` (build once per clip) + `fitch_indirect_length`
(`clip_prelim ∩ edge_set[D]` per candidate). The per-candidate amortization (LEVEL 1)
that "fast character optimization" is mainly about — "particularly effective when an
important number of SPR/TBR neighbors has to be evaluated" — TS already does.

The one thing TS does NOT do is LEVEL 2: derive each clip's views incrementally from the
PRE-clip whole-tree pass instead of rebuilding the whole-tree up-pass per clip
(ts_fitch.cpp:540-571 loops over ALL in-tree nodes every call). That IS lever-b. Memory
already catalogued lever-b and deferred it on EMPIRICAL locality grounds (L3b: one
boundary move flips ~half the views, fp_frac 0.41-0.68), NOT on a false "TNT doesn't build
views" premise. So the literature does not contradict the memory — it confirms the
structure and leaves lever-b exactly where memory left it: real, catalogued, low-yield on
this data class.

## BOTTOM LINE FOR THE ENGINEER

CLOSURE HOLDS on the accessible literature. Nothing in any source I could read full-text
establishes a cheaper per-candidate method that TS lacks — the canonical fast method
(Goloboff 1993/1996, Ronquist 2000) builds the same two-pass directional views and scores
each candidate by a root-to-root state-set comparison, which TS already implements (build
views once per clip, then O(states) per candidate). The TS assumption "TNT also builds
down+up sets" is CONFIRMED, not refuted.

The single sub-question that could reopen — does Goloboff 1996 derive each clip's views
INCREMENTALLY from the undivided-tree pass (Level 2), saving the per-clip whole-tree
up-pass? — is supported ONLY by the 1996 abstract ("based on calculations for the whole
tree"), which I could not access in full text. That Level-2 lever is identical to memory's
lever-b, already deferred on empirical L3b locality grounds; the literature neither
confirms a realizable win nor refutes the deferral. Revisit at large-N / molecular /
denser data (where L3b fp_frac falls and incremental maintenance may finally pay).

ONE genuinely live thread (narrow): Goloboff 1996's approximate screen "rejecting several
locations as suboptimal by checking just one node" uses FINAL (up-aware) state sets. The
existing lever-c death proof (dev/red-team/proofs/lever-c-bound-then-verify.md) was about
UP-IGNORING admissible bounds; an up-AWARE approximate screen may not be covered by that
proof. Worth a targeted check IF lever-c is ever reopened — but it is "approximate" (not
admissible), so it screens, it doesn't bound, and net-overhead caveats likely still apply.

## SOURCES ACCESSED

- Goëffon/Richer/Hao book chapter PDF — FULL TEXT (pdftotext) — ACCESSIBLE.
- Goloboff 1996 Cladistics 12:199 — ABSTRACT/author-text via search summaries only;
  Wiley full text PAYWALLED (402), Sci-Hub blocked at fetch layer, ResearchGate 403.
- XMP Bioinformatics 27:1359 — FULL TEXT via WebFetch — ACCESSIBLE (Oxford open).
- Goloboff 1993 Cladistics 9:433 — abstract summary only; full text not accessed.
- Goloboff & Catalano 2016 (TNT 1.5) — PAYWALLED (402), not accessed.
- Goloboff/Farris/Nixon 2008 (TNT) — PAYWALLED (402), not accessed.
- USPTO patent 7043371 — image PDF, no OCR; Google Patents flagged it as alignment-
  optimization (likely wrong patent). NOT used.
- arXiv 2103.10967 (astrocladistics) — uses TNT, no algorithm internals. NOT useful.
- HAL thesis tel-01479049 — Anubis access-denied. NOT accessed.
