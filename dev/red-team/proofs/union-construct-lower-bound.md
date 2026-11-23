# Is Goloboff (1996)'s "union construct" screen a sound lower bound on TBR/SPR insertion cost?

**Lane:** math-prover · **Target:** the Mission-B Goloboff union-construct screen (un-built) ·
**Harness:** `dev/red-team/union-construct-gate.R` (pure R, no engine needed).

**Verdict (one line):** **SOUND lower bound — but ONLY with TRUE MPR finals `M(A)` for the
ancestor term; with the engine's COLLAPSED `uppass_node` finals `fin(A)` it is UNSOUND
(over-rejects genuinely-improving moves).** So the screen is buildable-correct, and its
buildability condition is: the ancestor's most-parsimonious set must be the exact directional
edge set `E(A)=combine(prelim(A),up(A))=M(A)` (which the engine already computes per-clip), NOT
the cheap `tree.final_` collapse.

> Status: the formal exchange-lemma write-up was interrupted by an account session limit
> (resumes 4pm Europe/London). This file records the DECISIVE empirical gate result + the
> analytic structure; the full step-by-step exchange proof can be completed on resume. The
> empirical gate is exhaustive at small N and is on its own sufficient to gate the build.

## Claim

Fitch (unordered, symmetric). Clipped subtree with basal state set `r` reinserted in the main
subtree. For internal node `D` with ancestor `A`, define per character:
- `union_set(D)` = ∪ of leaf state sets over all descendants of `D` (a plain down-pass union).
- **union construct** `UC(D) = union_set(D) ∪ MP(A)`.
- added length of joining `r` to a set `S` = `w·[r ∩ S = ∅]` (indirect-length rule).

CLAIM (coverage/lower bound): for every descendant edge `D'` in `subtree(D)` (incl. `D`),
`M(D') ⊆ UC(D)`, where `M(D')=combine(P(D'),Up(D'))` is the exact MPR edge set (proven exact in
the sibling doc `union-of-finals-bound-proof.md` S1). By antitone-ness of `[r ∩ · = ∅]`, this
gives `added_UC(D) ≤ added_exact(D')` for all `r` — so if the UC estimate ≥ cutoff, every
location in `subtree(D)` is safely rejected.

## Result (empirical gate — `dev/red-team/union-construct-gate.R`)

Reference validated first: `E(D)=combine(P,Up)` equals the brute-force per-edge added length by
physically attaching a tip (0 mismatches / 5,600 edges), so `E`/`M` is the right exact object.

| regime | `vio_TRUE` (UC with `M(A)`) | `vio_COLL` (UC with `fin(A)`) |
|---|---|---|
| exhaustive nTip=4 k=3 | **0** | 1,656 |
| exhaustive nTip=5 k=3 | **0** | 10,344 |
| exhaustive nTip=5 k=4 | **0** | 25,632 |
| exhaustive nTip=6 k=4 | **0** | 111,984 |
| exhaustive nTip=7 k=3 | **0** | 31,704 |
| polymorphic nTip=6 k=4 (pPoly .4, pMiss .1) | **0** | 6,859 |
| polymorphic nTip=7 k=4 | **0** | 6,114 |
| polymorphic nTip=8 k=5 | **0** | 5,666 |

Plus a direct added-length cross-check (sample clip masks `r`, confirm `added_UC_true(D) ≤
added_exact(D')` for every descendant `D'`): **0 over-rejections / 180,850 checks, worst excess 0.**

- `vio_TRUE = 0` everywhere (exhaustive to 7 tips; polymorphic + missing to 8 tips / 5 states)
  ⇒ **UC with exact MPR finals `M(A)` is a sound coverage/lower bound**, including at the base
  edge (`D'=D`), the `A=root` case, and under polymorphism/missing data.
- `vio_COLL > 0` massively ⇒ **UC with the engine's `uppass_node` collapse `fin(A)` is NOT a
  lower bound.** Minimal counterexample (4 tips, one character, pattern `1,1,0,0`): the collapsed
  final drops an equally-parsimonious ancestor state, so `UC_coll` misses a state of some
  `M(D')` and the screen **over-rejects a genuinely free (zero-added-length) insertion** — i.e.
  it would discard an improving/neutral move, exactly the failure mode that sank the
  union-of-finals screen (`union-of-finals-bound-proof.md`), in the same direction.

## Why (analytic structure)

`union_set(D) ⊇ P(D)` always (union ⊇ the intersect-else-union down-pass). So the only states of
`M(D')` at risk of escaping `UC(D)` are the *up-side* states (`Up(D')\P(D')`, present only in the
union case `P∩Up=∅`). The exchange argument shows every such up-side state, for any `D'` in
`subtree(D)`, lies in `M(A)` (it is reconstructible at `A` in some MPR, because it arrives at `D'`
from *above*, i.e. through `A`). This is precisely the containment the gate confirms with 0
violations. It **fails for `fin(A)`** because `uppass_node` keeps only `fin(A)⊆P(A)` (never unions
a parent state in — sibling doc Step 2), dropping exactly the equally-parsimonious ancestor states
the up-side needs.

## Buildability condition (the headline for Mission B)

A correct union-construct screen must use `MP(A) = M(A) = E(A)` (the exact directional edge set,
already computed per-clip by `compute_insertion_edge_sets`), **never** `tree.final_[A]`. It also
needs new `union_set(D)` machinery (a per-node descendant-union, maintainable incrementally à la
Goloboff). Soundness is therefore achievable but is a *soundness obligation on the finals source*,
not automatic — being "up-aware" is not sufficient (the collapsed finals are also up-derived yet
unsound).

## Relation to the mission decision

This proof establishes the screen *can* be built soundly. Whether it *should* is a separate
question answered NEGATIVELY by the Gate-#3 head-to-head (`dev/profiling/mission-b-goloboff-gates.md`):
TS already reaches the optimum at ≤1.13× TNT's wall on the most kernel-heavy large dataset, so a
candidate-count-reducing screen has no time-to-optimum gap to close on this corpus.
