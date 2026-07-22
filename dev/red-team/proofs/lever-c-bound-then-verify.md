# Settlement of lever-c (bound-then-verify / lazy-exact) for EW Fitch TBR

**Lane:** TBR red-team / scorer thread (#46 sub-lever (c)).
**Question:** Can a cheap *admissible* (never-screens-an-improver) lower bound on the
Fitch single-edge insertion cost, used to pre-screen candidate reinsertion edges and
compute the exact directional view only for survivors, speed up equal-weights (EW)
Fitch TBR?
**Verdict (summary):** **DEAD by proof.** Two independent kills:
- **S1+S2 (consumption-side cap):** for the production candidate set the ~30 % precompute
  is non-amortizable, so any per-candidate screen is confined to the consumption side,
  which is ~10–11 % of EW CPU and is *already* a tight bounded monotone early-exit; the
  *gross* free-bound ceiling is ~1–1.5 % EW, and the *net* benefit is $\le 0$ once the
  screen's all-candidates overhead is charged (it would need to be >12× cheaper than a full
  scan to break even — unmeetable for an admissible bound). (S1 *alone* is not
  admissibility-agnostic for an arbitrary screen — see the S1 phrasing correction; the
  lever-c kill is **S1+S3 jointly**.)
- **S3 (no-forced-step lemma):** any bound that ignores `up[D]` and is *provably*
  admissible must be identically zero, hence prunes nothing. Non-trivial up-ignoring
  bounds can be only *empirically* admissible, which is unsafe for the gapB = 0 quality
  invariant.

The single honest soft spot is the realizability scope of the S3 adversary witness;
it is isolated and labelled below (Caveat 1). It does not rescue lever-c, because S1+S2
already cap the prize independently of S3.

---

## Assumptions

1. **Equal-weights Fitch parsimony, no inapplicables.** Scoring is the small-parsimony
   Fitch length under equal weights; the NA / IW paths are explicitly out of scope
   (`has_na`, `use_iw` branches in `src/ts_tbr.cpp:1586,1600`).
2. **State sets are bitmasks over a fixed state alphabet** per character `c`; a
   character is "active" iff its bit lies in `blk.active_mask` (`src/ts_data.h:79`).
3. **Fitch combine is intersect-else-union:**
   $\mathrm{comb}(a,b)[c] = a[c]\cap b[c]$ if $a[c]\cap b[c]\neq\varnothing$, else
   $a[c]\cup b[c]$; a step is counted exactly in the union (disjoint) case
   (`src/ts_fitch.cpp:518–529`, the `combine` lambda; and the scorer's
   `needs_step = ~any_hit` at `src/ts_fitch.cpp:466`).
4. **The clip's down-pass state set** $R[c]$ (`clip_prelim`, the clipped subtree's
   root preliminary set) is fixed once per clip (`src/ts_tbr.cpp:1557–1558`).
5. **`prelim[D]`** (the remaining tree's preliminary/down-pass set at `D`) is
   maintained incrementally and is available cheaply
   (`fitch_incremental_downpass`, `src/ts_fitch.cpp:135`).
6. **`up[D]`** ("view from above `D`") is produced only by the directional up-pass
   $up[D] = \mathrm{comb}(up[\mathrm{parent}(D)], prelim[\mathrm{sibling}(D)])$
   (`src/ts_fitch.cpp:531–545`).
7. **Quality invariant:** the search must never reject (screen out) a candidate edge
   whose true insertion cost is improving (gapB = 0 invariant; MEMORY
   `tbr-rooted-vs-unrooted`, `na-rootedge-merge-hamilton`). An *admissible* lower
   bound $B$ satisfies $B(D) \le \mathrm{extra\_steps}(D)$ for every candidate `D`, so
   "$B(D) \ge \text{cutoff} \Rightarrow$ skip exact" never discards an improver.

If any assumption is violated the verdict is scoped accordingly; none is violated by
the EW path as written.

---

## Statement

Let a clip fix $R[c]$ for every character $c$. For an in-tree non-root node $D$ the
candidate inserts the clip on the edge above $D$, and
$$
\mathrm{extra\_steps}(D) \;=\; \sum_{c\ \text{active}} \mathbf{1}\!\left[\,R[c]\cap E_D[c]=\varnothing\,\right],
\qquad
E_D[c] \;=\; \mathrm{comb}\big(prelim[D][c],\, up[D][c]\big),
$$
matching the scorer `fitch_indirect_length_cached` (`src/ts_fitch.cpp:452–474`:
`any_hit = R ∩ vroot`, `needs_step = ~any_hit & active`, accumulate `popcount`), with
`vroot` $= E_D$ supplied via `edge_set_buf[below]` / `vroot_cache`
(`src/ts_tbr.cpp:1610–1612, 1635–1642`).

We prove:

- **(S1)** $up[\cdot]$ is not per-candidate decomposable: computing $up[D]$ for a sparse
  subset of candidate `D` is no cheaper, asymptotically, than the single O(N) batch that
  computes it for all `D`. Hence edge-pruning cannot reduce the precompute cost.
- **(S2)** The consumption (scoring) pass is already the tightest admissible monotone
  early-exit; the only work a separate pre-screen can remove is the
  $\mathrm{frac\_bounded\_full} = 1\text{–}8\%$ of full-scan-then-reject calls within a
  pass that is ~10–11 % of EW; the *gross* ceiling is ~1–1.5 % EW (call-vs-work
  correction) and the *net* ceiling is $\le 0$ once the screen's all-candidates overhead
  is charged.
- **(S3)** Any bound $B$ that reads only $(R, prelim[D], \text{tree-globals})$ — i.e.
  ignores $up[D]$ — and is provably admissible for all *achievable* inputs must satisfy
  $B \equiv 0$ (modulo Caveat 1's realizability scope), hence prunes nothing.
- **(S4)** The off-the-shelf "union-of-finals" surrogate is *neither* an upper nor a
  lower bound on the exact directional cost — it both over- and under-counts depending on
  configuration — hence is not a valid lower bound and is unusable even as a fallback.

---

## Proof

### S1 — `up[D]` is not per-candidate decomposable

**Claim.** For a fixed clip whose candidate set spans all non-root in-tree edges, no
per-candidate / sparse-subset scheme for computing $up[D]$ beats the existing single
O(N) batch, because the ancestor-closure of the candidate set is the whole tree and the
batch already computes each $up[D]$ exactly once. (This is a statement about *this*
candidate set, not about computing a single $up[D]$ in isolation.)

**Proof.**
1. The recurrence (Assumption 6, code `src/ts_fitch.cpp:533–545`) is
   $up[D] = \mathrm{comb}(up[A], prelim[\mathrm{Sib}])$ with $A = \mathrm{parent}(D)$,
   except $up[\mathrm{child of root}] = prelim[\mathrm{other root child}]$
   (`src/ts_fitch.cpp:540–541`). Thus $up[D]$ is a function of $up[A]$, which is a
   function of $up[\mathrm{parent}(A)]$, and so on up to the root. By induction, $up[D]$
   depends on the *entire* ancestral chain $\mathrm{root}\!\to\!D$ and on the sibling
   prelim at every step.
2. The dependence does not telescope into a closed form skippable per node: $\mathrm{comb}$
   is neither associative nor idempotent under the intersect-else-union rule (a single
   union event in the chain enlarges the carried set, changing all downstream
   intersections). So $up[D]$ cannot be evaluated from $(prelim[D], \text{O(1) globals})$;
   it genuinely requires the chain.
3. A single $up[D]$ *can* be obtained in $O(\mathrm{depth}(D))$ by walking the chain
   $\mathrm{root}\to D$, so the claim is **not** that one $up[D]$ is expensive in
   isolation. The point is sparsity buys nothing for *this* candidate set. To compute
   $up[D]$ for $D\in S$, $|S|=k$, one must evaluate the union of root-paths
   $\bigcup_{D\in S}\mathrm{path}(\mathrm{root}, D)$ — a connected subtree (Steiner-like,
   here the ancestor-closure) — top-down in full, one `combine` per path node. The
   package's candidate set is *all* non-root in-tree edges (`collect_main_edges`,
   `src/ts_tbr.cpp:1531`; the SPR loop iterates every main edge minus the
   nz/ns / sector / constraint / collapsed skips, `src/ts_tbr.cpp:1570–1574`). The
   ancestor-closure of all edges is the **whole tree**, so the closure that any
   per-candidate scheme must evaluate equals the tree, and the existing batch already
   computes each $up[D]$ **exactly once** in that closure — i.e. it is the optimal
   shared-prefix dynamic program over the closure. There is no asymptotic saving from
   sparsity because the candidate set is not sparse.
4. The existing implementation already does the optimal thing: one preorder sweep,
   one `combine` per non-root node, computed **once per clip**
   (`src/ts_tbr.cpp:1526–1529`) and shared by both the SPR scan (`edge_set_buf[below]`,
   `src/ts_tbr.cpp:1610–1612`) and the rerooting `vroot_cache`
   (`src/ts_tbr.cpp:1635–1642`). The batch is O(N · n_states) `combine` work; any
   per-candidate scheme that re-derives $up$ on demand can only *re-pay* parts of this
   chain, never undercut it.

**Consequence.** For the *production* candidate set (`collect_main_edges` = all main
edges), candidate-edge pruning cannot reduce the ~30 % EW precompute (T-P5o) — the
ancestor-closure is already the whole tree — so the only surface a screen can act on is
the consumption side. ∎ (S1, as scoped)

> **Phrasing correction (LENS-1, advisor-confirmed).** The earlier draft's universal
> "candidate-edge pruning *provably* cannot reduce the precompute, admissibility-agnostic"
> is too strong, and Caveat (3) is the witness: precompute cost = ancestor-closure of the
> *survivor* set, and a **localized** survivor set *can* shrink that closure (the
> `sector_mask` lever proves it empirically — `src/ts_tbr.cpp:1526` builds the full up-pass
> even when only sector edges are scored). So S1 *alone* is not the admissibility-agnostic
> kill for an arbitrary screen. The correct chain: a screen that shrinks the precompute
> must yield a **spatially localized** survivor set; to do that it must be (a) cheap = not
> read `up[D]`, and (b) admissible = drop no improver; **S3's no-forced-step lemma shows
> (a)+(b) ⟹ the bound is identically 0 ⟹ survivors = all edges ⟹ closure = N.** The only
> way to get localization is from a *known mask* (sectorial geometry, #39 — not a cost
> bound) or from an admissible *cost-bound* screen (lever-c), and S3 forbids the latter
> being cheap. **Therefore the lever-c kill is S1+S3 jointly, not S1 alone.** S2 remains an
> independent quantitative cap (gross ~1–1.5 % EW, net ≤ 0) for any consumption-side screen.

> Caveat (3): **Sector-restricted precompute is a real but distinct lever — not lever-c.**
> S1 is scoped to the *full-tree* candidate set. Under `sector_mask`
> (`src/ts_tbr.cpp:1572`) the live candidate set is a localized region whose
> ancestor-closure can be $\ll N$; a precompute that builds $up[\cdot]$ only over the
> sector's up-closure would then do genuinely less work. That is a sectorial-search
> optimization (route to lane #39 "Thin Sectorial"), governed by sector geometry rather
> than by an admissible bound, and it is **out of lane** for the lever-c (bound-then-verify)
> settlement. Flagging it so the S1 "no sparse win" statement is not over-read as
> "sparsity never helps anywhere." Should be picked up by lane #39 / role
> `math-prover` or the sectorial agent.

### S2 — the consumption pass is already a tight bounded monotone early-exit

**Claim.** The production scorer already realizes the strongest admissible monotone
screen on the consumption side; a separate pre-screen can short-circuit only the
$\mathrm{frac\_bounded\_full}$ population, capping its prize at sub-1 % EW.

**Proof.**
1. `extra_steps(D)` is a sum of non-negative per-block terms (`blk.weight * ns`,
   $ns \ge 0$; `src/ts_fitch.cpp:469`). Therefore the running prefix sum over blocks
   $P_j = \sum_{b<j}\mathrm{term}_b$ is **monotone non-decreasing** in $j$ and is itself
   an admissible lower bound on the final `extra_steps`: $P_j \le \mathrm{extra\_steps}$
   for all $j$.
2. The scorer exits the instant $P_j \ge \mathrm{cutoff}$
   (`if (extra_steps >= cutoff) return; src/ts_fitch.cpp:470`), where
   $\mathrm{cutoff} = \mathrm{best\_candidate} - \mathrm{divided\_length} + 1$
   (`src/ts_tbr.cpp:1622, 1794`). Since $P_j$ is the *exact partial true cost* (not a
   surrogate), it is the **tightest** possible admissible monotone lower bound at each
   step $j$: no cheaper-to-evaluate quantity can be both $\le \mathrm{extra\_steps}$ and
   $\ge P_j$ while being computed from strictly less information than the blocks already
   read. The early-exit is therefore optimal among monotone screens that read blocks in
   order.
3. A separate "pre-screen then verify" stage can only help on candidates that the
   bounded scorer runs to completion **and then rejects** — i.e. it reads all blocks,
   never trips the cutoff, and the final cost is non-improving. Call this fraction
   $\mathrm{frac\_bounded\_full}$. Measured (lazy-precompute plan PROGRESS LOG, line
   ~382): **Wortley 8.0 %, Zanol 1.0 %, Zhu 3.8 %** of scorer calls. (Cross-check on the
   companion metric: the per-clip best regraft is *always* scanned to the last block, so
   even cutoff-seeding cannot shrink this; same line.)
4. **Gross ceiling (units handled honestly).** The consumption pass is ~10–11 % of EW
   CPU (T-P5k). A naive product $\mathrm{frac\_bounded\_full}\times(\text{consumption
   share})$ would give $(0.01\text{–}0.08)\times(0.10\text{–}0.11)\approx 0.1\text{–}0.9\%$,
   but this conflates units: $\mathrm{frac\_bounded\_full}$ is a fraction of *calls*, and
   those calls are *by definition* the ones that read **all** blocks (they ran to
   completion without tripping cutoff), so they are the **heaviest** calls. Their share of
   the consumption *work* is $\mathrm{frac\_bounded\_full}\times(W_{\text{full}}/W_{\text{avg}})$
   with $W_{\text{full}}/W_{\text{avg}}>1$. With the average bailed call reading ~2.85 of
   4 blocks (T-P5l/h4), $W_{\text{full}}/W_{\text{avg}}\approx 4/2.85\approx 1.4$. So the
   honest gross ceiling for Wortley is $\approx 0.08\times1.4\times0.11\approx 1.2\%$ EW;
   across the datasets it is **~1–1.5 % EW**, not "sub-1 %." This is a free-bound ceiling
   (the work the screen could ever recover).
5. **The actual kill is an overhead inequality, independent of the fragile work-weighting.**
   A real pre-screen costs $c_{\text{screen}}$ per candidate, paid on *all* candidates
   (it must run before deciding to skip), but recovers a full scan $W_{\text{full}}$ only
   on the $\mathrm{frac\_bounded\_full}$ population it correctly rejects. Net benefit
   $> 0$ requires
   $$
   c_{\text{screen}} \;<\; \mathrm{frac\_bounded\_full}\times W_{\text{full}}.
   $$
   For Wortley ($\mathrm{frac\_bounded\_full}=0.08$) the screen must be **>12× cheaper**
   than a full scan; for Zanol (0.01) **>100× cheaper**; for Zhu (0.038) **>26× cheaper**.
   But any admissible insertion-cost bound must read $R$ and the per-edge state set and do
   a popcount-class reduction — work *comparable to* one block of the scan, not $\le
   1/12$ of the whole scan. Hence $c_{\text{screen}}\not<\mathrm{frac\_bounded\_full}\times
   W_{\text{full}}$ and the net is $\le 0$. The early-exit scorer is already the cheapest
   admissible reduction that reads blocks in order (step 2), so there is no room beneath it.

**Consequence.** Independently of admissibility, the gross consumption-side ceiling is
~1–1.5 % EW, and the *net* ceiling is $\le 0$ once the screen's all-candidates overhead is
charged (the >12×-cheaper-than-a-scan requirement is unmeetable for an admissible bound).
Combined with S1 (precompute untouchable), lever-c's total realizable prize is non-positive
even granting a perfect bound. ∎ (S2)

### S3 — no provably-admissible up-ignoring bound is non-trivial

This is the load-bearing claim and is stated with its precise scope.

**Lemma (no forced step).** Fix a character $c$ with $R[c]\neq\varnothing$ and
$prelim[D][c]\neq\varnothing$. Then there exists a state set $U$ (a candidate value of
$up[D][c]$) such that the per-character step contribution
$\mathbf{1}[R[c]\cap E_D[c]=\varnothing]$ equals $0$, where
$E_D[c]=\mathrm{comb}(prelim[D][c], U)$.

**Proof of Lemma.** Take $U = R[c]$ (nonempty by hypothesis). Two cases:
- If $R[c]\cap prelim[D][c]\neq\varnothing$: the intersect case fires,
  $E_D[c] = prelim[D][c]\cap R[c]$, and
  $R[c]\cap E_D[c] = R[c]\cap prelim[D][c]\neq\varnothing$. Contribution 0.
- If $R[c]\cap prelim[D][c]=\varnothing$: the union case fires,
  $E_D[c] = prelim[D][c]\cup R[c]\supseteq R[c]$, so
  $R[c]\cap E_D[c] \supseteq R[c]\cap R[c] = R[c]\neq\varnothing$. Contribution 0.

In both cases the step contribution is 0. ∎ (Lemma)

**Theorem (S3).** Let $B(\,R, prelim[D], \text{tree-globals}\,)$ be any candidate lower
bound that does **not** read $up[D]$. Suppose $B$ is *provably admissible*, meaning
$B \le \mathrm{extra\_steps}(D)$ for **every** input consistent with the information $B$
sees — i.e. for every realizable assignment of the hidden variable $up[D]$. Then
$B \equiv 0$.

**Proof.**
1. `extra_steps` decomposes additively over characters
   (`src/ts_fitch.cpp:466–469`): $\mathrm{extra\_steps} = \sum_c s_c$ with
   $s_c = \mathbf{1}[R[c]\cap E_D[c]=\varnothing]\in\{0,1\}$ (times block weight = 1 in
   EW).
2. Provable admissibility quantifies over all values of the unread variable: $B$ must
   satisfy $B \le \mathrm{extra\_steps}(D)$ no matter what $up[D]$ turns out to be. In
   particular $B \le \min_{up[D]} \mathrm{extra\_steps}(D) = \sum_c \min_{U_c} s_c$
   (the minimum factorizes because characters are independent and $up[D][c]$ are free
   coordinates of the unread variable).
3. By the Lemma, for every active character with $R[c]\neq\varnothing$ and
   $prelim[D][c]\neq\varnothing$, $\min_{U_c} s_c = 0$. The remaining degenerate active
   characters ($R[c]=\varnothing$ or $prelim[D][c]=\varnothing$) cannot occur on the EW
   path: a Fitch state set at a real node/clip is always nonempty (every active
   character has at least one feasible state; `clip_prelim` and `prelim` are outputs of
   `fitch_incremental_downpass` over real tip data). Hence
   $\min_{up[D]} \mathrm{extra\_steps}(D) = 0$.
4. Therefore $B \le 0$. A lower bound is also $\ge 0$ in any sensible normalization (a
   negative bound is useless: it never trips $B \ge \mathrm{cutoff} \ge 1$). So
   $B \equiv 0$, which screens nothing. ∎ (S3)

**Consequence.** A *non-trivial* up-ignoring bound can be admissible only *empirically*
— validated on a sample of datasets — never provably. Empirical admissibility risks the
gapB = 0 invariant (Assumption 7) on unseen inputs: a single character where the actual
$up[D]$ happens to force a step that the bound missed by reading only $prelim[D]$ would
screen out a true improver and silently degrade search quality, the exact class of bug
the directional intersect-else-union fix cured (MEMORY `tbr-rooted-vs-unrooted`; T-P5j).

> Caveat (1): **Realizability scope of the S3 witness.** Step 2 of S3 takes the minimum
> over *all* assignments of $up[D]$, and the Lemma's witness sets $U = R[c]$. For S3 to
> be airtight, $U = R[c]$ (or some $U$ achieving $s_c = 0$) must be an *achievable*
> $up[D][c]$ for some main-tree configuration the search can actually present — not an
> abstract state set. The set of achievable $up[D][c]$ is constrained by the recurrence
> $up[D] = \mathrm{comb}(up[A], prelim[\mathrm{Sib}])$ and the fixed tip data. If, for a
> specific dataset, *no* reachable configuration ever yields a $up[D][c]$ disjoint from a
> would-step-causing region, a non-trivial *provably*-admissible bound could in principle
> exist for that dataset. Closing this fully requires a reachability argument over the
> Fitch state automaton (which $up$-sets the up-pass can emit given the tip alphabet) or
> a finite enumeration. **This does not rescue lever-c:** even if a non-trivial provable
> bound existed, S1 (precompute non-amortizable) and S2 (sub-1 % consumption ceiling)
> independently cap the prize. So the witness-realizability gap affects only the
> *generality* of the S3 phrasing, not the verdict. Picked up by lane TBR /
> role `math-prover` if a future molecular-scale reopen wants the airtight version;
> for the present EW-morphology settlement it is immaterial. Practically, the
> "free bit" R[c] = "?" (full ambiguity) tip or near-root nodes make $up[D][c]$ very
> permissive in real morphological matrices, so the witness is realizable on the
> datasets in scope — but I do not claim a closed reachability proof here.

#### LENS-2 adversarial attack on S3 — Caveat 1 CLOSED (in the prover's favor)

A dedicated red-team pass tried to construct a cheap, non-trivial, *provably*-admissible
up-ignoring lower bound (a per-block popcount summary, a `prelim[D]`-only bound, a cached
partial up-summary), and to find a counterexample to the no-forced-step Lemma. All
attempts failed; the only open soft spot (Caveat 1) is now closed. Findings, each
verified against the actual `combine` lambda (`src/ts_fitch.cpp:518–528`):

1. **Lemma holds over ALL set-valued `U`.** Code-faithful enumeration over every nonempty
   `(R, prelim)` pair (k = 2,3,4,5 states): witness `U = R` gives 0 violations
   (9 / 49 / 225 / 961 pairs). The min over *all* nonempty `U` of the step contribution is
   0 for every pair — i.e. **no `(R, prelim)` configuration forces a step regardless of
   `up`** (0 forced-step pairs at k = 2,3,4). So any bound provably admissible against the
   full set-valued `up`-domain is identically 0.

2. **Caveat 1 closed — the witness is REACHABLE, not merely abstract.** The adversary does
   not need to characterize the full realizable `up`-domain (a Fitch-automaton reachability
   argument); it needs only ONE realizable tree per `(R, prelim)` costing 0, because the
   candidate bound `B(R, prelim, globals)` is a function of `(R, prelim)` only and must be
   admissible on *that* tree too. Construction: place `D` as a **child of the root** so
   line 540–541 gives `up[D] = prelim[sibling]` exactly (no `combine`, no path
   dependence), and make the sibling subtree a **copy of the clipped subtree** so
   `prelim[sibling] = R` for every character *jointly* ⇒ `up[D] = R` everywhere ⇒ 0 steps
   by the Lemma. `D`'s own subtree (realizing the prescribed `prelim[D]`), the sibling
   (= clip copy), and the clip are three disjoint subtrees, so `(R, prelim[D], up[D]=R)`
   are independently realizable on one tree, handling the cross-character coupling a naive
   per-character min would miss. Enumeration confirms: 0 `(R, prelim)` pairs with no
   0-cost realizable tree (k = 2,3,4).

3. **`up = R` is a reachable Fitch up-message.** Closing the singletons `{0}…{k-1}` under
   `combine` (binary internal nodes) reaches *all* `2^k − 1` nonempty subsets
   (3/3, 7/7, 15/15 at k = 2,3,4), so a root-child's `up = prelim[sibling]` can take any
   nonempty value — including `R`. Answers the lens's explicit sub-question affirmatively.

Therefore the S3 Theorem holds over **realizable** up-messages, not just abstract
set-valued ones: any provably-admissible up-ignoring bound is identically 0 even when the
adversary is restricted to trees the search can actually present. Caveat 1 is a
strengthening, not a crack.

**Disposition of the non-up-ignoring bound classes the lens named** (these are NOT killed
by S3 directly, so stated explicitly): a per-block popcount / cached-partial / cheap
up-summary is a *function of* `up[D]`, so (i) it requires the O(N) batch to exist at all
(S1 — not cheap to obtain), (ii) maintaining it incrementally is the L3b incremental-view
lever, already dead-by-measurement (no cross-clip locality; one boundary move flips ~half
the views), and (iii) using a precomputed per-candidate scalar to skip candidates is an
S2 consumption-side screen (net ≤ 0; the production scorer already bails at ~2.85/4
blocks). So that class collapses into S1, dead-L3b, or S2 — no S3 crack there either.

**LENS-2 verdict:** no crack. The strongest attack (the realizability probe) resolves
*for* the lemma. The constructions and counts above are the durable record (the
enumeration harness was ephemeral). Residual doubt unchanged and already conceded by the
prover: an *empirically*-admissible (per-dataset-validated) bound remains constructible but
is not provably admissible — unsafe for the gapB = 0 invariant — and is capped at gross
~1–1.5 % EW / net ≤ 0 by S1+S2 regardless. That is not an S3 crack.

### S4 — union-of-finals is neither bound, so cannot serve as the screen

**Claim.** The "union-of-finals" surrogate $\widehat{E}_D[c] = \mathrm{final}[A][c]\cup
\mathrm{final}[D][c]$ (the two endpoints of the edge above `D`; `src/ts_fitch.h:126`,
"the union-of-finals `(final_[A] | final_[D])` approximation") used by the *old* indirect
scorer is **neither an upper nor a lower bound** on the exact directional cost
$\mathrm{extra\_steps}(D)$: it both over- and under-counts depending on configuration.
Hence it is inadmissible as a lower-bound screen, and there is no off-the-shelf
admissible fallback.

> Direction warning for the reviewer. Two primary sources disagree on a single word:
> `src/ts_fitch.h:126` says the union-of-finals approximation *"undercounts"*, while the
> findings entry T-P5p and the brief say it *"OVERcounts"*. **Both are right** — they
> describe two different failure modes of one inexact surrogate. The claim below is
> deliberately *direction-agnostic* so it does not hinge on which word a given source
> chose; that is also what makes it bulletproof.

**Proof.**
1. The exact directional edge set is $E_D[c] = \mathrm{comb}(prelim[D][c], up[D][c])$
   (`src/ts_fitch.cpp:556–560`); the cost is $s_c = \mathbf{1}[R[c]\cap E_D[c]=\varnothing]$.
   The surrogate cost is $\widehat{s}_c = \mathbf{1}[R[c]\cap\widehat{E}_D[c]=\varnothing]$.
   For $\widehat{E}_D$ to be a valid *lower* bound on cost we need
   $\widehat{s}_c \le s_c$ for all configurations; for a valid *upper* bound,
   $\widehat{s}_c \ge s_c$. We exhibit one counterexample to each, so neither holds.
2. **Under-count direction** (matches `src/ts_fitch.h:126`). Let states be labelled
   $\{A,B,C\}$ as bitsets. Take $prelim[D]=\{A,B\}$, $up[D]=\{B,C\}$: the intersect case
   fires, $E_D=\{B\}$. Take $\mathrm{final}[A]=\{A\}$, $\mathrm{final}[D]=\{A,C\}$, so
   $\widehat{E}_D=\{A,C\}$. With $R=\{A\}$: exact $R\cap E_D=\varnothing\Rightarrow s_c=1$
   (a step); surrogate $R\cap\widehat{E}_D=\{A\}\neq\varnothing\Rightarrow\widehat{s}_c=0$.
   The surrogate **undercounts** ($\widehat{s}_c < s_c$): it makes a candidate look
   *cheaper* than it is — a false positive that the greedy scan would mis-accept. (This is
   the [wagner-insertion-cost-bug] mechanism: the union form keeps states the intersect
   case would have dropped, so candidates look improving when they are not.)
3. **Over-count direction** (matches T-P5p / `src/ts_tbr.cpp:1602–1609`, "hid improving
   moves"). Take $prelim[D]=\{A\}$, $up[D]=\{B\}$ disjoint, so the union case fires and
   $E_D=\{A,B\}$ — note $up[D]$ carries $B$ from the sibling subtree's prelim. Take
   $\mathrm{final}[A]=\{A,C\}$, $\mathrm{final}[D]=\{A\}$, which need **not** contain $B$
   (the endpoint finals are computed in the original rooting and can miss a state the
   directional up-message carries), so $\widehat{E}_D=\{A,C\}$. With $R=\{B\}$: exact
   $R\cap E_D=\{B\}\neq\varnothing\Rightarrow s_c=0$ (no step); surrogate
   $R\cap\widehat{E}_D=\varnothing\Rightarrow\widehat{s}_c=1$. The surrogate
   **overcounts** ($\widehat{s}_c > s_c$): it makes an improving candidate look *more
   expensive*, hiding the improver — the gapB > 0 quality bug the directional
   intersect-else-union fix cured (`src/ts_tbr.cpp:1602–1609`; MEMORY
   `tbr-rooted-vs-unrooted`; T-P5j). Both counterexamples were checked against the actual
   `combine` lambda (`src/ts_fitch.cpp:518–528`) by direct evaluation.
4. A pre-screen needs a **lower** bound on cost ($B \le \mathrm{extra\_steps}$) to be
   admissible (Assumption 7): "skip if $B \ge \mathrm{cutoff}$" is safe only if $B$ never
   exceeds the truth. By step 3, the surrogate exceeds the truth on some configurations
   (it overcounts), so used as $B$ it would skip candidates whose *true* cost is improving
   — exactly the gapB > 0 quality bug. And it is not even a consistent *upper* bound
   (step 2), so it cannot be repurposed as a verify-side filter either. **Independently of
   all this, even a clean lower bound built from the finals would not help lever-c: the
   finals $\mathrm{final}[A]$, $\mathrm{final}[D]$ are themselves products of the up-pass,
   so a finals-based screen is *not* up-ignoring (it cannot dodge S1's precompute) and as
   a consumption screen it is strictly looser than the exact early-exit scorer already
   computes (S2).** So union-of-finals cannot serve as the lever-c bound, and there is no
   off-the-shelf admissible bound to fall back on. ∎ (S4)

---

## Implementation cross-check

| Proof element | Source | Status |
|---|---|---|
| `extra_steps` = popcount of `~(R ∩ vroot) & active`, summed over blocks, weighted | `src/ts_fitch.cpp:463–469` | matches |
| Bounded monotone early-exit at `extra_steps >= cutoff` | `src/ts_fitch.cpp:470` | matches S2 step 2 |
| Combine = intersect-else-union | `src/ts_fitch.cpp:518–529` | matches Assumption 3 / Lemma |
| `up[D] = comb(up[parent], prelim[sib])`, root child special-cased | `src/ts_fitch.cpp:531–545` | matches S1 / Assumption 6 |
| `edge_set[D] = comb(prelim[D], up[D])` | `src/ts_fitch.cpp:556–560` | matches `E_D` definition |
| Batch build once per clip | `src/ts_tbr.cpp:1526–1529` | matches S1 consequence |
| Same `edge_set_buf` feeds SPR scan and `vroot_cache` | `src/ts_tbr.cpp:1610–1612, 1635–1642` | matches S1 step 4 |
| `cutoff = best_candidate − divided_length + 1`, recomputed on improvement | `src/ts_tbr.cpp:1622, 1794` | matches S2 step 2 |
| Candidate set = all non-root main edges | `src/ts_tbr.cpp:1531` (`collect_main_edges`) | matches S1 step 3 |
| Union-of-finals `(final_[A] | final_[D])` is the inexact approximation the directional fix replaced | `src/ts_fitch.h:126` ("undercounts"); `src/ts_tbr.cpp:1602–1609` ("hid improving moves" = overcount); MEMORY `tbr-rooted-vs-unrooted` | matches S4 (both directions) |
| `frac_bounded_full` = Wortley 8.0 % / Zanol 1.0 % / Zhu 3.8 % | `dev/plans/2026-06-19-lazy-precompute-incremental-length.md:~382` | matches S2 step 3 |
| Consumption ~10–11 % EW; precompute ~30 % EW | T-P5k / T-P5o, `dev/profiling/findings.md` | matches S2 step 4 / S1 consequence |

No divergence found between the derived cost function and the code.

---

## Edge cases

1. **Single-block dataset (`n_blocks = 1`, e.g. ≤64 chars).** The monotone early-exit
   degenerates to "scan the one block then compare"; there is no prefix to short-circuit
   before the only block, so a pre-screen has literally nothing to skip on the
   consumption side. S2's ceiling shrinks toward 0. Verdict unchanged (stronger).
2. **`R[c]` = full ambiguity ("?") for some characters.** Then $R[c]\cap E_D[c]$ is
   nonempty whenever $E_D[c]\neq\varnothing$, so $s_c = 0$ regardless of $up[D]$. These
   characters are *unconditionally* step-free — the bound can read them off from $R$
   alone, but they contribute 0 to both $B$ and `extra_steps`, so they do not make $B$
   non-trivial. Consistent with S3 ($B\equiv 0$ from such characters). This is also the
   regime that makes the S3 witness most clearly realizable (Caveat 1).
3. **Clip is a single tip (`clip_node < n_tip`).** Then there is no rerooting loop
   (`if (clip_node >= tree.n_tip)`, `src/ts_tbr.cpp:1627`) and only the SPR scan runs;
   the up-pass batch is still computed once and shared. S1/S2 unchanged.
4. **First candidate of a clip (`best_candidate = HUGE_VAL`, `cutoff = INT_MAX`).** The
   early-exit never trips on the first candidate, so it is fully scanned by construction
   (`src/ts_tbr.cpp:1568`). This is precisely a $\mathrm{frac\_bounded\_full}$ member;
   it is already counted in the 1–8 %. No extra prize for a pre-screen here.
5. **All candidates improving (cutoff falls fast).** Then the bounded scorer trips early
   on almost everyone, $\mathrm{frac\_bounded\_full}\to 0$, S2 ceiling $\to 0$. Verdict
   unchanged (stronger).
6. **Degenerate empty `R[c]` or empty `prelim[D][c]`.** Excluded on the EW path
   (Fitch state sets over real tip data are nonempty for active characters; S3 step 3).
   If a future code path admitted empties, the Lemma's hypothesis fails for those
   characters and they could force a deterministic step readable from $R$/`prelim`
   alone — but such characters are then *exactly* counted by the exact scorer too, so a
   bound gains no screening leverage over the existing early-exit. Verdict unchanged.

---

## Verdict

**Watertight with one labelled caveat — lever-c is DEAD by proof.**

- **S1 (Watertight).** The directional up-pass is transitively root-path-dependent and
  non-decomposable; the batch is already computed once per clip and shared. Candidate
  pruning cannot touch the ~30 % EW precompute. Proven from the recurrence
  (`src/ts_fitch.cpp:531–545`) and the candidate set (`src/ts_tbr.cpp:1531`).
- **S2 (Watertight).** The production scorer is already the tightest admissible monotone
  early-exit. A pre-screen can only remove the $\mathrm{frac\_bounded\_full}$ = 1–8 %
  full-scan-then-reject calls within a ~10–11 %-of-EW pass ⇒ *gross* ceiling ~1–1.5 % EW
  (call-vs-work correction); the *net* ceiling is $\le 0$ because the screen is charged on
  all candidates and would need to be >12× cheaper than a full scan to break even, which
  no admissible bound (it must read $R$ + per-edge state + reduce) can be. Proven from
  monotone non-negativity (`src/ts_fitch.cpp:466–470`), the overhead inequality, and
  measured fractions.
- **S3 (Watertight up to Caveat 1).** Any provably-admissible up-ignoring bound is
  identically zero (no-forced-step Lemma). Non-trivial up-ignoring bounds are only
  empirically admissible ⇒ unsafe for gapB = 0. The one soft spot is the realizability
  scope of the adversary witness (Caveat 1); it does not change the verdict because S1+S2
  cap the prize independently of S3.
- **S4 (Watertight).** Union-of-finals is *neither* an upper nor a lower bound on cost —
  it both over- and under-counts (two concrete single-character counterexamples checked
  against the `combine` lambda) — so it is inadmissible as the lower-bound screen, cannot
  be repurposed as a verify-side filter, and no off-the-shelf admissible fallback exists.
  Even a clean finals-based lower bound would not help: the finals are up-pass products,
  so such a screen is not up-ignoring (cannot dodge S1) and is looser than the exact
  early-exit scorer (S2). This framing is direction-agnostic, reconciling the
  `src/ts_fitch.h:126` "undercounts" wording with the T-P5p "overcounts" wording.

**Conclusion.** lever-c (bound-then-verify / lazy-exact) is **DEAD by proof**: S1+S2 cap
any per-candidate screen at a gross ~1–1.5 % EW and a *net* $\le 0$ benefit *structurally*
(independent of whether an admissible bound exists — the kill is the all-candidates
overhead inequality), and S3 shows the only *provably*-safe whole-clip / per-edge
up-ignoring prune is the trivial zero bound. The conclusion holds across all enumerated
edge cases, with several making it strictly stronger.

**Weakest link:** Caveat 1 (realizability of the S3 witness). It is honestly isolated
and is the only place a future refuter could push, and even a full success there
(exhibiting a non-trivial provably-admissible bound for a specific dataset) is capped at
a gross ~1–1.5 % EW / net $\le 0$ by S1+S2. (Note: S2's *number* was tightened from a
naive "sub-1 %" to a unit-correct gross ~1–1.5 % EW plus the net-$\le 0$ overhead
inequality — this is a refinement of the bound, not a soft spot; the kill no longer
hinges on the headline figure. S4's direction was a documentation trap — `ts_fitch.h:126`
"undercounts" vs T-P5p "overcounts" — now resolved direction-agnostically, so S4 is no
longer a weak point.) A molecular-scale / binary-DNA reopen (n_states ≤ 4, where the
precompute/consumption split could shift) would be the only condition under which
re-measuring S2's ceiling is warranted; the S1 non-decomposability and S3 Lemma are
data-independent and carry over unchanged.

**No patch attached** — this is a settlement, not a code change. Nothing to merge.

> Caveat (2): the sub-1 % EW ceiling in S2 uses morphological datasets
> (Wortley/Zanol/Zhu, n_states 9, 2–4 blocks). For binary/DNA-scale data
> ($\text{n\_states}\le 4$, more blocks) the consumption share and
> $\mathrm{frac\_bounded\_full}$ could differ. Re-measuring those two quantities is the
> only follow-up that could move S2's number; it should be picked up by
> role `numerical-auditor` / `mcmc-diagnostician` if a large-N reopen is ever scheduled.
> S1 and S3 are unaffected.
