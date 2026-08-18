# Is `fitch_indirect_length` (union of endpoint finals) a sound lower bound on Fitch insertion cost?

**Lane:** math-prover · **Target:** `src/ts_fitch.cpp` union-of-finals insertion approximation
**Verdict (one line):** **Implementation bug (documented invariant is false).** For the finals the
engine actually deploys, the union-of-finals approximation is a sound **upper** bound (it
**over-counts**), the exact opposite of the in-code claim at `src/ts_fitch.cpp:406-410`. It is
therefore **not** a sound lower bound and must not be used as a pruning/rejection screen.

---

## Theorem (informal title)

Consider the approximate insertion scorer `fitch_indirect_length` (`src/ts_fitch.cpp:397-439`),
which charges an extra Fitch step for inserting a leaf with state set $r$ on the tree edge
$(A,D)$ whenever $r$ is disjoint from the **union of the two endpoints' final state sets**,
$U = \mathrm{final}(A)\cup\mathrm{final}(D)$. I prove three facts, per character:

1. **(S1 — exact scorer is correct.)** The exact directional scorer
   `fitch_indirect_length_cached` over the edge set $E(D)=\mathrm{combine}(\mathrm{prelim}(D),
   \mathrm{up}(D))$ computes exactly the true added Fitch length: $E(D)=M(D)$, the
   most-parsimonious-reconstruction (MPR) state set of the virtual root on edge $(A,D)$.

2. **(S2 — deployed union OVER-counts.)** With the finals the engine actually stores in
   `tree.final_` (the `uppass_node` rule, `src/ts_fitch.cpp:41-70`),
   $$U=\mathrm{final}(A)\cup\mathrm{final}(D)\ \subseteq\ E(D)\quad\text{for every character,}$$
   so the union approximation's added length is **$\ge$** the exact added length. It is a sound
   **upper** bound. This is the reverse of the containment asserted in the code comment.

3. **(S3 — union of *true* finals UNDER-counts.)** If instead the finals were the true MPR sets
   $M(\cdot)$, then $M(A)\cup M(D)\supseteq E(D)$, so that union under-counts and **is** a sound
   lower bound. This is what the code comment's argument silently assumes — but the engine does
   not compute MPR finals.

I also show **(S4)** that neither variant is exact on binary / unambiguous data (an explicit
all-singleton-tip, two-state, fully-resolved counterexample), refuting the hoped-for
"exact-on-binary" property.

Everything below is *per character*; the reported length is the weighted sum
$\sum_c w_c\,[\,r_c\cap S_c=\varnothing\,]$ over characters $c$, and each per-character containment
propagates monotonically to that sum (§ Proof, step 8).

---

## Assumptions

1. **Data / model.** Unordered (Fitch) parsimony, one character at a time; states are encoded as
   non-empty bitmasks over $\{0,\dots,k-1\}$; weights $w_c\ge 0$ are per-character multipliers
   (`blk.weight`, `blk.upweight_mask`). All results are stated per character and summed.
2. **Tree.** A rooted binary tree on the *clipped* remnant, with a **degree-2 root** (the code
   treats `root = n_tip` as a 2-child vertex; `src/ts_fitch.cpp:585-586`). Every non-root node $D$
   has a parent $A$ and a unique sibling $\mathrm{Sib}$.
3. **Fitch `combine`.** For state sets $a,b$:
   $\mathrm{combine}(a,b)=a\cap b$ if $a\cap b\neq\varnothing$, else $a\cup b$
   (`src/ts_fitch.cpp:563-574`; `intersect-else-union`).
4. **Down-pass (prelim).** $\mathrm{prelim}(v)=$ observed set at a tip, else
   $\mathrm{combine}(\mathrm{prelim}(c_1),\mathrm{prelim}(c_2))$ over the two children.
   Write $P(v):=\mathrm{prelim}(v)$.
5. **Directional up-pass.** For non-root $D$ with parent $A$, sibling $\mathrm{Sib}$:
   $\mathrm{up}(D)=\mathrm{combine}(\mathrm{up}(A),P(\mathrm{Sib}))$, with the special case
   $\mathrm{up}(D)=P(\mathrm{Sib})$ when $A$ is the root (`src/ts_fitch.cpp:576-590`).
   Write $U\!p(D):=\mathrm{up}(D)$.
6. **Edge set.** For non-root $D$: $E(D)=\mathrm{combine}(P(D),U\!p(D))$
   (`src/ts_fitch.cpp:601-609`). This is the exact scorer's virtual-root set.
7. **Deployed finals.** $\mathrm{final}(\cdot)$ is the array `tree.final_` populated by
   `uppass_node` (`src/ts_fitch.cpp:41-70`): $\mathrm{final}(\text{root})=P(\text{root})$; and for
   non-root $D$ with parent $A$,
   $$\mathrm{final}(D)=\begin{cases}\mathrm{final}(A)\cap P(D) & \text{if } \mathrm{final}(A)\cap P(D)\neq\varnothing,\\ P(D) & \text{otherwise.}\end{cases}$$
   Write $\mathrm{fin}(v):=\mathrm{final}(v)$. **Crucially, this rule never unions in a state from
   the parent.**
8. **Freshness.** `tree.prelim` and `tree.final_` are the fresh, mutually consistent full
   down-pass / up-pass of the *current* clipped tree (no staleness). The over-counting below is a
   property of the *formula*, present even under this ideal assumption; staleness can only add
   further, unsigned error (§ Edge cases).
9. **Insertion semantics.** Attaching a leaf with state set $r$ on edge $(A,D)$ adds a step iff $r$
   is disjoint from the reconstructible set at the attachment point. Exact scorer uses $E(D)$;
   union scorer uses $U=\mathrm{fin}(A)\cup\mathrm{fin}(D)$.

---

## Statement (formal)

For a fixed character and a non-root node $D$ with parent $A$:

- **(S1)** $E(D)=M(D)$, where $M(D)=\{x:\exists$ an MPR of the clipped tree assigning state $x$ to
  the virtual root of edge $(A,D)\}$. Hence for any insert set $r$,
  $[\,r\cap E(D)=\varnothing\,]$ equals the true indicator of an added step.
- **(S2)** $\mathrm{fin}(A)\cup\mathrm{fin}(D)\subseteq E(D)$. Hence
  $[\,r\cap(\mathrm{fin}(A)\cup\mathrm{fin}(D))=\varnothing\,]\ \ge\ [\,r\cap E(D)=\varnothing\,]$,
  i.e. $\text{added}_{\text{union}}\ge\text{added}_{\text{exact}}$ (over-count / upper bound).
- **(S3)** $M(A)\cup M(D)\supseteq E(D)$, hence with true MPR finals
  $\text{added}_{\text{union-true}}\le\text{added}_{\text{exact}}$ (under-count / lower bound).
- **(S4)** There exist fully-resolved binary trees with **all tips unambiguous singletons** on
  which (S2) is a *strict* over-count and (S3) is a *strict* under-count; so neither union is exact
  on binary data.

---

## Proof

Throughout, fix a single character with state universe of size $k$; all sets are subsets of
$\{0,\dots,k-1\}$; all quantities $P,U\!p,E,\mathrm{fin},M$ are that character's slices.

### Step 0 — Fitch two-value structure and the MPR set (used for S1, S2, S3)

For a single unordered character define the standard subtree/rest cost functions
(Fitch 1971; Hartigan 1973; Felsenstein 2004, *Inferring Phylogenies*, ch. 2):
$$g_v(x)=\min(\text{steps below } v \mid v=x),\qquad h_v(x)=\min(\text{steps above } v \mid v=x).$$
The Fitch down-/up-pass establishes the **two-value property**
$$g_v(x)=\gamma_v+[\,x\notin P(v)\,],\qquad h_v(x)=\eta_v+[\,x\notin U\!p(v)\,],$$
where $\gamma_v,\eta_v$ are constants (the minimum below/above), $[\cdot]$ is the 0/1 indicator,
and $P(v),U\!p(v)$ are exactly the down-/up-pass sets of Assumptions 4–5. The total length with
$v=x$ is $g_v(x)+h_v(x)$; the tree length is $L^\*=\min_x(g_v+h_v)$ and
$$M(v)=\{x:g_v(x)+h_v(x)=L^\*\}.$$
Because $[x\notin P]+[x\notin U\!p]\in\{0,1,2\}$ and equals $0$ iff $x\in P\cap U\!p$,
$$L^\*=\gamma_v+\eta_v+[\,P(v)\cap U\!p(v)=\varnothing\,],\qquad
M(v)=\begin{cases}P(v)\cap U\!p(v) & P\cap U\!p\neq\varnothing,\\ P(v)\cup U\!p(v) & P\cap U\!p=\varnothing.\end{cases}$$
That is exactly $M(v)=\mathrm{combine}(P(v),U\!p(v))$.

### Step 1 — (S1): the exact scorer equals the MPR set

By Step 0, $E(D)=\mathrm{combine}(P(D),U\!p(D))=M(D)$. Placing the virtual root on edge $(A,D)$ is
the standard rerooting whose MPR set is $M(D)$; inserting a leaf $r$ there adds a step iff no MPR
of the augmented tree keeps length, iff $r\cap M(D)=\varnothing$ (Swofford & Maddison 1987, *Math.
Biosci.* 87:199–229). So `fitch_indirect_length_cached` over $E(D)$ is exact. $\square$

> Empirical backstop: the parallel exactness gate reports $E(D)=$ full-rescore truth on
> **6066/6066** candidates (fresh finals); my own random/ground-truth harness found **0** mismatches
> over $\sim\!3\times10^5$ insertions.

### Step 2 — deployed finals are contained in prelim: $\mathrm{fin}(v)\subseteq P(v)$

Immediate from Assumption 7: $\mathrm{fin}(\text{root})=P(\text{root})$, and for non-root $D$,
$\mathrm{fin}(D)$ is either $\mathrm{fin}(A)\cap P(D)\subseteq P(D)$ or $P(D)$. The `uppass_node`
kernel (`src/ts_fitch.cpp:59-63`) builds `new_val` only from `prelim(node) & final(anc)` or
`prelim(node)` — it never contributes a state from `final(anc)` that is absent from
`prelim(node)`. So $\mathrm{fin}(v)\subseteq P(v)$, always non-empty. $\square$

### Step 3 — Lemma (extend-to-child): $\mathrm{fin}(A)\subseteq M(D)$, given $\mathrm{fin}(A)\subseteq M(A)$

Take $x\in\mathrm{fin}(A)\subseteq M(A)$: there is an MPR $R$ of the clipped tree (length $L^\*$)
with $A=x$. Let $D_R$ be $R$'s state at $D$; since $R$ is optimal, the subtree below $D$ is
optimized for $D=D_R$, contributing $g_D(D_R)$, and the edge $(A,D)$ costs $[x\neq D_R]$. Build
$R''$: keep $R$ everywhere except reassign the subtree below $D$ optimally for $D=x$ (cost
$g_D(x)$) and set edge $(A,D)$ to $[x\neq x]=0$. Then
$$\text{len}(R'')-L^\*=g_D(x)-[x\neq D_R]-g_D(D_R)
=\big([x\notin P(D)]-[D_R\notin P(D)]\big)-[x\neq D_R].$$
If $x=D_R$ the bracket is $0$. If $x\neq D_R$ then $[x\neq D_R]=1$ and
$[x\notin P(D)]-[D_R\notin P(D)]\le 1-0=1$, so the difference is $\le 0$. But $R''$ is a valid
reconstruction, so $\text{len}(R'')\ge L^\*$, forcing $\text{len}(R'')=L^\*$. Hence $x\in M(D)$.
$\square$

### Step 4 — Lemma (case-ii): $\mathrm{fin}(A)\cap P(D)=\varnothing\Rightarrow P(D)\cap U\!p(D)=\varnothing$, given $\mathrm{fin}(A)\subseteq M(A)$

Contrapositive: assume $P(D)\cap U\!p(D)\neq\varnothing$, pick $z\in P(D)\cap U\!p(D)$; show
$\mathrm{fin}(A)\cap P(D)\neq\varnothing$. Two cases on $A$'s down-pass:

- **$A$ is an intersection node** ($P(D)\cap P(\mathrm{Sib})\neq\varnothing$): then
  $P(A)=P(D)\cap P(\mathrm{Sib})\subseteq P(D)$, so by Step 2
  $\mathrm{fin}(A)\subseteq P(A)\subseteq P(D)$; since $\mathrm{fin}(A)\neq\varnothing$,
  $\mathrm{fin}(A)\cap P(D)=\mathrm{fin}(A)\neq\varnothing$. ✔

- **$A$ is a union node** ($P(D)\cap P(\mathrm{Sib})=\varnothing$, so $P(A)=P(D)\sqcup P(\mathrm{Sib})$):
  - If $A=\text{root}$: $U\!p(D)=P(\mathrm{Sib})$, so $z\in P(\mathrm{Sib})$ while $z\in P(D)$
    contradicts disjointness — this configuration cannot arise (vacuous). ✔
  - If $A\neq\text{root}$: $U\!p(D)=\mathrm{combine}(U\!p(A),P(\mathrm{Sib}))$ and $z\notin P(\mathrm{Sib})$.
    If $U\!p(A)\cap P(\mathrm{Sib})\neq\varnothing$ then $U\!p(D)=U\!p(A)\cap P(\mathrm{Sib})\subseteq
    P(\mathrm{Sib})\ni z$, contradiction; hence $U\!p(A)\cap P(\mathrm{Sib})=\varnothing$ and
    $U\!p(D)=U\!p(A)\cup P(\mathrm{Sib})$, giving $z\in U\!p(A)$. Then $z\in P(D)\subseteq P(A)$ and
    $z\in U\!p(A)$, so $P(A)\cap U\!p(A)\neq\varnothing$ and by Step 0 $M(A)=P(A)\cap U\!p(A)$.
    Using $U\!p(A)\cap P(\mathrm{Sib})=\varnothing$,
    $$P(A)\cap U\!p(A)=(P(D)\sqcup P(\mathrm{Sib}))\cap U\!p(A)=P(D)\cap U\!p(A)\subseteq P(D).$$
    By hypothesis $\mathrm{fin}(A)\subseteq M(A)=P(A)\cap U\!p(A)\subseteq P(D)$; since
    $\mathrm{fin}(A)\neq\varnothing$, $\mathrm{fin}(A)\cap P(D)\neq\varnothing$. ✔

So $\mathrm{fin}(A)\cap P(D)=\varnothing\Rightarrow P(D)\cap U\!p(D)=\varnothing$. And by Step 0,
$P(D)\cap U\!p(D)=\varnothing\Rightarrow M(D)=P(D)\cup U\!p(D)\supseteq P(D)$. $\square$

### Step 5 — $\mathrm{fin}(v)\subseteq M(v)$ for every $v$

Induction down the tree. Base: $\mathrm{fin}(\text{root})=P(\text{root})=M(\text{root})$ (at the
degree-2 root $U\!p$ is empty; $M(\text{root})=P(\text{root})$). Step: assume
$\mathrm{fin}(A)\subseteq M(A)$ for parent $A$; prove $\mathrm{fin}(D)\subseteq M(D)$.
- If $\mathrm{fin}(A)\cap P(D)\neq\varnothing$: $\mathrm{fin}(D)=\mathrm{fin}(A)\cap P(D)\subseteq
  \mathrm{fin}(A)\subseteq M(D)$ by Step 3.
- If $\mathrm{fin}(A)\cap P(D)=\varnothing$: $\mathrm{fin}(D)=P(D)\subseteq M(D)$ by Step 4. $\square$

### Step 6 — (S2): $\mathrm{fin}(A)\cup\mathrm{fin}(D)\subseteq E(D)$

$\mathrm{fin}(A)\subseteq M(D)$ (Step 3) and $\mathrm{fin}(D)\subseteq M(D)$ (Step 5), so
$U=\mathrm{fin}(A)\cup\mathrm{fin}(D)\subseteq M(D)=E(D)$ (Step 1). $\square$

> Exhaustive/random backstop: $U\subseteq E$ held on **958,016 / 958,016** edges (rooted binary
> trees up to 8 tips, up to 5 states); the two sub-lemmas $\mathrm{fin}(A)\subseteq E(D)$ and
> $\mathrm{fin}(D)\subseteq E(D)$, and the case-ii lemma, each had **0** violations.

### Step 7 — (S3): with true finals, union is a superset of $E(D)$

If $\mathrm{final}(\cdot)$ were the true MPR sets $M(\cdot)$, then $\mathrm{final}(D)=M(D)=E(D)$
(Step 1), so $M(A)\cup M(D)\supseteq M(D)=E(D)$. This is the containment the code comment relies on;
it is correct **for MPR finals**. (Indeed $M(D)$ alone already equals $E(D)$, so the union with
$M(A)$ is a merely looser — still valid — lower bound.) $\square$

### Step 8 — from per-character containment to the reported length

For insert set $r$ and any candidate set $S$, `fitch_indirect_length*` adds
$w_c\cdot[\,r_c\cap S_c=\varnothing\,]$ for character $c$ (`src/ts_fitch.cpp:432-435`, `504-...`).
The indicator $[\,r\cap S=\varnothing\,]$ is **antitone** in $S$: shrinking $S$ can only turn a $0$
into a $1$, never the reverse. By Step 6, $U\subseteq E(D)$ per character, so
$[\,r\cap U=\varnothing\,]\ge[\,r\cap E(D)=\varnothing\,]$; summing with non-negative weights,
$$\text{added}_{\text{union}}=\sum_c w_c[\,r_c\cap U_c=\varnothing\,]\ \ge\
\sum_c w_c[\,r_c\cap E_c=\varnothing\,]=\text{added}_{\text{exact}}.$$
Symmetrically (Step 7), the true-finals union gives $\le$. This proves (S2) and (S3). $\blacksquare$

---

## Implementation cross-check

| Proof object | Source (`cpp-search` @ `d8c5998`) | Notes |
|---|---|---|
| Union scorer $[\,r\cap(\mathrm{fin}(A)\cup\mathrm{fin}(D))=\varnothing\,]$ | `src/ts_fitch.cpp:397-439` | reads `tree.final_[a_base]`, `tree.final_[d_base]` (`:428-429`) via `any_hit_reduce3` (checks $r\cap(f_A\cup f_D)\neq\varnothing$); step iff `~any_hit` (`:432`). |
| **False invariant** | `src/ts_fitch.cpp:406-410` | Comment: "union … is a **superset** of the true directional Fitch edge set … hence it **UNDER-counts** (never over-counts)." **Both clauses are wrong for `tree.final_`**: Step 6 gives $U\subseteq E$ (a *subset*), so it over-counts. |
| Deployed finals rule (Assumption 7) | `src/ts_fitch.cpp:41-70` (`uppass_node`) | `final(D)=final(A)∩prelim(D)` if non-empty else `prelim(D)` (`:59-63`); **never unions parent state in** ⇒ $\mathrm{fin}\subseteq\mathrm{prelim}$ (Step 2). This is *not* the MPR up-pass. |
| `combine` = intersect-else-union | `src/ts_fitch.cpp:563-574` | matches Assumption 3 / Step 0. |
| Directional up-pass $U\!p(D)$ | `src/ts_fitch.cpp:576-590` | `up[D]=combine(up[parent],prelim[sib])`; root ⇒ `up=prelim[sib]` (`:585-586`). |
| Edge set $E(D)=\mathrm{combine}(P(D),U\!p(D))$ | `src/ts_fitch.cpp:601-609` | consumed by exact scorer. |
| Exact scorer (S1) | `src/ts_fitch.cpp:475-...` (`fitch_indirect_length_cached`) over `vroot=edge_set[D]` from `compute_insertion_edge_sets` (`:521-618`). | proven exact = $M(D)$. |
| `any_hit_reduce` / `any_hit_reduce3` semantics | `src/ts_simd.h` (`any_hit_reduce*`) | per-character "$\exists s: a_s\cap b_s\neq0$" reduction; used both in `uppass_node` (`:52`) and the union scorer (`:426`). |

**Divergence found:** the comment at `src/ts_fitch.cpp:406-410` documents an invariant
(*superset ⇒ under-count ⇒ lower bound*) that is **false for the finals the function reads**. The
comment's reasoning would be valid only if `tree.final_` held MPR sets; `uppass_node` produces a
strict subset of the MPR set (Steps 2, 5), inverting the containment.

---

## Edge cases

- **$k=2$ (binary), all tips unambiguous singletons, fully resolved — (S4).**
  Caterpillar $\big(((\,t_0,t_1)_{4},\,t_2)_{5},\,t_3\big)_{6}$ with singleton tips
  $t_0{=}\{0\},t_1{=}\{1\},t_2{=}\{0\},t_3{=}\{1\}$ (pattern $0,1,0,1$, homoplastic, length 2).
  Down-pass: $P(4)=\{0,1\}$, $P(5)=\{0\}$, $P(6)=\{0,1\}$. Deployed finals:
  $\mathrm{fin}(6)=\{0,1\},\mathrm{fin}(5)=\{0\},\mathrm{fin}(4)=\{0\},\mathrm{fin}(t_0)=\{0\}$.
  Up-pass: $U\!p(5)=\{1\}$, $U\!p(4)=\{0,1\}$, $U\!p(t_0)=\{1\}$. Edge above $t_0$:
  $E(t_0)=\mathrm{combine}(\{0\},\{1\})=\{0,1\}$.
  Insert leaf $r=\{1\}$ on edge $(4,t_0)$:
  $\text{added}_{\text{exact}}=[\{1\}\cap\{0,1\}=\varnothing]=0$, but
  $U=\mathrm{fin}(4)\cup\mathrm{fin}(t_0)=\{0\}$ so
  $\text{added}_{\text{union}}=[\{1\}\cap\{0\}=\varnothing]=1$. **Full rescore = +0** (attaching
  $\{1\}$ there keeps length 2). Union **over-counts by 1**. Every tip is a definite binary state
  and the tree is fully resolved — so "unambiguous binary ⇒ exact" is **false**. The cause: under
  homoplasy the internal prelim sets ($P(4)=\{0,1\}$) and up sets are non-singleton, and the
  simplified up-pass collapses $\mathrm{fin}(4)=\mathrm{fin}(t_0)=\{0\}$ (an ACCTRAN-like single
  choice), dropping the equally-parsimonious state $1$ that $E(t_0)$ retains.
- **Constant / compatible (homoplasy-free) character.** Then every $P,U\!p,M$ is a singleton and
  $\mathrm{fin}(v)=P(v)=M(v)$, so $U=E(D)$ and both scorers are exact. Exactness holds *only* in the
  no-homoplasy regime, not the "unambiguous" regime.
- **Equality census.** Over 140,044 sampled edges, $E=U$ on $\approx 82\%$; $U\subsetneq E$ on the
  rest; $E\subsetneq U$ and incomparable pairs **never** occurred — consistent with $U\subseteq E$.
  Even restricted to *singleton finals* ($|\mathrm{fin}(A)|=|\mathrm{fin}(D)|=1$), $E\neq U$ in
  21,757 / 104,181 cases, so singleton finals do **not** imply exactness. Per-character exactness
  for *all* $r$ holds iff $E(D)=U$, i.e. iff $M(D)\subseteq\mathrm{fin}(A)\cup\mathrm{fin}(D)$.
- **Degree-2 root special case.** Covered by the $A=\text{root}$ branches of Steps 4–5 and
  `src/ts_fitch.cpp:585-586`.
- **`upweight`/weights.** Enter only as non-negative multipliers (Step 8); do not affect
  containment.
- **Staleness (Assumption 8 relaxed).** The over-count in (S2) is intrinsic to the formula with
  *fresh* finals; it is **not** a staleness artifact. Stale finals (endpoints not re-passed for the
  current clip) add further error of *unsigned* direction on top — so a prior "$\sim$2.75%
  over-count" observation is explained by the formula itself, and staleness cannot be relied on to
  turn it into a lower bound.

---

## Verdict

**Implementation bug — documented invariant is false; the approximation is not a sound lower
bound.**

- The exact scorer `fitch_indirect_length_cached` over $E(D)=\mathrm{combine}(P,U\!p)$ is correct
  (S1; gate 6066/6066).
- The deployed union scorer `fitch_indirect_length` **over-counts**: $U=\mathrm{fin}(A)\cup
  \mathrm{fin}(D)\subseteq E(D)$ (S2, proven; 958,016/958,016 edges), so it is a sound **upper**
  bound, matching the empirical gate's P2 ("over-counts, unsafe"). The in-code comment at
  `src/ts_fitch.cpp:406-410` asserts the opposite containment and the opposite bound direction; it
  is wrong **for the finals the function reads**. It would be correct only for true MPR finals
  (S3 / gate P3), which `uppass_node` does not compute.
- **Deployment consequence.** Do **not** use `fitch_indirect_length` as a pruning/rejection screen
  ("skip candidates whose bound $\ge$ current best"): being an upper bound, it can discard genuinely
  improving moves. Its only safe use is the current one — a heuristic *ranker* whose chosen move is
  **re-scored exactly** afterwards (the `temper` caller); even there it systematically
  over-penalises attachment points where the simplified final dropped an MPR state.
- Not exact on binary/unambiguous data (S4): explicit all-singleton, two-state, fully-resolved
  counterexample; the gate independently found over-counts on 34.8% and true-finals under-counts on
  27.6% of binary candidates.

**Recommendations (reported, not patched — the change is non-trivial):**
1. **Fix the comment** at `src/ts_fitch.cpp:406-414` to state the true direction: with the deployed
   `uppass_node` finals, $\mathrm{final}(A)\cup\mathrm{final}(D)\subseteq E(D)$, so the union
   **over-counts** and is an **upper** bound — safe only for approximate ranking with exact
   re-score, never for pruning. *(This is a documentation-only trivial fix; a one-line patch is
   noted below.)*
2. If a *sound lower bound* is ever wanted (e.g. for admissible pruning), either (a) score against
   the exact edge set $E(D)$ directly (already available, exact, S1), or (b) build the union from
   **true MPR finals** (a proper Fitch/Sankoff up-pass, `final(D)=combine(prelim(D),up(D))`), not
   the `uppass_node` collapse. Option (a) is strictly preferable since $E(D)$ is both exact and
   already computed.

> Caveat (1): whether the deployed `uppass_node` simplification is *also* a defect at its other
> call sites (single-tree final-state reconstruction, ACCTRAN/DELTRAN labelling) is out of lane
> here — this proof only concerns its use inside `fitch_indirect_length`. Should be picked up by a
> follow-up review of `fitch_uppass` consumers / role `math-prover`.
> Caveat (2): the practical search-quality impact of the over-count in the `temper` ranker (does it
> measurably change accepted trajectories or time-to-optimum?) is empirical. Should be picked up by
> role `mcmc-diagnostician` / a heavy A/B on `temper`.

---

### Trivial-fix patch (documentation only)

A one-hunk comment correction is captured at
`dev/red-team/patches/union-of-finals-comment.patch` (**not committed, not merged** — the
orchestrator decides). It replaces the inverted "superset ⇒ under-counts (never over-counts)" note
with the proven direction. No behavioural code changes.
