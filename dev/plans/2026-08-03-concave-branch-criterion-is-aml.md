# The concave-branch criterion is Ancestral Maximum Likelihood — decision note

**Date:** 2026-08-03. **Branch:** `claude/maximum-information-tree`.
**Outcome: do not build it as an optimality criterion.** One open experiment
(§7) and one cheap reusable by-product (§6) survive.

## 0. Where this came from: the published criterion, refuted

Nishimaki & Sato (2026, *J Mol Evol* 94:726–735, `TUS-Satolab/PhyloInfo`) score a
rooted tree as the **mean over its $2(n-1)$ edges of $\mathrm{MI}(\text{parent};\text{child})/S(\text{parent})$**,
with ancestral sequences reconstructed by a distance-weighted plurality vote over
*all* leaves ($w_{ij} = L_{\max} - L_{ij}$). **Do not build a search on it.** Four
independent killers, all verified numerically:

1. **Consensus collapse.** $w_{ij}$ decays *linearly* in path length while leaf
   counts grow *exponentially* with depth, so distant shells dominate the weight
   mass. A node's own two daughters hold 75/50/31.8/19.6/11.7/6.8/3.9/1.2 % of its
   reconstruction weight at $n=4/8/16/32/64/128/256/1024$. At the paper's own
   34–121 taxa that is ~4–12%: all reconstructions converge on one global
   consensus, internal edges pin at 1.0, and the residual score is *exactly
   topology-independent*. Not repairable by rooting, normalisation or chance
   correction.
2. **Invariant-site inflation.** Holding informative content fixed at 4 sites, the
   edge score runs 0.311 → 0.894 as the invariant fraction runs 0 → 0.94.
   Invariant fraction 0.733 reproduces the paper's headline 0.83354, so the
   published numbers largely measure *alignment conservation*; the score is
   gameable by padding and is not comparable across matrices. Also explains the
   tiny 0.0033 NJ/ML/BI spread.
3. **Label invariance.** MI is invariant to any bijection of either alphabet, so a
   child that is the exact *complement* of its parent scores 1.000 — identical to
   zero changes — while 6/60 scattered changes scores 0.522. Wrong-direction rate
   4.5% overall, **12.6% on already-divergent edges**, which is the region a
   hill-climber explores. Step 1's plurality tie-break is not label-invariant and
   no deterministic tie-break can be, so under unit lengths (the morphological
   case) the score depends on arbitrary 0/1/2 coding.
4. **What it actually measures.** At leading order the objective reduces to
   $\tfrac12 + \tfrac12\,\mathrm{mean}_j\,\mathrm{MI}(\text{recon parent}_j;\text{leaf}_j)/H(\text{recon parent}_j)$
   — a *local sister-agreement statistic*. That is why the paper's SPR-neighbourhood
   ranking succeeds and why it implies nothing about global maximisation.

For an $m$-state matrix the criterion decomposes exactly as
$\log(m-1)\,L(T) + \sum_e k H_b(d_e/k)$ — parsimony length *plus* a concave
per-branch penalty. Chance-correcting it (Cohen's $\kappa$) fixes the pathology but
$\kappa$ is linear in the change count ($R^2 = 0.999999$), so it collapses onto
weighted parsimony: the repair and the novelty vanish together. The remaining exit
was the concave penalty itself — which is what the rest of this note is about, and
which also fails. Prior art for the general programme: Wheeler & Varón 2025,
*Cladistics* 41(2):193–211, doi 10.1111/cla.12603 (PMDL).

## 1. The criterion

Give every edge its own free substitution probability under a symmetric
$m$-state (Mk) process, maximize the likelihood over topology, edge
probabilities **and ancestral states**, and take negative log-likelihood:

$$\text{minimize}\quad \Phi(T) \;=\; \min_{\text{reconstruction}} \sum_{e \in E(T)} f_m(d_e)$$

where $d_e$ is the number of characters changing on edge $e$, $k$ the number of
characters, and

$$f_m(d) \;=\; \begin{cases} k\,H_b\!\left(\dfrac{d}{k}\right) + d\log(m-1), & \dfrac{d}{k} \le \dfrac{m-1}{m}\\[2ex] k\log m, & \text{otherwise}\end{cases}$$

with $H_b$ the binary entropy in nats. Verified: $f_m(d)$ equals the negated
maximized log-likelihood of a symmetric Mk edge whose change probability is
constrained to the identifiable range $p \le (m-1)/m$, to
$\mathbf{5.7\times10^{-14}}$ across 30 $(m,k,d)$ cells.

### Corrected clamp

An earlier draft wrote $f_m(d) = k H_b\!\left(\min\!\left(\frac{d}{k},\frac{m-1}{m}\right)\right) + d\log(m-1)$.
**That is wrong for $m>2$**: it clamps the entropy term but lets the linear term
keep growing. The constrained MLE sits at $p=(m-1)/m$, where
$1-p = p/(m-1) = 1/m$ and so $-\ell = k\log m$ exactly — a *constant*. The
piecewise form above is continuous at the join (verified algebraically and
numerically) and has the clean reading that **an edge can never cost more than
$k\log m$**, the cost of describing the child from scratch. For $m=2$ the two
forms coincide because $\log(m-1)=0$, which is why the error went unnoticed.

### Shape (verified, $k=200$)

Monotone non-decreasing on $0..k$ for $m=2,3,5$; concave below the clamp; caps
exactly at $k\log m$ (138.629 / 219.722 / 321.888). Clamp bites at $d=101/134/161$.

## 2. This is a named criterion: AML

$\Phi$ is **Ancestral Maximum Likelihood**, studied since Barry & Hartigan
(1987). Not novel. The relevant results:

- **NP-hard.** Addario-Berry, Chor, Hallett, Lagergren, Panconesi & Steel
  (2004), *Ancestral maximum likelihood of evolutionary trees is hard*. This
  settles the tractability question raised below in §4 — negatively.
- **Approximable.** Alon, Chor, Pardi & Rapoport (2010), IEEE/ACM TCBB
  7:183–187. Reformulate two-state MP and two-state AML as Steiner tree in a
  hypercube; **16/9** approximation for AML, asymptotically 1.55 for MP. They
  use the same two-state entropy function.
- **Inconsistent.** *Shrinkage Effect in Ancestral Maximum Likelihood*
  (arXiv:0802.0914): AML "can 'shrink' short edges in a tree, resulting in a
  tree that has no internal resolution as the sequence length grows."
- **Positioned.** The literature places AML "in between" MP and ML, with some
  properties of each. Goldman (1990, *Syst. Zool.* 39:345–361) is the ancestor:
  with $p$ *fixed and equal across branches*, maximising over reconstructions
  gives exactly parsimony. Freeing $p$ per branch is what produces $\Phi$.

## 3. Why it is the transpose of implied weighting

Let $x_{je}=1$ iff character $j$ changes on edge $e$. Then

| criterion | form | separability |
|---|---|---|
| equal weights | $\sum_j\sum_e x_{je}$ | both indices |
| implied weighting | $\sum_j f\!\left(\sum_e x_{je}\right)$ | concave over *edges within a character*; the DP separates by character ⇒ **easy** |
| this criterion | $\sum_e f\!\left(\sum_j x_{je}\right)$ | concave over *characters within an edge*; the DP separates by character ⇒ **coupled ⇒ hard** |

IW and AML apply the same concave transform along opposite axes of the same
matrix. IW is cheap for exactly the reason AML is NP-hard. This is the crispest
statement of the difficulty and it is worth keeping even though the criterion is
not being built.

## 4. Algorithmics

**MM / concave linearisation.** $f$ concave ⇒ its tangent is a global
over-estimator ⇒ a valid majoriser. Linearising at the current $d^{(t)}$ gives
per-edge weights

$$\lambda_e \;=\; f_m'(d_e) \;=\; \log\frac{(k-d_e)(m-1)}{d_e}$$

— the **empirical log-odds against change on that edge** — and minimising the
majoriser is exactly **edge-weighted Sankoff**, solvable per character by the
standard DP. Iterating gives monotone descent
$\Phi(d^{(t+1)}) \le M(d^{(t+1)}\,|\,d^{(t)}) \le M(d^{(t)}\,|\,d^{(t)}) = \Phi(d^{(t)})$.
Verified: tangent-below-$f$ violations 0, subadditivity violations 0;
$\lambda = 5.293/4.595/3.664/2.197/1.099$ at $d = 1/2/5/20/50$ ($k=200$, $m=2$).

$\lambda_e \ge 0$ iff $d_e \le k(m-1)/m$ — **the identifiability clamp is exactly
non-negativity of the MM weight.** Pleasing, and a useful internal check.

**Exact minima sit at pattern-class vertices.** The feasible set of $d$-vectors
is a Minkowski sum over pattern classes; extreme points of a Minkowski sum are
sums of extreme points, and the polytope of a class of $n_c$ identical
characters has vertices $n_c p$ for $p$ a single reconstruction. A concave
minimum sits at a vertex, and that vertex is integral. Hence **the exact optimum
assigns every character of a pattern class to one reconstruction**, so exactness
scales with the number of *distinct patterns*, not with $k$. Used to make the
4-taxon tests below exact-by-construction rather than heuristic.

**No usable lower bound.** Subadditivity plus monotonicity give
$\Phi \ge f_m(L(T)) $ and $\Phi \ge \frac{f_m(k)}{k}L(T) = L(T)\log 2$ for $m=2$;
both are ~7–14× loose at realistic sizes, so branch-and-bound is out and there
is no exactness certificate.

**$f'(0)=\infty$.** An edge with no changes has infinite resistance to acquiring
one, so MM freezes and can never move a change onto an empty edge.
Krichevsky–Trofimov / Jeffreys smoothing $\hat p = (d+\frac12)/(k+1)$ removes the
singularity ($\lambda(0)=\log(2k+1)$), is the proper MDL two-part code rather
than a hack, and adds a $\frac12\log k$ per-edge parameter cost that correctly
penalizes resolution. **This is not merely an algorithmic detail — it is the
mechanism of the published shrinkage pathology** (§2), so the smoothed variant
was tested as a candidate repair (§5).

## 5. Measured behaviour

Exact 4-taxon tests. True tree $((1,3),(2,4))$ with tips 1 and 2 long and
**non-sister** — the Felsenstein zone. Internal $p = 0.03$.

**(a) Analytic zone map** ($m=2$, $k=200$; $s$ convergent changes shared by the
long tips, $u$ private each; true tree costs $2f(s+u)$, wrong tree
$f(s)+2f(u)$):

| $s=u$ | $d/k$ | concave picks | parsimony picks |
|---|---|---|---|
| 2, 5, 10, 20 | 0.010–0.100 | wrong | wrong |
| 25, 30, 40, 50 | 0.125–0.250 | **TRUE** | wrong |

Crossover $q^\ast = 0.1215$: concavity flips the answer once $\gtrsim 12\%$ of
characters change on the long pendants. Parsimony is wrong at every level.

**(b) Simulated under Mk, topology recovered (8 reps/cell):**

| $k$ | long-branch $p$ | parsimony | AML | AML+KT |
|---|---|---|---|---|
| 200 | 0.10 | **8/8** | 4/8 | 6/8 |
| 200 | 0.30 | 0/8 | **5/8** | **5/8** |
| 1000 | 0.10 | 8/8 | 8/8 | 8/8 |
| 1000 | 0.30 | 0/8 | 5/8 | **7/8** |
| 5000 | 0.10 | 8/8 | 8/8 | 8/8 |
| 5000 | 0.30 | 0/8 | **7/8** | 6/8 |

Two findings, opposite in sign:

- **AML rescues the Felsenstein zone.** Parsimony is 0/8 at *every* $k$ —
  textbook inconsistency. AML reaches 5–7/8.
- **AML damages the easy regime at morphological $k$.** At $k=200$, $p=0.10$
  parsimony is 8/8 and AML is 4/8 — *chance*. It recovers by $k=1000$. KT
  smoothing partly repairs it (6/8) and helps 3 of 6 cells overall.

**(c) Internal-edge bias** (median changes placed on the true short internal
edge; expected $0.03k$): at $p=0.30$, $k=5000$, truth 150, parsimony 127.5, AML
**413.5**. Concavity *extremises* — it pushes each edge toward 0 or toward
saturation. The published shrinkage result is the $\to 0$ half; this is the
$\to$ saturation half. Same mechanism, so AML edge lengths are bimodal and
badly biased in whichever direction local geometry favours.

Caveat: 6–8 reps per cell. The 0/8-vs-5/8 LBA contrast is solid; the
8/8-vs-4/8 easy-regime contrast is suggestive and would need more reps to quote.

## 6. What survives

**A diagnostic, not a criterion.** One MM step from the ordinary parsimony
reconstruction yields $\lambda_e = \log[(k-d_e)(m-1)/d_e]$ per edge — a
principled, zero-new-kernel flag for rate-outlier branches on an existing MPT.
It reuses the reconstruction the search already computes and answers "which
branches in this tree are saturated?" with a likelihood-ratio interpretation.

## 7. The one open question

The four-taxon test **cannot** compare against implied weighting: with binary
characters on four taxa every $h_j \in \{0,1\}$, so IW is a monotone function of
parsimony length and is *provably identical to equal weights*. So these results
do not show AML beats IW — nor that it is redundant with it.

**The decisive experiment is AML vs implied weighting vs equal weights on a
$\ge 6$-taxon Felsenstein-zone simulation**, which is the smallest case where
both correction mechanisms have room to act. Goloboff (2019, *Cladistics*
**36**:1, "Likelihood approximations of implied weights parsimony can be
selected over the Mk model by the Akaike information criterion") already gives
IW its likelihood/AIC footing, so the comparison is well posed and, framed as
"branch-IW vs character-IW", is a genuinely unoccupied question. It is a paper,
not a package feature.

## 8. Why the line closes

The chain runs: NS-2026 information metric → chance-corrected agreement →
weighted parsimony → concave-branch penalty → **AML** → NP-hard *and*
statistically inconsistent → and the principled repair for inconsistency is to
*integrate out* ancestral states rather than maximize over them, which lands on
ordinary Mk ML, already implemented consistently elsewhere. The incidental-
parameters problem (Neyman–Scott 1948) is structural here: the number of
ancestral states grows with $k$, so no amount of smoothing makes joint
estimation consistent.

There is no room on this axis for a criterion that is simultaneously novel,
computable and consistent.
