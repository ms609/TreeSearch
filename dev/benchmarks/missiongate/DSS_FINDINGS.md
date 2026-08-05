# DSS (diameter-limited local sectors) — mission-gate verdict: NO BUILD (mission-dead)

Date 2026-07-16. Branch `claude/dss-probe` (off `cpp-search` 2faf86ca). Mission B (SPEED).
Sibling to rev5 (`FINDINGS.md`, [[sector-resolve-status]]); same shared-T0 missiongate
protocol, same 4 datasets × 3 TNT-`mult` T0 seeds. **Confirmed on Hamilton at n=15.**

## Question (the mission gate)
TNT `dss N D` restricts each sector to a topological neighbourhood of diameter D around a
pivot, instead of an arbitrary clade in a size band. Hypothesis: smaller, locally-bounded
sectors cut per-sector rebuild cost, so sectorial escape becomes cheaper per candidate —
cheap enough to FLIP rev5's ranking (ratchet −1.34 vs sectorial −0.25 steps/Mcand on Zanol)
and DISPLACE the ratchet that is 66–83% of production wall. If it can't beat ratchet, the
thread is mission-dead and that negative is the deliverable.

## Method: config-first, ZERO new C++ (rev4/rev5 discipline; advisor-directed)
The primary driver of "smaller sector" is TIP-COUNT, controlled by the EXISTING
`sectorMaxSize` knob — which rev5 NEVER swept small (baseline held maxS at 0.65n, selectem at
1.0n; only picks/ras/rounds swept). Diameter-limiting is a STRICTLY-more-restrictive geometry
ON TOP of size, so if shrinking maxS doesn't climb toward ratchet, the diameter build cannot
rescue it. So we probed the size axis first, no new source: `smallmax` (maxS ∈ {12,20,30},
min 6, collapse off, ras=3, auto picks, rounds {1,4,16,64}) vs `selectem_large` (maxS=n,
collapse 0.4n, ras=3, AUTO picks) vs `ratchet_ref` (default, single-start), from the fixed TNT
T0, RSS the only escape engine. **CAVEAT on `selectem_large`:** it uses auto picks (~3/round for
the large band), which UNDER-provisions rev5's `selectem` (explicit pk ∈ {5,20}); it is a rough
large-sector reference, NOT a faithful rev5 reproduction (my Zhu selectem_large median −1 vs
rev5's selectem −6). This only weakens the large-sector arm — it does not affect the
smallmax-vs-ratchet verdict below (smallmax auto picks give MANY picks for small bands, e.g.
~16/round at maxS=12, so smallmax is if anything over-provisioned). (A true diameter *ball* needs multi-HTU reduced datasets —
rejected: breaks the root_ok/reinsertion/HTU-anchor invariants the file rests on, for a lever
the priors argue against.) Scripts: `dss_size_probe.R`, `dss_size_aggregate.R`,
`dss_reach_zanol.R`; Hamilton `hamilton_dss_{build,array}.sh`.

## The steps/Mcand axis is BIASED and cannot adjudicate small-vs-large
`candidates_evaluated` is a raw count; it doesn't know a maxS=12 candidate is scored on a
~13-tip REDUCED dataset. Measured per-candidate wall (total wall / total cand): small-sector
candidates are ~2× MORE expensive per wall than ratchet's (per-pick reduced-dataset build +
3× RAS overhead over fewer candidates), so steps/Mcand and steps/WALL rank the arms
DIFFERENTLY, and steps/Mcand systematically FLATTERS small sectors — the direction of the
apparent result. rev5 could lean on steps/Mcand (its selectem was near ratchet's candidate
scale AND the bias ran against its negative). We adjudicate on WALL (the true mission axis,
[[mission-wallclock-to-optimum]]) + REACH.

## HAMILTON n=15 RESULTS (authoritative; array 17898332, build 17898331, HEAD 023b0389)
Escape<0 = better than T0. st/s = escape/wall (wall axis). st/Mc = escape/Mcand [BIASED].

| dataset (target) | ratchet esc_med / reach / st·s⁻¹ / st·Mc⁻¹ | best smallmax st·s⁻¹ (esc_med) | best smallmax st·Mc⁻¹ | beats ratchet@wall? |
|---|---|---|---|---|
| Zanol2014 (−13) | −12 / −13 / **−25.4** / **−1.90** | −17.8 (esc −5) | −1.78 | NO |
| Zhu2013 (−8)    | −3 / −9 / −13.8 / −0.63 | −6.2 (esc −1) | −0.44 | NO |
| Wortley2006 (−5)| −7 / −15 / −86.4 / −6.76 | faster but shallower median | −9.8 | NO |
| Giles2015 (−4)  | −3 / −5 / −13.2 / −0.55 | ~−5 (esc 0) | −0.44 | NO |

**VERDICT: 0 / 4 datasets — no small-sector config has better MEDIAN escape than ratchet at
≤ its wall** (esc_med vs esc_med, same 3 seeds pooled for both arms → apples-to-apples). On
the HARD target (Zanol) ratchet dominates on EVERY axis, including the biased steps/Mcand
(−1.90 beats every smallmax; best smallmax −1.78). On easy datasets small sectors are
faster-per-step (higher st/s) but reach a SHALLOWER median, so they never satisfy the gate
(reach ratchet's depth cheaper). Small sectors' median escape is uniformly worse than
ratchet's. NB the per-config `reach`/`esc_best` column is NOT load-bearing: the aggregate's
"target" pools 3 seeds with different starts, so treat esc_best only as a loose best-of-15,
not a clean target comparison. The verdict rests on the esc_med-vs-ratchet-at-wall gate.

## The n=3 local pilot's apparent LEAD was noise + metric bias — do not trust it
The local pilot (n=3, Zanol s1 / Zhu s1) showed smallmax maxS=12 r=1 = −2.39 st/Mcand
"beating" ratchet −1.51, and a Zhu maxS=30 wall win. At n=15 BOTH evaporate: Zanol maxS=12 r=1
falls to −1.19 st/Mc (median escape −4, not the pilot's −8); Zhu ratchet median is −3 and every
smallmax median is −1. Two independent reasons the pilot misled — (i) steps/Mcand is biased
toward small sectors, (ii) n=3 single-seed is far too thin (rev5 flagged Zhu ratchet as a
coin-flip). Fixing the metric AND getting proper n both had to happen; the negative survives
both. (This is why we did NOT build or dispatch on the n=3 "lead".)

## Verdict: NEGATIVE — dss mission-dead, NO build
On the wall axis and on reach, small-local sectors do not beat ratchet (0/4), and on the hard
target ratchet wins even on the biased metric. Diameter-limiting is strictly MORE local →
shallower/less-competitive → moves the wrong way. **No C++ diameter knob is warranted.** This
CONFIRMS and EXTENDS rev5 ("sectorial Pareto-dominated by ratchet; 0/4 escape cheaper") to the
small-maxS regime rev5 never measured — the whole size axis is now covered, and rev5's verdict
holds across it. **NOT tested here:** rev5's own handoff (large-band `selectem` at pk 5-20,
ras≥3 being lower-VARIANCE than a single ratchet start on ratchet-unreliable Zhu/Giles). That
is a different config (large band, explicit picks) and a variance question; my `selectem_large`
used auto picks (~3) so it under-provisions that config and cannot speak to it. The handoff
remains an OPEN budget/recipe question for another session, exactly as rev5 left it — the dss/
speed question does not bear on it. Sectorial-escape / dss (small-sector) thread CLOSED.
