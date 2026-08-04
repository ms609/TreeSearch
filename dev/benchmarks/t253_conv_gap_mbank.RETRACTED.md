# ⛔ `t253_conv_gap_mbank.csv` IS RETRACTED (2026-08-03)

A sidecar rather than a header, because a comment line inside the CSV would break every
reader that currently parses it.

**Do not use `t253_conv_gap_mbank.csv` as evidence.** Its rows are derived from the
`t252_mbank_*` CSVs (2026-03-27), produced by an engine whose Wagner addition was worse on
**every matrix in the sample** — so the "convergence gap" it tabulates is substantially the
addition bug rather than a property of the datasets.

Measured, Hamilton job 18183151 (`t253_wagner_era_probe.R`, results in
`t253_wagner_era_decision.csv`): bare `AdditionTree`, no search, identical preprocessing for
both engines, 3 seeds each, every tree re-scored by a single scorer, and a provenance
assertion confirming the two arms really used different libraries.

| | |
|---|---|
| matrices where the March engine built a longer tree | **25 of 25** |
| matrices >10% longer | **23 of 25** |
| median ratio (March ÷ current) | **1.365** |
| range | 1.050 → 3.232 |

The probe was built to decide *annotate one row* vs *retract the analysis*. 25 of 25 answers
retract: `project4284` (3.23×) is the extreme of a universal effect, not an outlier to be
excused.

Full reasoning, and what survives, in the retraction box at the top of
`t253_gap_characterization.md`. Note that its `t265` half — 8 named datasets, TNT vs
TreeSearch at 120 s — does **not** depend on the t252 engine and stands.
