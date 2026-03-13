# TreeSearch – Posit Assistant Memory

## Project

R package for phylogenetic tree search and character analysis.
Source: `C:/Users/pjjg18/GitHub/TreeSearch`

## LocalConcordance

### Function

`LocalConcordance(dataset, sigma, ...)` in `R/LocalConcordance.R`.
Computes Gaussian-weighted mean NMI between alignment sites, returning a
matrix (rows = alignment columns, cols = sigma values).

### Key parameters

| Parameter | Default | Notes |
|---|---|---|
| `sigma` | `c(0.5, 0.8, 1, 1.2, 1.5, 2, 3, 5)` | σ ≈ 1.5 ≈ PhyIN `-d 2`; weight at d=σ is exp(−½) ≈ 61% of d=1 weight |
| `internal_gaps` | `TRUE` | Recode internal gaps as real state before NMI; rescues informative indels |
| `block_size` | `NULL` | Integer: activates rolling-window noise scan (PhyIN `-b` analogue) |
| `threshold` | `0.5` | Min proportion of noisy sites in block to trigger trim (PhyIN `-p`) |
| `noise_level` | `0` | NMI threshold for classifying a site as noisy (natural with `normalize=TRUE`) |

### Sigma ↔ PhyIN -d correspondence

PhyIN uses a uniform hard-cutoff window; our Gaussian is "soft".
Most direct analogue: **σ ≈ 1.5 ≈ PhyIN `-d 2`** (weight at d=2 is ~50% of d=1,
negligible beyond d=3).  σ=1 ≈ d=1; σ=2 ≈ d=2–3.

### Case study dataset

Harmochirina jumping spiders (Maddison UCE data).

| Object | Description |
|---|---|
| `harmo` | Single UCE locus: 48 taxa × 2138 sites, 1122 unique patterns |
| `harmo_trim` | PhyIN-filtered version of `harmo`: 1444 sites (694 removed, 32.5%) |
| `harmo_full` | Full concatenated alignment: 48 taxa × 122,136 sites |
| `kept` | Logical (length 122,136): PhyIN kept/removed for full alignment |
| `phyin_removed` | Logical (length 2138): per-site PhyIN decisions for `harmo` locus (reconstructed via sequential column matching of `harmo_trim` against `harmo`) |
| `trim_15` | Logical (length 2138): LC block-scan trim at σ=1.5, block=10, threshold=0.5 |

Load data: `harmo` was loaded from the Harmochirina UCE FASTA; reconstruct
`phyin_removed` via sequential column matching (see session history).

### Key findings

- LC (block scan) and PhyIN have **zero overlap** in trimmed sites.
- LC trims the **sparse left flank** (cols ~1–650): NA/near-zero NMI due to
  only 1–2 taxa with data. PhyIN's first block (cols 1–100, 1–2 taxa) is
  almost certainly a **coverage filter**, not a compatibility failure —
  two-taxon sites cannot exhibit all four gametes.
- PhyIN trims **scattered sites in the informative core** (cols 500–2100);
  these have median NMI ≈ 0.9 at σ=1.5 — high local concordance but
  globally incompatible characters.
- `internal_gaps = TRUE` rescues ~122 sites (6%) that without gap recoding
  would be flagged as noisy; their NMI jumps from median −0.06 → +0.62.

### Useful visualisation

Four-panel base-graphics plot: **NMI trace on top**, then PhyIN strip,
LC strip, alignment image. Use `layout()` with identical `mar[2]` and `mar[4]`
across all panels — this is the only reliable way to guarantee pixel-perfect
column alignment. Key elements:

```r
# Measure left margin in inches from actual label widths BEFORE layout:
op <- par(mar = c(0, 0, 0, 0)); plot.new()
left_mai  <- max(strwidth(taxa_labels, units = "inches", cex = 0.42)) + par("cin")[1] * 0.5
right_mai <- par("cin")[2] * 1.5
par(op)

layout(matrix(1:4), heights = c(3.5, 1, 1, 13))
# Use par(mai = ...) in inches (NOT par(mar = ...) in lines) for all panels —
# this is the only reliable way to guarantee pixel-perfect column alignment.
# CRITICAL: use xaxs = "i" in every plot.window() call.
#   image() defaults to xaxs = "i" (no axis expansion).
#   plot.window() defaults to xaxs = "r" (4% expansion each side).
#   Without xaxs = "i", strip and NMI panels are shifted ~4% right vs alignment.
#
# NMI panel: plot.window(xlim, ylim, xaxs="i", yaxs="i")
#   shade NA regions (grey92) before drawing traces;
#   draw σ_lo and σ_hi faintly (#BBBBBB), main σ in grey20.
# Strips: plot.window(xlim, ylim, xaxs="i"); rle() the logical trim vector,
#   then rect() per run.
# Alignment: image(x=sites, y=taxa, z=t(aln_num)[, nTaxa:1])
#   with state_cols = c("grey85","#4DAF4A","#377EB8","#FF7F00","#E41A1C")
# Drop taxa with <1% coverage before plotting.
```

### Tests

37 tests in `tests/testthat/test-LocalConcordance.R`, all passing.
`stats::filter` must be qualified to avoid masking by `dplyr::filter`.
