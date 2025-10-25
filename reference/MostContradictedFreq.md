# Frequency of most common contradictory split

`MostContradictedFreq()` counts the occurrences of the single split that
most frequently contradicts each split in `tree`.

This function was written during a code sprint: its documentation and
test cases have not yet been carefully scrutinized, and its
implementation may change without notice. Please alert the maintainer to
any issues you encounter.

## Usage

``` r
MostContradictedFreq(tree, forest)
```

## Arguments

- tree:

  A tree of class [`phylo`](https://rdrr.io/pkg/ape/man/read.tree.html).

- forest:

  a list of trees of class `phylo`, or a `multiPhylo` object; or a
  `Splits` object.

## Value

`MostContradictedFreq()` returns, for each split in `tree`, the number
of times that its most common contradictory split occurs in `forest`.

## Details

Goloboff et al. (2003) propose comparing the frequency of a split in a
resampled population with the frequency of the most common contradictory
split. This measure contributes to the "groups present / contradicted"
score.

## See also

[`PresCont()`](https://ms609.github.io/TreeSearch/reference/PresCont.md)
calculates the "groups present / contradicted" score.

Other split support functions:
[`JackLabels()`](https://ms609.github.io/TreeSearch/reference/JackLabels.md),
[`Jackknife()`](https://ms609.github.io/TreeSearch/reference/Jackknife.md),
[`MaximizeParsimony()`](https://ms609.github.io/TreeSearch/reference/MaximizeParsimony.md),
[`PresCont()`](https://ms609.github.io/TreeSearch/reference/PresCont.md),
[`SiteConcordance`](https://ms609.github.io/TreeSearch/reference/SiteConcordance.md)
