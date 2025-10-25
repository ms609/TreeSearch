# Tree suboptimality

How suboptimal is a tree?

## Usage

``` r
Suboptimality(trees, proportional = FALSE)
```

## Arguments

- trees:

  list of trees, to include an optimal tree

- proportional:

  logical stating whether to normalise results to lowest score

## Value

`Suboptimality()` returns a vector listing, for each tree, how much its
score differs from the optimal (lowest) score.
