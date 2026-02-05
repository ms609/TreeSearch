# Generate colour to depict the amount and quality of observations

Generate colour to depict the amount and quality of observations

## Usage

``` r
QACol(amount, quality)

QCol(amount, quality)

QALegend(
  where = c(0.1, 0.3, 0.1, 0.3),
  n = 5,
  Col = QACol,
  xlab = "Amount →",
  ylab = "Quality →",
  ...
)
```

## Arguments

- amount:

  Numeric vector of values between 0 and 1, denoting the relative amount
  of information

- quality:

  Numeric vector of values between -1 and 1, denoting the quality of
  observations, where 0 is neutral.

- where:

  Location of legend, passed to `par(fig = where)`

- n:

  Integer vector giving number of cells to plot in swatch for `quality`
  and `amount`.

- Col:

  Function that takes vectors `amount` and `quality` and returns a
  vector of colours. QCol colours by data quality (concordance); QACol
  by quality and amount of information.

- xlab:

  Character giving a label for the x axis.

- ylab:

  Character giving a label for the y axis.

- ...:

  Additional parameters to
  [`mtext()`](https://rdrr.io/r/graphics/mtext.html).

## Value

`QACol()` returns an RGB hex code for a colour, where lighter colours
correspond to entries with a higher `amount`; unsaturated colours denote
a neutral `quality`; and red/cyan colours denote low/high `quality`.

`QCol()` returns an RGB hex code for a colour, where darker, unsaturated
colours denote a neutral `quality`; and red/cyan colours denote low/high
`quality`. `amount` is ignored.

## See also

Other utility functions:
[`ClusterStrings()`](https://ms609.github.io/TreeSearch/dev/reference/ClusterStrings.md),
[`QuartetResolution()`](https://ms609.github.io/TreeSearch/dev/reference/QuartetResolution.md),
[`WhenFirstHit()`](https://ms609.github.io/TreeSearch/dev/reference/WhenFirstHit.md)

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
amount <- runif(80, 0, 1)
quality <- runif(80, -1, 1)
plot(amount, quality, col = QACol(amount, quality), pch = 15)
abline(h = 0)
```
