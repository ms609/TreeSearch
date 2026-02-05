# Number of trees with *m* steps

Calculate the number of trees in which Fitch parsimony will reconstruct
*m* steps, where *a* leaves are labelled with one state, and *b* leaves
are labelled with a second state.

## Usage

``` r
Carter1(m, a, b)

Log2Carter1(m, a, b)

LogCarter1(m, a, b)
```

## Arguments

- m:

  Number of steps.

- a, b:

  Number of leaves labelled `0` and `1`.

## Details

Implementation of theorem 1 from Carter et al. (1990)

## References

Carter M, Hendy M, Penny D, Székely LA, Wormald NC (1990). “On the
distribution of lengths of evolutionary trees.” *SIAM Journal on
Discrete Mathematics*, **3**(1), 38–47.
[doi:10.1137/0403005](https://doi.org/10.1137/0403005) .

See also:

Steel MA (1993). “Distributions on bicoloured binary trees arising from
the principle of parsimony.” *Discrete Applied Mathematics*, **41**(3),
245–261.
[doi:10.1016/0166-218X(90)90058-K](https://doi.org/10.1016/0166-218X%2890%2990058-K)
.

Steel M, Charleston M (1995). “Five surprising properties of
parsimoniously colored trees.” *Bulletin of Mathematical Biology*,
**57**(2), 367–375.
[doi:10.1016/0092-8240(94)00051-D](https://doi.org/10.1016/0092-8240%2894%2900051-D)
.

(Steel M, Goldstein L, Waterman MS (1996). “A central limit theorem for
the parsimony length of trees.” *Advances in Applied Probability*,
**28**(4), 1051–1071.
[doi:10.2307/1428164](https://doi.org/10.2307/1428164) . )

## See also

Other profile parsimony functions:
[`PrepareDataProfile()`](https://ms609.github.io/TreeSearch/dev/reference/PrepareDataProfile.md),
[`StepInformation()`](https://ms609.github.io/TreeSearch/dev/reference/StepInformation.md),
[`WithOneExtraStep()`](https://ms609.github.io/TreeSearch/dev/reference/WithOneExtraStep.md),
[`profiles`](https://ms609.github.io/TreeSearch/dev/reference/profiles.md)

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
# The character `0 0 0 1 1 1`
Carter1(1, 3, 3) # Exactly one step
#> [1] 9
Carter1(2, 3, 3) # Two steps (one extra step)
#> [1] 54

# Number of trees that the character can map onto with exactly _m_ steps
# if non-parsimonious reconstructions are permitted:
cumsum(sapply(1:3, Carter1, 3, 3))
#> [1]   9  63 105

# Three steps allow the character to map onto any of the 105 six-leaf trees.
```
