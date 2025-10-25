# Set and get the character weightings associated with a Morphy object.

`MorphyWeights()` details the approximate and exact weights associated
with characters in a `Morphy` object; `SetMorphyWeights()` edits them.

## Usage

``` r
MorphyWeights(morphyObj)

SetMorphyWeights(weight, morphyObj, checkInput = TRUE)
```

## Arguments

- morphyObj:

  Object of class `morphy`, perhaps created with
  [`PhyDat2Morphy()`](https://ms609.github.io/TreeSearch/reference/PhyDat2Morphy.md).

- weight:

  A vector listing the new weights to be applied to each character

- checkInput:

  Whether to sanity-check input data before applying. Defaults to `TRUE`
  to protect the user from crashes.

## Value

`MorphyWeights()` returns a data frame with two named rows and one
column per character pattern: row 1, `approx`, is a list of integers
specifying the approximate (integral) weights used by MorphyLib; row 2,
`exact`, is a list of numerics specifying the exact weights specified by
the user.

`SetMorphyWeights()` returns the Morphy error code generated when
applying `weight`.

## See also

Other Morphy API functions:
[`GapHandler()`](https://ms609.github.io/TreeSearch/reference/GapHandler.md),
[`MorphyErrorCheck()`](https://ms609.github.io/TreeSearch/reference/MorphyErrorCheck.md),
[`PhyDat2Morphy()`](https://ms609.github.io/TreeSearch/reference/PhyDat2Morphy.md),
[`SingleCharMorphy()`](https://ms609.github.io/TreeSearch/reference/SingleCharMorphy.md),
[`UnloadMorphy()`](https://ms609.github.io/TreeSearch/reference/UnloadMorphy.md),
[`is.morphyPtr()`](https://ms609.github.io/TreeSearch/reference/is.morphyPtr.md),
[`mpl_apply_tipdata()`](https://ms609.github.io/TreeSearch/reference/mpl_apply_tipdata.md),
[`mpl_attach_rawdata()`](https://ms609.github.io/TreeSearch/reference/mpl_attach_rawdata.md),
[`mpl_attach_symbols()`](https://ms609.github.io/TreeSearch/reference/mpl_attach_symbols.md),
[`mpl_delete_Morphy()`](https://ms609.github.io/TreeSearch/reference/mpl_delete_Morphy.md),
[`mpl_first_down_recon()`](https://ms609.github.io/TreeSearch/reference/mpl_first_down_recon.md),
[`mpl_first_up_recon()`](https://ms609.github.io/TreeSearch/reference/mpl_first_up_recon.md),
[`mpl_get_charac_weight()`](https://ms609.github.io/TreeSearch/reference/mpl_get_charac_weight.md),
[`mpl_get_gaphandl()`](https://ms609.github.io/TreeSearch/reference/mpl_get_gaphandl.md),
[`mpl_get_num_charac()`](https://ms609.github.io/TreeSearch/reference/mpl_get_num_charac.md),
[`mpl_get_num_internal_nodes()`](https://ms609.github.io/TreeSearch/reference/mpl_get_num_internal_nodes.md),
[`mpl_get_numtaxa()`](https://ms609.github.io/TreeSearch/reference/mpl_get_numtaxa.md),
[`mpl_get_symbols()`](https://ms609.github.io/TreeSearch/reference/mpl_get_symbols.md),
[`mpl_init_Morphy()`](https://ms609.github.io/TreeSearch/reference/mpl_init_Morphy.md),
[`mpl_new_Morphy()`](https://ms609.github.io/TreeSearch/reference/mpl_new_Morphy.md),
[`mpl_second_down_recon()`](https://ms609.github.io/TreeSearch/reference/mpl_second_down_recon.md),
[`mpl_second_up_recon()`](https://ms609.github.io/TreeSearch/reference/mpl_second_up_recon.md),
[`mpl_set_charac_weight()`](https://ms609.github.io/TreeSearch/reference/mpl_set_charac_weight.md),
[`mpl_set_num_internal_nodes()`](https://ms609.github.io/TreeSearch/reference/mpl_set_num_internal_nodes.md),
[`mpl_set_parsim_t()`](https://ms609.github.io/TreeSearch/reference/mpl_set_parsim_t.md),
[`mpl_translate_error()`](https://ms609.github.io/TreeSearch/reference/mpl_translate_error.md),
[`mpl_update_lower_root()`](https://ms609.github.io/TreeSearch/reference/mpl_update_lower_root.md),
[`mpl_update_tip()`](https://ms609.github.io/TreeSearch/reference/mpl_update_tip.md),
[`summary.morphyPtr()`](https://ms609.github.io/TreeSearch/reference/summary.morphyPtr.md)

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
tokens <- matrix(c(
  0, 0, 0, 1, 1, 2,
  0, 0, 0, 0, 0, 0), byrow = TRUE, nrow = 2L,
  dimnames = list(letters[1:2], NULL))
pd <- TreeTools::MatrixToPhyDat(tokens)
morphyObj <- PhyDat2Morphy(pd)
MorphyWeights(morphyObj)
#>        [,1] [,2] [,3]
#> approx 3    2    1   
#> exact  3    2    1   
if (SetMorphyWeights(c(1, 1.5, 2/3), morphyObj) != 0L) message("Errored")
MorphyWeights(morphyObj)
#>        [,1] [,2] [,3]     
#> approx 6    9    4        
#> exact  1    1.5  0.6666667
morphyObj <- UnloadMorphy(morphyObj)
```
