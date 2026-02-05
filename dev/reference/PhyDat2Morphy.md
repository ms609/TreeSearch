# Initialize a Morphy object from a `phyDat` object

Creates a new Morphy object with the same size and characters as the
`phyDat` object. Once finished with the object, it should be destroyed
using
[`UnloadMorphy()`](https://ms609.github.io/TreeSearch/dev/reference/UnloadMorphy.md)
to free the allocated memory.

## Usage

``` r
PhyDat2Morphy(phy, gap = "inapplicable", weight = attr(phy, "weight"))
```

## Arguments

- phy:

  An object of phangorn class `phyDat`.

- gap:

  An unambiguous abbreviation of `inapplicable`, `ambiguous` (=
  `missing`), or `extra state`, specifying how gaps will be handled.

- weight:

  Numeric giving weight to apply to each character. Will be recycled.

## Value

`PhyDat2Morphy()` returns a pointer to an initialized Morphy object.

## See also

Other Morphy API functions:
[`GapHandler()`](https://ms609.github.io/TreeSearch/dev/reference/GapHandler.md),
[`MorphyErrorCheck()`](https://ms609.github.io/TreeSearch/dev/reference/MorphyErrorCheck.md),
[`MorphyWeights()`](https://ms609.github.io/TreeSearch/dev/reference/MorphyWeights.md),
[`SingleCharMorphy()`](https://ms609.github.io/TreeSearch/dev/reference/SingleCharMorphy.md),
[`UnloadMorphy()`](https://ms609.github.io/TreeSearch/dev/reference/UnloadMorphy.md),
[`is.morphyPtr()`](https://ms609.github.io/TreeSearch/dev/reference/is.morphyPtr.md),
[`mpl_apply_tipdata()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_apply_tipdata.md),
[`mpl_attach_rawdata()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_attach_rawdata.md),
[`mpl_attach_symbols()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_attach_symbols.md),
[`mpl_delete_Morphy()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_delete_Morphy.md),
[`mpl_first_down_recon()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_first_down_recon.md),
[`mpl_first_up_recon()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_first_up_recon.md),
[`mpl_get_charac_weight()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_get_charac_weight.md),
[`mpl_get_gaphandl()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_get_gaphandl.md),
[`mpl_get_num_charac()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_get_num_charac.md),
[`mpl_get_num_internal_nodes()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_get_num_internal_nodes.md),
[`mpl_get_numtaxa()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_get_numtaxa.md),
[`mpl_get_symbols()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_get_symbols.md),
[`mpl_init_Morphy()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_init_Morphy.md),
[`mpl_new_Morphy()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_new_Morphy.md),
[`mpl_second_down_recon()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_second_down_recon.md),
[`mpl_second_up_recon()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_second_up_recon.md),
[`mpl_set_charac_weight()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_set_charac_weight.md),
[`mpl_set_num_internal_nodes()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_set_num_internal_nodes.md),
[`mpl_set_parsim_t()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_set_parsim_t.md),
[`mpl_translate_error()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_translate_error.md),
[`mpl_update_lower_root()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_update_lower_root.md),
[`mpl_update_tip()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_update_tip.md),
[`summary.morphyPtr()`](https://ms609.github.io/TreeSearch/dev/reference/summary.morphyPtr.md)

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
data("Lobo", package="TreeTools")
morphyObj <- PhyDat2Morphy(Lobo.phy)
# Set object to be destroyed at end of session or closure of function
# on.exit(morphyObj <- UnloadMorphy(morphyObj), add = TRUE)

# Do something with pointer
# ....

# Or, instead of on.exit, manually destroy morphy object and free memory:
morphyObj <- UnloadMorphy(morphyObj)
```
