# Reconstructs the first (downpass) nodal reconstructions

Reconstructs the preliminary nodal set for all characters for a
particular node. This function is called over a postorder sequence of
internal nodes where left and right descendants are known. Because this
function needs to be fairly high-performance, it does not do much
checking for parameter validity, thus unsafe usage of this function
might not be caught. It is up to calling functions to ensure that the
appropriate parameters have been set before use.

## Usage

``` r
mpl_first_down_recon(node_id, left_id, right_id, morphyobj)
```

## Arguments

- node_id:

  The index of the node being reconstructed.

- left_id:

  The index of the left descendant.

- right_id:

  The index of the right descendant.

- morphyobj:

  An instance of the Morphy object.

## Value

The integral parsimony length (right now)

## See also

Other Morphy API functions:
[`GapHandler()`](https://ms609.github.io/TreeSearch/dev/reference/GapHandler.md),
[`MorphyErrorCheck()`](https://ms609.github.io/TreeSearch/dev/reference/MorphyErrorCheck.md),
[`MorphyWeights()`](https://ms609.github.io/TreeSearch/dev/reference/MorphyWeights.md),
[`PhyDat2Morphy()`](https://ms609.github.io/TreeSearch/dev/reference/PhyDat2Morphy.md),
[`SingleCharMorphy()`](https://ms609.github.io/TreeSearch/dev/reference/SingleCharMorphy.md),
[`UnloadMorphy()`](https://ms609.github.io/TreeSearch/dev/reference/UnloadMorphy.md),
[`is.morphyPtr()`](https://ms609.github.io/TreeSearch/dev/reference/is.morphyPtr.md),
[`mpl_apply_tipdata()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_apply_tipdata.md),
[`mpl_attach_rawdata()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_attach_rawdata.md),
[`mpl_attach_symbols()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_attach_symbols.md),
[`mpl_delete_Morphy()`](https://ms609.github.io/TreeSearch/dev/reference/mpl_delete_Morphy.md),
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

Martin Brazeau
