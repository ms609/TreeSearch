# Cluster similar strings

Calculate string similarity using the Levenshtein distance and return
clusters of similar strings.

## Usage

``` r
ClusterStrings(x, maxCluster = 12)
```

## Arguments

- x:

  Character vector.

- maxCluster:

  Integer specifying maximum number of clusters to consider.

## Value

`NameClusters()` returns an integer assigning each element of `x` to a
cluster, with an attribute `med` specifying the median string in each
cluster, and `silhouette` reporting the silhouette coefficient of the
optimal clustering. Coefficients \< 0.5 indicate weak structure, and no
clusters are returned. If the number of unique elements of `x` is less
than `maxCluster`, all occurrences of each entry are assigned to an
individual cluster.

## See also

Other utility functions:
[`QACol()`](https://ms609.github.io/TreeSearch/dev/reference/QACol.md),
[`QuartetResolution()`](https://ms609.github.io/TreeSearch/dev/reference/QuartetResolution.md),
[`WhenFirstHit()`](https://ms609.github.io/TreeSearch/dev/reference/WhenFirstHit.md)

## Author

[Martin R. Smith](https://smithlabdurham.github.io/)
(<martin.smith@durham.ac.uk>)

## Examples

``` r
ClusterStrings(c(paste0("FirstCluster ", 1:5),
                 paste0("SecondCluster.", 8:12),
                 paste0("AnotherCluster_", letters[1:6])))
#>  [1] 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3
#> attr(,"silhouette")
#> [1] 0.911867
#> attr(,"med")
#> [1] "FirstCluster 1"   "SecondCluster.10" "AnotherCluster_a"
```
