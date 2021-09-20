## Test environments

* local Windows 10 install, R 4.1.1

* [Github Actions](https://github.com/ms609/TreeSearch/actions):
  - Ubuntu 20.04
    - R 3.6.3
    - R release (tests, examples & vignettes run with valgrind)
    - R devel
  - Mac OS X 10.15.7, R release
  - Microsoft Windows Server 2019 10.0.17763, R release
  
* R-hub, with `rhub::check_for_cran()` and `devtools::check_win_devel()`


## R CMD check results
There were no ERRORs or WARNINGs.

There was one NOTE:

> Found the following (possibly) invalid URLs:
>   URL: http://www.ncbi.nlm.nih.gov/pubmed/12116579 (moved to  https://pubmed.ncbi.nlm.nih.gov/12116579/)
>     From: man/TreeLength.Rd
>     Status: 200
>     Message: OK

This URL is generated from the pmid field of a bibtex entry (by Rdpack?) so
cannot be edited manually.

## Downstream dependencies

Reverse dependencies have been checked using "revdepcheck" on
[GitHub Actions](https://github.com/ms609/TreeSearch/actions/workflows/revdep.yml).

There are currently two downstream dependencies for this package
(and I maintain both).
