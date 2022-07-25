## Test environments

* local Windows 10 install, R devel

* [Github Actions](https://github.com/ms609/TreeSearch/actions):
  - Ubuntu 20.04
    - R 4.1
    - R release (tests, examples & vignettes run with valgrind & ASan)
    - R devel
  - Mac OS X 10.15.7, R release
  - Microsoft Windows Server 2019 10.0.17763, R release
  
* R-hub, with `rhub::check_for_cran()` and `devtools::check_win_devel()`


## R CMD check results
There were no ERRORs or WARNINGs or (relevant) NOTEs.


## Downstream dependencies

Reverse dependencies have been checked using "revdepcheck" on
[GitHub Actions](https://github.com/ms609/TreeSearch/actions/workflows/revdep.yml).

There are currently two downstream dependencies for this package
(and I maintain both).
