## Test environments

* local Windows 10 install, R 4.0.3

* [Github Actions](https://github.com/ms609/TreeSearch/actions):
  - Ubuntu 20.04 LTS, R 3.6.3, release and devel
  - Mac OS X 10.15.7, R release
  - Microsoft Windows Server 2019 10.0.17763, R release
  
* R-hub, with `rhub::check((platform = platforms()$name)`


## R CMD check results
There were no ERRORs, WARNINGs or NOTEs.



## Downstream dependencies
There are currently two downstream dependencies for this package (and I maintain both).

`revdepcheck::revdep_check()` identified no changes to worse in:

v CongreveLamsdell2016 1.0.2             -- E: 0     | W: 0     | N: 0    
OK: 1