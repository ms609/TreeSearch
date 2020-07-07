## Test environments
* local Windows 10 install, R 4.0.2
* Ubuntu 16.04.6 LTS, R release and devel, via [Travis CI](https://travis-ci.org/ms609/TreeSearch)
* Mac OS X 10.13.6, R release, via Travis

Because 'TreeTools' 1.1.0 is not yet available on the rhub servers, it has not
been possible to test with

* win-builder, with `check_win_devel()`, R devel
* R-hub, with `rhub::check_for_cran()`

## R CMD check results
There were no ERRORs or WARNINGs.

There was one NOTE:

> Found the following (possibly) invalid URLs:
>   URL: http://doi.org/10.1093/sysbio/syy083
>     From: man/CharacterLength.Rd
>           man/inapplicable.citations.Rd
>           man/inapplicable.datasets.Rd
>           man/inapplicable.phyData.Rd
>     Status: Error
>     Message: libcurl error code 56:
>       	OpenSSL SSL_read: SSL_ERROR_SYSCALL, errno 104
> 
> Found the following (possibly) invalid DOIs:
>   DOI: 10.1093/sysbio/syy083
>     From: DESCRIPTION
>     Status: libcurl error code 56:
>     	OpenSSL SSL_read: SSL_ERROR_SYSCALL, errno 104
>     Message: Error
>   DOI: 10.2307/2412182
>     From: DESCRIPTION
>     Status: libcurl error code 56:
>     	OpenSSL SSL_read: SSL_ERROR_SYSCALL, errno 104
>     Message: Error

I have verified that the DOIs and URLs are correct.

## Downstream dependencies
There is currently one downstream dependency for this package.

`revdepcheck::revdep_check()` identified no changes to worse in:

v CongreveLamsdell2016 1.0.2             -- E: 0     | W: 0     | N: 0    
OK: 1