#' @param concavity Determines the degree to which extra steps beyond the first
#' are penalized.  Specify a numeric value to use implied weighting
#' (Goloboff 1993, where `concavity` specifies _k_ in _k_ / _e_ + _k_).
#'  in implied weighting. A value of 10 is recommended;
#' TNT sets a default of 3, but this is too low in some circumstances
#' (Goloboff _et al_. 2018, Smith, 2019).
#' Better still explore the sensitivity of results under a range of
#' concavity values, e.g. `k = 2 ^ (1:7)`.
#' Specify `Inf` to weight each additional step equally.
#' Specify `'profile'` to employ profile parsimony (Faith & Trueman 2001).
