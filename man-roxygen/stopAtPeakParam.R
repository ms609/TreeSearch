#' @param stopAtPeak Logical specifying whether to terminate search once a 
#' subsequent iteration recovers a sub-optimal score.
#' Will be overridden if a passed function has an attribute `stopAtPeak` set by 
#' `attr(FunctionName, 'stopAtPeak') <- TRUE`.
