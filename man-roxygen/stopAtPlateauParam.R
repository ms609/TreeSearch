#' @param stopAtPlateau Integer. If > 0, tree search will terminate if the score
#' has not improved after `stopAtPlateau` iterations.
#' Will be overridden if a passed function has an attribute `stopAtPlateau` set
#' by `attr(FunctionName, "stopAtPlateau") <- TRUE`.
