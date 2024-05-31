.DivideBetween <- function(x, n) {
  if (length(n) && n < 2) {
    x
  } else {
    rep(x %/% n, n) + tabulate(seq_len(x %% n), n)
  }
}

.MaxEntropy <- function(x) {
  if (is.character(x)) {
    x <- unname(table(x))
  }
  zones <- length(x)
  if (zones < 2) {
    0
  } else {
    x <- sort.int(x, decreasing = TRUE)
    regions <- lapply(ceiling(x / 2), integer)
    nRegions <- integer(zones)
    
    maxRegions <- floor(x / 2)
    stragglers <- x %% 2
    
    for (i in rev(seq_len(zones)[-1])) {
      newRegions <- maxRegions[[i]]
      nRegions[seq_len(i)] <- newRegions
      maxRegions <- maxRegions - newRegions
      if (stragglers[[i]]) {
        nRegions[[i]] <- nRegions[[i]] + 1
      }
    }
    if (maxRegions[[1]] > 0) {
      nRegions[[1]] <- nRegions[[1]] + 1
    }
    
    regions <- lapply(seq_len(zones),
                      function(z) .DivideBetween(x[[z]], nRegions[[z]]))
    # Return:
    #"attr<-"(.Entropy(unlist(regions)), "regions", regions)
    if (isTRUE(getOption("TODO"))) message(paste(unlist(regions), collapse = " "))
    .Entropy(unlist(regions))
  }
}
