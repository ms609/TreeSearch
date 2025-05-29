#' Matrix reconciliation
#' @param file1,file2 Paths to matrices to compare
#' @param matchTaxa named character vector; each named taxon in `file1` will be
#'  treated as synonymous with the corresponding taxon in `file2`.
#' @importFrom TreeTools ReadCharacters ReadNotes
#' @importFrom TreeDist LAPJV
#' @importFrom stringdist stringdistmatrix
MatRec <- function(
    file1 = "C:/Users/pjjg18/GitHub/crico-phylo/TreeSearch/mbank_X30710_2025-05-27-1307.nex",
    file2 = "C:/Users/pjjg18/GitHub/larva/phylo/TreeSearch/mbank_X27231_2024-03-22-1117.nex",
    #matchTaxa = character(0)
    matchTaxa = c('Tubiluchus_lemburgi' = 'Tubiluchus', 'Cricocosmia_n._sp.' = 'Cricocosmia', 'Actinarctus_neretinus' = 'Actinarctus', 'Halobiotus_crispae' = 'Halobiotus')
) {
  mat1 <- ReadCharacters(file1)
  states1 <- attr(mat1, "state.labels")
  tax1 <- rownames(mat1)
  unmatched <- !names(matchTaxa) %in% tax1
  if (any(unmatched)) {
    stop("The following taxa in `matchTaxa` are not present in `file1`: ",
         paste(names(matchTaxa)[unmatched], collapse = ", "))
  }
  
  mat2 <- ReadCharacters(file2)
  states2 <- attr(mat2, "state.labels")
  tax2 <- rownames(mat2)
  unmatched <- !matchTaxa %in% tax2
  if (any(unmatched)) {
    stop("The following taxa in `matchTaxa` are not present in `file2`: ",
         paste(matchTaxa[unmatched], collapse = ", "))
  }
  tax2[match(matchTaxa, tax2)] <- names(matchTaxa)
  rownames(mat2) <- tax2
  
  taxAnotB <- setdiff(tax1, tax2)
  taxBnotA <- setdiff(tax2, tax1)
  
  if (length(taxAnotB) && length(taxBnotA)) {
    lev <- stringdistmatrix(taxAnotB, taxBnotA, method = "jw")
    taxMatch <- LAPJV(lev)[["matching"]]
    
    message("Found unmatched taxa; possible matching: \n matchTaxa = c(",
            paste(paste0("'", taxAnotB[!is.na(taxMatch)], "' = '",
                         taxBnotA[taxMatch[!is.na(taxMatch)]], "'"),
                  collapse = ", "),
            ")")
    if (length(taxAnotB) > length(taxBnotA)) {
      message(" leaving unmatched: \n",
              paste(sort(taxAnotB[is.na(taxMatch)]), collapse = ", "))
    } else if (length(taxAnotB) > length(taxBnotA)) {
      message(" leaving unmatched: \n",
              paste(sort(taxBnotA[-taxMatch]), collapse = ", "))
    }
  }
  
  commonTaxa <- intersect(tax1, tax2)
  mat1 <- mat1[commonTaxa, , drop = FALSE]
  mat2 <- mat2[commonTaxa, , drop = FALSE]
  
  # Identical names = identical characters
  names1 <- colnames(mat1)
  names2 <- colnames(mat2)
  charMatch <- match(names1, names2)
  
  notes1 <- ReadNotes(file1)
  notes2 <- ReadNotes(file2)
  charNote1 <- vapply(notes1, function(x) paste0(x[[1]], ""), character(1))
  charNote2 <- vapply(notes2, function(x) paste0(x[[1]], ""), character(1))
  
  .Unmatched1 <- function() {
    is.na(charMatch)
  }
  .Unmatched2 <- function() {
    !seq_along(names2) %in% charMatch[!is.na(charMatch)]
  }
  .Index2 <- function() {
    seq_along(names2)[.Unmatched2()]
  }
  # Match characters with identical notes
  charMatch[.Unmatched1()] <- .Index2()[match(charNote1[.Unmatched1()], charNote2[.Unmatched2()])]
  
  
  nameDist <- stringdistmatrix(names1[.Unmatched1()], names2[.Unmatched2()],
                              method = "jw")
  noteDist <- stringdistmatrix(charNote1[.Unmatched1()],
                               charNote2[.Unmatched2()],
                               method ="jw")
  
  stateDist <- stringdistmatrix(
    method = "jw",
    vapply(states1[.Unmatched1()], paste, collapse = " | ", character(1)),
    vapply(states2[.Unmatched2()], paste, collapse = " | ", character(1))
  )
  
  clusters <- unique(c("?", as.character(mat1[, .Unmatched1()]),
                       as.character(mat2[, .Unmatched2()])))
  
  tokenSim <- apply(matrix(match(mat1[, .Unmatched1()], clusters) - 1,
                           nrow = nrow(mat1), ncol = sum(.Unmatched1())), 2,
                    function(col1) {
    apply(matrix(match(mat2[, .Unmatched2()], clusters) - 1,
                 nrow = nrow(mat2), ncol = sum(.Unmatched2())), 2,
          .TokenCompare, col1)
                    })
  
  # If tokens are similar, prefer the match by tokenWeight
  tokenWeighting <- 2
  tokenSim[is.nan(tokenSim) | is.na(tokenSim)] <- 0
  tokenWeight <- tokenWeighting * tokenSim
  tokenWeight[tokenSim < 0.85] <- 1
  
  ret <- charMatch
  mode(ret) <- "list"
  ret[.Unmatched1()] <- apply(nameDist * noteDist * stateDist / t(tokenWeight),
                              1, order, simplify = FALSE)
  ret
}


#' Compare the mutual information of two columns
.TokenCompare <- function(a, b) {
  ambig <- a == 0 | b == 0
  oA <- a[!ambig]
  hA <- .Entropy(tabulate(oA))
  oB <- b[!ambig]
  hB <- .Entropy(tabulate(oB))
  hAB <- .Entropy(table(oA, oB))
  2 * (hA + hB - hAB) / (hA + hB)
}

#' @importFrom TreeTools ReadAsPhyDat ReadNotes
#' @importFrom ape read.tree read.nexus
#' @importFrom utils menu
ViewRec <- function(file1, file2, tree, matchTaxa,
                    matching = MatRec(file1, file2, matchTaxa),
                    startAt = 1
                    ) {
  
  pd1 <- ReadAsPhyDat(file1)
  pd2 <- ReadAsPhyDat(file2)
  ch1 <- gsub(",|\\(|\\)", "", ReadCharacters(file1))
  ch2 <- gsub(",|\\(|\\)", "", ReadCharacters(file2))
  note1 <- ReadNotes(file1)
  note2 <- ReadNotes(file2)
  
  tax2 <- names(pd2)
  tax2[match(matchTaxa, tax2)] <- names(matchTaxa)
  names(pd2) <- tax2
  rownames(ch2) <- tax2
  
  commonTaxa <- intersect(names(pd1), names(pd2))
  if (is.character(tree)) {
    tree <- tryCatch(read.tree(tree, keep.multi = TRUE)[[1]], error = function(e) {
      read.nexus(tree, force.multi = TRUE)[[1]]})
  }
  newTreeLabels <- TipLabels(tree)
  newTreeLabels[match(matchTaxa, newTreeLabels)] <- names(matchTaxa)
  tree[["tip.label"]] <- newTreeLabels
  tree <- KeepTip(tree, commonTaxa)
  treeLabels <- TipLabels(tree)
  
  par(mfrow = c(1, 2), cex = 0.7, mar = rep(0, 4))
  i <- startAt
  option <- 1
  chosen <- `length<-`(integer(0), length(matching))
  while (i <= length(matching)) {
    j <- matching[[i]][[option]]
    pad <- paste0("ch %", ceiling(log10(max(ncol(ch1), ncol(ch2)))), "s: ")
    message("Matching ", sprintf(pad, i), colnames(ch1)[[i]], "\n",
            "- with - ", sprintf(pad, j), colnames(ch2)[[j]])
    
    sameState <- ch1[treeLabels, i] == ch2[treeLabels, j]
    
    if (any(!sameState)) {
      message("States differ for ", sum(!sameState),
              if(sum(!sameState) > 1) " taxa" else " taxon")
      notesI <- substr(note1[[i]][[2]][treeLabels[!sameState]], 1, 24)
      notesJ <- substr(note2[[j]][[2]][treeLabels[!sameState]], 1, 24)
      notesI[is.na(notesI)] <- ""
      notesJ[is.na(notesJ)] <- ""
      
      message(paste("\n  ", 
                    treeLabels[!sameState], "=", 
                    ch1[treeLabels[!sameState], i], notesI, "|",
                    ch2[treeLabels[!sameState], j], notesJ))
      PlotCharacter(tree, pd1, i, tip.col = 2 - sameState)
      PlotCharacter(tree, pd2, j, direction = "l", tip.col = 2 - sameState)
    } else {
      plot.new()
      text(1, 0.5, paste("Character codings are identical", 
                         switch(i %% 4 + 1, "\U1F44D", "\U263A", "\U1F389",
                                "\U1F973")),
           xpd = NA, col = 3, cex = 2.5, font = 2)
      plot.new()
    }
    optionsToShow <- min(length(matching[[i]]), 6)
    lastOption <- option
    option <- menu(c(colnames(ch2)[head(matching[[i]], optionsToShow)],
                     "[Show all options]"),
                   title = "Match a differrent character (0 for next character):",
                   graphics = FALSE)
    if (option == optionsToShow + 1) {
      matching[[i]] <- c(matching[[i]], setdiff(seq_along(colnames(ch2)),
                                                 matching[[i]]))
      option <- match(menu(colnames(ch2),
                           title = "Match a differrent character (0 for next character):"),
                      matching[[i]])
    } else {
      chosen[[i]] <- matching[[i]][[lastOption]]
    }
    
    if (is.na(option) || option == 0) {
      i <- i + 1
      option <- 1
    }
  }
  # Return:
  chosen
}
