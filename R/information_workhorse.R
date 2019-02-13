# Interfaces to this function are in information_interface.R

# Static functions
# Note that memoising often introduces a substantial overhead
#' @importFrom memoise memoise memoize
Btwn <- (function(i) if (max(i) == min(i)) In(i[1]) else paste0('Btwn ', min(i), max(i)))
In <- (function(i) paste0('In ', i))
TDm <- (function(i, others) {
  switch(length(unique(c(i, others))),
         In(i),
         if (match(i, others, nomatch=0)) {
           In(i)
         } else {
           Btwn(c(i, others[1]))
         },
         paste0('TDm ', i, '.', min(others), max(others)))
})
DDm <- (function(i, ii, farCols=NULL) {
  if (i == ii) {
    return (In(i))
  } else {
    if (!is.null(farCols)) {
      if (farCols[1] == farCols[2]) {
        return (TDm(i, c(farCols[1], ii)))
      } else if (i %in% farCols) {
        return (In(i))
      } else if (ii %in% farCols) {
        return(Btwn(c(i, ii)))
      }
    }
    return(paste0('DDm ', i, '.', ii))
  }
})
DDmDDm <- (function (i, ii) {
  if (i[1] == i[2]) {
    TDm(i[1], ii)
  } else if (ii[1] == ii[2]) {
    TDm(ii[1], i)
  } else if (length(missing <- which(!1:4 %in% c(i,ii))) > 0) {
    tab <- table(c(i,ii))
    In(names(tab[tab > 1]))
  } else {
    'DDmDDm'
  }
})
#Cell <- memoise(function(colno, rowno, edges) {ncol(edges) * (rowno - 1) + colno})
Cell <- function(colno, rowno, edges) ncol(edges) * (rowno - 1L) + colno
ColNo <- function(colnames) match(colnames, edgeNames)
MeetEdges <- function (edges) meetEdges[min(edges), max(edges)]    # edges as the index of the edge type in edgeNames

{
  # Constants and functions for function cw2... ()
  edgeNames <- c(paste("In", 1:4), paste("Btwn", c(12, 13, 14, 23, 24, 34)), paste('TDm', c(1.23, 1.24, 1.34, 2.13, 2.14, 2.34, 3.12, 3.14, 3.24, 4.12, 4.13, 4.23)), paste('DDm', c(1.2, 1.3, 1.4, 2.1, 2.3, 2.4, 3.1, 3.2, 3.4, 4.1, 4.2, 4.3)), 'DDmDDm')
  colourOn <- boundaryOf <- matrix(as.logical(c(
    0,0,0,0, 1,1,1,0,0,0, 1,1,1, 1,1,0, 1,1,0, 1,1,0,  1,1,1, 1,0,0, 1,0,0, 1,0,0, 0,
    0,0,0,0, 1,0,0,1,1,0, 1,1,0, 1,1,1, 1,0,1, 1,0,1,  1,0,0, 1,1,1, 0,1,0, 0,1,0, 0,
    0,0,0,0, 0,1,0,1,0,1, 1,0,1, 1,0,1, 1,1,1, 0,1,1,  0,1,0, 0,1,0, 1,1,1, 0,0,1, 0,
    0,0,0,0, 0,0,1,0,1,1, 0,1,1, 0,1,1, 0,1,1, 1,1,1,  0,0,1, 0,0,1, 0,0,1, 1,1,1, 0
  )), nrow=4, byrow=TRUE)
  memberOf <- matrix(as.logical(c(
    1,0,0,0, 1,1,1,0,0,0, 1,1,1, 0,0,0, 0,0,0, 0,0,0,  1,1,1, 0,0,0, 0,0,0, 0,0,0, 0,
    0,1,0,0, 1,0,0,1,1,0, 0,0,0, 1,1,1, 0,0,0, 0,0,0,  0,0,0, 1,1,1, 0,0,0, 0,0,0, 0,
    0,0,1,0, 0,1,0,1,0,1, 0,0,0, 0,0,0, 1,1,1, 0,0,0,  0,0,0, 0,0,0, 1,1,1, 0,0,0, 0,
    0,0,0,1, 0,0,1,0,1,1, 0,0,0, 0,0,0, 0,0,0, 1,1,1,  0,0,0, 0,0,0, 0,0,0, 1,1,1, 0
  )), nrow=4, byrow=TRUE)
  endsOf <- matrix(c(1, 1, 2, 2, 3, 3, 4, 4, 1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4, 1, 23, 1, 24, 1, 34, 2, 13, 2, 14, 2, 34, 3, 12, 3, 14, 3, 24, 4, 12, 4, 13, 4, 23, 1, 20, 1, 30, 1, 40, 2, 10, 2, 30, 2, 40, 3, 10, 3, 20, 3, 40, 4, 10, 4, 20, 4, 30, 5, 5), nrow=2)
  tripleDummyEnds <- matrix(c(rep(0, (4 + 6) * 3), c(1,2,3, 1,2,4, 1,3,4, 2,1,3, 2,1,4, 2,3,4, 3,1,2, 3,1,4, 3,2,4, 4,1,2, 4,1,3, 4,2,3)), nrow=3)
  for (i in 1:4) colourOn[i, i] <- TRUE
  edgeOf <-   matrix(c( # Duplicates boundaryOf.  #TODO verify & delete, unless useful as which(boundaryOf).
    c(Btwn(c(1,2)), Btwn(c(1,3)), Btwn(c(1,4)),
      TDm(1,c(2,3)), TDm(1,c(2,4)), TDm(1,c(3,4)),
      DDm(1,2), DDm(1,3), DDm(1,4)
    ),
    c(Btwn(c(2,1)), Btwn(c(2,3)), Btwn(c(2,4)),
      TDm(2,c(1,3)), TDm(2,c(1,4)), TDm(2,c(3,4)),
      DDm(2,1), DDm(2,3), DDm(2,4)
    ),
    c(Btwn(c(3,2)), Btwn(c(3,1)), Btwn(c(3,4)),
      TDm(3,c(2,1)), TDm(3,c(2,4)), TDm(3,c(1,4)),
      DDm(3,2), DDm(3,1), DDm(3,4)
    ),
    c(Btwn(c(4,2)), Btwn(c(4,3)), Btwn(c(4,1)),
      TDm(4,c(2,3)), TDm(4,c(2,1)), TDm(4,c(3,1)),
      DDm(4,2), DDm(4,3), DDm(4,1)
    )
  ), ncol=4)
  assign('paths', matrix(c(1,2,3,4, 1,2,4,3, 1,3,2,4, 1,3,4,2,    1,4,2,3, 1,4,3,2, 2,1,3,4, 2,1,4,3,    3,1,2,4,    3,1,4,2,    4,1,2,3,    4,1,3,2), nrow=4), envir=.GlobalEnv)
  assign('listPaths', split(t(paths), seq(12)))
  assign('Ys', matrix(c(   1,2,3,4, 1,3,2,4, 1,4,2,3, 2,1,3,4,    2,3,1,4, 2,4,1,3, 3,1,2,4, 3,2,1,4,    3,4,1,2,    4,1,2,3,    4,2,1,3,    4,3,1,2), nrow=4), envir=.GlobalEnv)
  assign('list.Ys', split(t(Ys), seq(12)))
  EdgeTypes <- function (cases) rep(cases, times=c(4, 6, 12, 12, 1))
  assign('edge.types', EdgeTypes(1:5), envir=.GlobalEnv)
  mask4 <- matrix(0, nrow=31, ncol=35)
  colnames(mask4) <- edgeNames
  rownames(mask4) <- c(paste0('path', apply(paths, 2, paste, collapse='')), paste0('y', apply(Ys, 2, paste, collapse='')),
                       paste0('hub', 1:4), 'ddm2', 'ddm3', 'ddm4')
  mask4[1 , c(Btwn(paths[1:2, 1 ]), Btwn(paths[2:3, 1 ]), Btwn(paths[3:4, 1 ]))] <- 1
  mask4[2 , c(Btwn(paths[1:2, 2 ]), Btwn(paths[2:3, 2 ]), Btwn(paths[3:4, 2 ]))] <- 1
  mask4[3 , c(Btwn(paths[1:2, 3 ]), Btwn(paths[2:3, 3 ]), Btwn(paths[3:4, 3 ]))] <- 1
  mask4[4 , c(Btwn(paths[1:2, 4 ]), Btwn(paths[2:3, 4 ]), Btwn(paths[3:4, 4 ]))] <- 1
  mask4[5 , c(Btwn(paths[1:2, 5 ]), Btwn(paths[2:3, 5 ]), Btwn(paths[3:4, 5 ]))] <- 1
  mask4[6 , c(Btwn(paths[1:2, 6 ]), Btwn(paths[2:3, 6 ]), Btwn(paths[3:4, 6 ]))] <- 1
  mask4[7 , c(Btwn(paths[1:2, 7 ]), Btwn(paths[2:3, 7 ]), Btwn(paths[3:4, 7 ]))] <- 1
  mask4[8 , c(Btwn(paths[1:2, 8 ]), Btwn(paths[2:3, 8 ]), Btwn(paths[3:4, 8 ]))] <- 1
  mask4[9 , c(Btwn(paths[1:2, 9 ]), Btwn(paths[2:3, 9 ]), Btwn(paths[3:4, 9 ]))] <- 1
  mask4[10, c(Btwn(paths[1:2, 10]), Btwn(paths[2:3, 10]), Btwn(paths[3:4, 10]))] <- 1
  mask4[11, c(Btwn(paths[1:2, 11]), Btwn(paths[2:3, 11]), Btwn(paths[3:4, 11]))] <- 1
  mask4[12, c(Btwn(paths[1:2, 12]), Btwn(paths[2:3, 12]), Btwn(paths[3:4, 12]))] <- 1
  mask4[13, c(Btwn(c(1, 2)), TDm(2, c(3, 4)), TDm(3, c(2, 4)), TDm(4, c(3, 2)))] <- 1
  mask4[14, c(Btwn(c(1, 3)), TDm(3, c(2, 4)), TDm(2, c(3, 4)), TDm(4, c(2, 3)))] <- 1
  mask4[15, c(Btwn(c(1, 4)), TDm(4, c(2, 3)), TDm(2, c(4, 3)), TDm(3, c(2, 4)))] <- 1
  mask4[16, c(Btwn(c(2, 1)), TDm(1, c(3, 4)), TDm(3, c(1, 4)), TDm(4, c(3, 1)))] <- 1
  mask4[17, c(Btwn(c(2, 3)), TDm(3, c(1, 4)), TDm(1, c(3, 4)), TDm(4, c(1, 3)))] <- 1
  mask4[18, c(Btwn(c(2, 4)), TDm(4, c(1, 3)), TDm(1, c(4, 3)), TDm(3, c(1, 4)))] <- 1
  mask4[19, c(Btwn(c(3, 1)), TDm(1, c(2, 4)), TDm(2, c(1, 4)), TDm(4, c(2, 1)))] <- 1
  mask4[20, c(Btwn(c(3, 2)), TDm(2, c(1, 4)), TDm(1, c(2, 4)), TDm(4, c(1, 2)))] <- 1
  mask4[21, c(Btwn(c(3, 4)), TDm(4, c(1, 2)), TDm(1, c(4, 2)), TDm(2, c(1, 4)))] <- 1
  mask4[22, c(Btwn(c(4, 1)), TDm(1, c(2, 3)), TDm(2, c(1, 3)), TDm(3, c(2, 1)))] <- 1
  mask4[23, c(Btwn(c(4, 2)), TDm(2, c(1, 3)), TDm(1, c(2, 3)), TDm(3, c(1, 2)))] <- 1
  mask4[24, c(Btwn(c(4, 3)), TDm(3, c(1, 2)), TDm(1, c(3, 2)), TDm(2, c(1, 3)))] <- 1
  mask4[25, c(Btwn(c(1, 1)), Btwn(c(1, 2)), Btwn(c(1, 3)), Btwn(c(1, 4)))] <- 1
  mask4[26, c(Btwn(c(2, 1)), Btwn(c(2, 2)), Btwn(c(2, 3)), Btwn(c(2, 4)))] <- 1
  mask4[27, c(Btwn(c(3, 1)), Btwn(c(3, 2)), Btwn(c(3, 3)), Btwn(c(3, 4)))] <- 1
  mask4[28, c(Btwn(c(4, 1)), Btwn(c(4, 2)), Btwn(c(4, 3)), Btwn(c(4, 4)))] <- 1
  mask4[29, c(DDm(1, 2), DDm(2, 1), DDm(3, 4), DDm(4, 3), 'DDmDDm')] <- 1
  mask4[30, c(DDm(1, 3), DDm(3, 1), DDm(2, 4), DDm(4, 2), 'DDmDDm')] <- 1
  mask4[31, c(DDm(1, 4), DDm(4, 1), DDm(3, 2), DDm(2, 3), 'DDmDDm')] <- 1
  
  mask3 <- matrix(0, nrow=4, ncol=35)
  colnames(mask3) <- edgeNames
  rownames(mask3) <- c(paste0('path', apply(matrix(c(2,1,3, 1,2,3, 1,3,2), 3, 3), 2, paste, collapse='')), 'tri')
  mask3[1 , c(Btwn(c(1, 2)), Btwn(c(1, 3)))] <- 1
  mask3[2 , c(Btwn(c(1, 2)), Btwn(c(2, 3)))] <- 1
  mask3[3 , c(Btwn(c(1, 3)), Btwn(c(2, 3)))] <- 1
  mask3[4 , c(TDm(1, c(2, 3)), TDm(2, c(1, 3)), TDm(3, c(2, 1)))] <- 1
  
  
  
  meetEdges <- matrix(0, ColNo(Btwn(3:4)), length(edgeNames))
  QDm <- memoise(function(i, i1, i2, i3) {
    if (i %in% c(i1, i2, i3)) In(i) else 'DDmDDm'
  })
  Btwn2 <- memoise(function(i, j) {
    if (all(i %in% j)) {
      Btwn(i) # Illegal configuration!
    } else if (any(i %in% j)) {
      In(i[i %in% j])
    } else {
      'DDmDDm'
    }
  })
  Overlap <- function (i, j) In(i[i %in% j][1]) # Will not fail given incompatible i & j
  for (i in 1:4) meetEdges[i, ] <- ColNo(c(
    Btwn(c(i, 1)), Btwn(c(i, 2)), Btwn(c(i, 3)), Btwn(c(i, 4)),
    TDm (i, c(1,2)), TDm (i, c(1,3)), TDm (i, c(1,4)), TDm (i, c(2,3)), TDm (i, c(2,4)), TDm (i, c(3,4)),
    QDm(i, 1, 2, 3), QDm(i, 1, 2, 4), QDm(i, 1, 3, 4),
    QDm(i, 2, 1, 3), QDm(i, 2, 1, 4), QDm(i, 2, 3, 4),
    QDm(i, 3, 1, 2), QDm(i, 3, 1, 4), QDm(i, 3, 2, 4),
    QDm(i, 4, 1, 2), QDm(i, 4, 1, 3), QDm(i, 4, 2, 3),
    rep(In(i), 13))
  )
  k <- 4
  for (i in 1:3) for (j in (1+i):4) meetEdges[k <- k + 1, 4 + 1:(6*3)] <- ColNo(c(
    Btwn2(c(i, j), c(1, 2)), Btwn2(c(i, j), c(1, 3)), Btwn2(c(i, j), c(1, 4)),
    Btwn2(c(i, j), c(2, 3)), Btwn2(c(i, j), c(2, 4)), Btwn2(c(i, j), c(3, 4)),
    Overlap(c(i, j), c(1, 2, 3)), Overlap(c(i, j), c(1, 2, 4)), Overlap(c(i, j), c(1, 3, 4)),
    Overlap(c(i, j), c(2, 1, 3)), Overlap(c(i, j), c(2, 1, 4)), Overlap(c(i, j), c(2, 3, 4)),
    Overlap(c(i, j), c(3, 1, 2)), Overlap(c(i, j), c(3, 1, 4)), Overlap(c(i, j), c(3, 2, 4)),
    Overlap(c(i, j), c(4, 1, 2)), Overlap(c(i, j), c(4, 1, 3)), Overlap(c(i, j), c(4, 2, 3))
  ))
}


# Tile editing functions
AddBoundaryToTile <- function (cols, edges) {
  # Edges is a single row corresponding to the tile in question
  if (any(soloTips <- edges[In(cols)] == 0)) {
    # not including -1, which have yet to be attached to the tree
    edges <- NULL # This would be a triple junction, and thus counted sparately
  } else {
    edges[In(cols)]   <- edges[In(cols)] + 1
    edges[Btwn(cols)] <- edges[Btwn(cols)] + 1
  }
  edges
}
AddTripleDummyTile <- function(cols, edges) {
  edges[In(cols)] <- edges[In(cols)] + 1
  for (i in 1:3) edges[TDm(cols[i], cols[-i])] <- edges[TDm(cols[i], cols[-i])] + 1
  edges
}
AddTilePaths <- function (perm, edges) {
  AddBoundaryToTile(perm[c(1,2)], AddBoundaryToTile(perm[c(2,3)], AddBoundaryToTile(perm[c(3,4)], edges)))
}
AddTileYs <- function(perm, edges) {
  AddBoundaryToTile(perm[1:2], AddTripleDummyTile(perm[2:4], edges))
}
AddTileStars <- function(hub, edges) {
  points <- (1:4)[-hub]
  AddBoundaryToTile(c(hub, points[1]), AddBoundaryToTile(c(hub, points[2]), AddBoundaryToTile(c(hub, points[3]), edges)))
}
AddTileDummyPair <- function(sisOf.4, tileEdges) {
  others <- (1:3)[-sisOf.4]
  sisters <- matrix(c(4, sisOf.4, sisOf.4, 4, others, rev(others)), nrow=2)
  for (i in 1:4) {
    tileEdges[DDm(sisters[1, i], sisters[2, i])] <-
      tileEdges[DDm(sisters[1, i], sisters[2, i])] + 1
  }
  tileEdges[In(1:4)]  <- tileEdges[In(1:4)]  + 1
  tileEdges["DDmDDm"] <- tileEdges["DDmDDm"] + 1
  tileEdges
}

ColoursOnTile <- function(tile, tileNumber) {
  if (tileNumber <= 2) tile[1:2] <- tile[1:2] + 3
  which(c(any(tile[c(1:2 %in% 1, boundaryOf[1, -2:-1])] > 0),
          any(tile[c(1:2 %in% 2, boundaryOf[2, -2:-1])] > 0)
          )
  )
}

# Edge editing functions
EditEdges <- function (edits, editingEdges) {
  # edits = vector, with every three concurrent entries corresponding to add/subtract, rows, and cols
  edits <- matrix(edits, nrow=3)
  add <- as.integer(edits[1, ])
  tiles <- edits[2, ]
  already.numeric <- tiles %in% 1:40
  tiles[!already.numeric] <- ColNo(tiles[!already.numeric])
  cols <- edits[3, ]
  already.numeric <- cols %in% 1:40
  cols[!already.numeric] <- ColNo(cols[!already.numeric])
  cells <- Cell(as.integer(tiles), as.integer(cols), editingEdges)
  if (any(editingEdges[cells][add < 0] <= 0)) return (NULL)
  for (i in 1:length(cells)) editingEdges[cells[i]] <- editingEdges[cells[i]] + add[i]
  editingEdges
}
InToIn <- function (tiles, colour, edges) {
  # In this and the following, tiles contains
  #   [1] Location of 1st branch to break
  #   [2] Location of 2nd branch to break
  #   [3] Location of left half of 1st branch to create
  #   [4] Location of right half of 1st branch to create
  #   [5] Location of left half of 2nd branch to create
  #   [6] Location of right half of 2nd branch to create
  #   [7] Location of new branch
  EditEdges(c(-1, tiles[1], In(colour),
              -1, tiles[2], In(colour),
              +1, tiles[3], In(colour),
              +1, tiles[4], In(colour),
              +1, tiles[5], In(colour),
              +1, tiles[6], In(colour),
              +1, tiles[7], In(colour)
  ), edges)
}
InToSingleton <- function (tiles, colour, edges) {
  if (sum(edges[tiles[2], 1:4]) == -3-3-3-1 && sum(edges[tiles[2], -(1:4)]) == 0 && edges[tiles[2], colour] == -1) {
    # Singleton is present
    EditEdges(c(-1, tiles[1], In(colour),
                00, tiles[2], In(colour), # No edge to break
                +1, tiles[3], In(colour),
                +1, tiles[4], In(colour),
                00, tiles[5], In(colour), # As didn't break edge above
                +1, tiles[6], In(colour),
                +1, tiles[7], In(colour)
    ), edges)
  } else {
    NULL
  }
}
In1ToIn2 <- function (tiles, col1, col2, edges) {
  EditEdges(c(-1, tiles[1], In(col1),
              -1, tiles[2], In(col2),
              +1, tiles[3], In(col1),
              +1, tiles[4], In(col1),
              +1, tiles[5], In(col2),
              +1, tiles[6], In(col2),
              +1, tiles[7], Btwn(c(col1, col2))
  ), edges)
}
InToBtwn <- function (tiles, inCol, btwnCols, edges) {
  EditEdges(c(-1, tiles[1], In(inCol),
              -1, tiles[2], Btwn(btwnCols),
              +1, tiles[3], In(inCol),
              +1, tiles[4], In(inCol),
              +1, tiles[5], TDm(btwnCols[1], c(btwnCols[2], inCol)),
              +1, tiles[6], TDm(btwnCols[2], c(btwnCols[1], inCol)),
              +1, tiles[7], TDm(inCol, btwnCols)
  ), edges)
}
InToTdm <- function (tiles, inCol, tdmAttach, tdmOther, edges) {
  if (inCol == tdmAttach) {
    EditEdges(c(-1, tiles[1], In(inCol),
                -1, tiles[2], TDm(inCol, tdmOther),
                +1, tiles[3], In(inCol),
                +1, tiles[4], In(inCol),
                +1, tiles[5], In(inCol),
                +1, tiles[6], TDm(inCol, tdmOther),
                +1, tiles[7], In(inCol)
    ), edges)
  } else if (inCol %in% tdmOther) {
    tdmOtherLoc1 <- which(edges[, TDm(tdmOther[1], c(tdmAttach, tdmOther[2]))] == 1)
    tdmOtherLoc2 <- which(edges[, TDm(tdmOther[2], c(tdmAttach, tdmOther[1]))] == 1)
    if (length(tdmOtherLoc1) == 0 || length(tdmOtherLoc2) == 0) return (NULL)
    stopifnot(length(tdmOtherLoc1) == 1)
    stopifnot(length(tdmOtherLoc2) == 1)
    EditEdges(c(-1, tiles[1], In(inCol),
                +1, tiles[3], In(inCol),
                +1, tiles[4], In(inCol),
                
                -1, tiles[2], TDm(tdmAttach, tdmOther),
                +1, tiles[5], In(inCol),
                +1, tiles[6], Btwn(c(inCol, tdmAttach)),
                
                -1, tdmOtherLoc1, TDm(tdmOther[1], c(tdmAttach, tdmOther[2])),
                +1, tdmOtherLoc1, Btwn(c(inCol, tdmOther[1])),
                
                -1, tdmOtherLoc2, TDm(tdmOther[2], c(tdmAttach, tdmOther[1])),
                +1, tdmOtherLoc2, Btwn(c(inCol, tdmOther[2])),
                
                +1, tiles[7], In(inCol)
    ), edges)
  } else {
    tdmOtherLoc1 <- which(edges[, TDm(tdmOther[1], c(tdmAttach, tdmOther[2]))] == 1)
    tdmOtherLoc2 <- which(edges[, TDm(tdmOther[2], c(tdmAttach, tdmOther[1]))] == 1)
    if (length(tdmOtherLoc1) == 0 || length(tdmOtherLoc2) == 0) return (NULL)
    stopifnot(length(tdmOtherLoc1) == 1)
    stopifnot(length(tdmOtherLoc2) == 1)
    EditEdges(c(-1, tiles[1], In(inCol),
                +1, tiles[3], In(inCol),
                +1, tiles[4], In(inCol),
                
                +1, tiles[7], DDm(inCol, tdmAttach),
                
                -1, tiles[2], TDm(tdmAttach, tdmOther),
                +1, tiles[5], DDm(tdmAttach, inCol),
                +1, tiles[6], 'DDmDDm',
                
                -1, tdmOtherLoc1, TDm(tdmOther[1], c(tdmAttach, tdmOther[2])),
                +1, tdmOtherLoc1, DDm(tdmOther[1], tdmOther[2]),
                -1, tdmOtherLoc2, TDm(tdmOther[2], c(tdmAttach, tdmOther[1])),
                +1, tdmOtherLoc2, DDm(tdmOther[2], tdmOther[1])
                
    ), edges)
  }
}
InToDdm <- function (tiles, inCol, sisCol, edges) {
  EditEdges(c(-1, tiles[1], In(inCol),
              -1, tiles[2], DDm(inCol, sisCol),
              +1, tiles[3], In(inCol),
              +1, tiles[4], In(inCol),
              +1, tiles[5], In(inCol),
              +1, tiles[6], DDm(inCol, sisCol),
              +1, tiles[7], In(inCol)
  ), edges)
}
ResolveDdm <- function (tiles, inCol, targetCol, edges) {
  if (sum(dummies <- edges[tiles[2], 0:2 + ColNo(DDm(targetCol, min((1:4)[-targetCol])))]) != 1) return (NULL)
  sisCol <- (1:4)[-targetCol][as.logical(dummies)]
  whereDummyOther.half <- which(edges[, DDm(targetCol, sisCol)] == 1)
  whereDummyEdge <- which(edges[, 'DDmDDm'] == 1)
  otherDummyEnds <- (1:4)[-c(targetCol, sisCol)]
  whereDummy.12   <- which(edges[, DDm(otherDummyEnds[1], otherDummyEnds[2])] == 1)
  whereDummy.21   <- which(edges[, DDm(otherDummyEnds[2], otherDummyEnds[1])] == 1)
  if (inCol == targetCol) {
    EditEdges(c(# New edges on easy tile
      -1, tiles[1], In(inCol),
      +1, tiles[3], In(inCol),
      +1, tiles[4], In(inCol),
      -1, tiles[2], DDm(targetCol, sisCol),
      +1, tiles[5], DDm(targetCol, sisCol),
      +1, tiles[6], In(targetCol),
      +1, tiles[7], In(targetCol)
    ), edges)
  } else {
    EditEdges(c(# New edges on easy tile
      -1, tiles[1], In(inCol),
      +1, tiles[3], In(inCol),
      +1, tiles[4], In(inCol),
      
      # Join tiles
      +1, tiles[7], In(inCol),
      
      # Delete double dummy edges and add replacements
      -1, tiles[2], DDm(targetCol, sisCol),
      +1, tiles[5], Btwn(c(inCol, targetCol)),
      +1, tiles[6], In(inCol),
      
      -1, whereDummyOther.half, DDm(sisCol, targetCol),
      +1, whereDummyOther.half, Btwn(c(inCol, sisCol)),
      -1, whereDummyEdge, 'DDmDDm',
      +1, whereDummyEdge, TDm(inCol, otherDummyEnds),
      -1, whereDummy.12, DDm(otherDummyEnds[1], otherDummyEnds[2]),
      +1, whereDummy.12, TDm(otherDummyEnds[1], c(inCol, otherDummyEnds[2])),
      -1, whereDummy.21, DDm(otherDummyEnds[2], otherDummyEnds[1]),
      +1, whereDummy.21, TDm(otherDummyEnds[2], c(inCol, otherDummyEnds[1]))
    ), edges)
  }
}
InToDdmDdm <- function (tiles, inCol, edges) {
  otherCols <- (1:4)[-inCol]
  ddm.decode <- matrix(ColNo(c(DDm(1,2), DDm(1,3), DDm(1,4),
                               DDm(2,1), DDm(2,3), DDm(2,4),
                               DDm(3,1), DDm(3,2), DDm(3,4),
                               DDm(4,1), DDm(4,2), DDm(4,3))), nrow=3)
  sisCol <- otherCols[colSums(edges[, ddm.decode[,inCol]]) > 0]
  if (!length(sisCol)) return (NULL)
  farCols <- (1:4)[-c(inCol, sisCol)]
  inLoc <- which(edges[, DDm(inCol, sisCol)] == 1)
  sisLoc <- which(edges[, DDm(sisCol, inCol)] == 1)
  farLoc1  <- which(edges[, DDm(farCols[1], farCols[2])] == 1)
  farLoc2  <- which(edges[, DDm(farCols[2], farCols[1])] == 1)
  stopifnot(length(c(inLoc, sisLoc, farLoc1, farLoc2)) == 4)
  EditEdges(c(-1, tiles[1], In(inCol),
              +1, tiles[3], In(inCol),
              +1, tiles[4], In(inCol),
              
              -1, tiles[2], 'DDmDDm',
              +1, tiles[5], In(inCol),
              +1, tiles[6], TDm(inCol, farCols),
              
              -1, inLoc, DDm(inCol, sisCol),
              +1, inLoc, In(inCol),
              -1, sisLoc, DDm(sisCol, inCol),
              +1, sisLoc, Btwn(c(sisCol, inCol)),
              -1, farLoc1, DDm(farCols[1], farCols[2]),
              +1, farLoc1, TDm(farCols[1], c(inCol, farCols[2])),
              -1, farLoc2, DDm(farCols[2], farCols[1]),
              +1, farLoc2, TDm(farCols[2], c(inCol, farCols[1])),
              +1, tiles[7], In(inCol)
  ), edges)
}
BtwnToBtwn <- function (tiles, btwn1Cols, btwn2Cols, edges) {
  EditEdges(c(-1, tiles[1], Btwn(btwn1Cols),
              -1, tiles[2], Btwn(btwn2Cols),
              +1, tiles[3], DDm(btwn1Cols[1], btwn1Cols[2], btwn2Cols),
              +1, tiles[4], DDm(btwn1Cols[2], btwn1Cols[1], btwn2Cols),
              +1, tiles[5], DDm(btwn2Cols[1], btwn2Cols[2], btwn1Cols),
              +1, tiles[6], DDm(btwn2Cols[2], btwn2Cols[1], btwn1Cols),
              +1, tiles[7], DDmDDm(btwn1Cols, btwn2Cols)
  ), edges)
}
BtwnToTdm <- function (tiles, sharedCol, btwnCol, tdmCols, edges) {
  # sharedCol is present on both tiles
  # btwnCol is the other end of the edge on tile 1
  # tdmCols are the three colours on the triple dummy, starting with the colour of the edge that we'll break
  tdmOtherLoc1 <- which(edges[, TDm(tdmCols[2], tdmCols[c(1, 3)])] == 1)
  tdmOtherLoc2 <- which(edges[, TDm(tdmCols[3], tdmCols[c(1, 2)])] == 1)
  if (length(c(tdmOtherLoc1, tdmOtherLoc2)) < 2) return(NULL)
  EditEdges(c(-1, tiles[1], Btwn(c(sharedCol, btwnCol)),
              +1, tiles[3], In(sharedCol),
              +1, tiles[4], Btwn(c(sharedCol, btwnCol)),
              
              -1, tiles[2], TDm(tdmCols[1], tdmCols[(1:3)[-1]]),
              +1, tiles[5], Btwn(c(sharedCol, tdmCols[1])),
              +1, tiles[6], TDm(sharedCol, tdmCols[2:3]),
              
              -1, tdmOtherLoc1, TDm(tdmCols[2], tdmCols[c(1, 3)]),
              +1, tdmOtherLoc1, TDm(tdmCols[2], c(sharedCol, tdmCols[3])),
              -1, tdmOtherLoc2, TDm(tdmCols[3], tdmCols[c(1, 2)]),
              +1, tdmOtherLoc2, TDm(tdmCols[3], c(sharedCol, tdmCols[2])),
              +1, tiles[7], In(sharedCol)
  ), edges)
}

#' Link two tiles
LinkTwoTiles <- function (edgesList) {
  
  if (is.null(edgesList)) return (NULL)
  if (!is.list(edgesList)) edgesList <- list(edgesList)
  # We are joining two tiles, rather than a dummy point
  oneToTwo <- ColNo(c(In(tile1), In(tile2), In(tile1), In(tile1), In(tile2), In(tile2), Btwn(c(tile1, tile2))))
  twoToOne <- oneToTwo[c(2,1,5,6,3,4,7)]
  
  output <- unique(RemoveNullsFromList(unlist(lapply(edgesList, function (edges) {
    if (is.null(edges)) return(NULL)
    tile1Colours <- ColoursOnTile(edges[1, ], 1L)
    tile2Colours <- ColoursOnTile(edges[2, ], 2L)
    if (length(tile1Colours) == 0 || length(tile2Colours) == 0) return(NULL)
    sharedColour <- tile1Colours[tile1Colours %in% tile2Colours]
    #
    #
    #
    #
    #
    #
    #
    # We are here.
    # Simplify this function to eliminate options that don't exist for 2 
    # tiles, 2 colours
    #
    #
    #
    #
    #
    #
    if (length(sharedColour)) {
      # Colour in common
      unsharedColour <- (1:2)[-sharedColour]
      newList <- c(list(InToIn(oneToTwo, sharedColour, edges)),
                   list(InToSingleton(oneToTwo, sharedColour, edges)),
                   list(
                     InToBtwn(oneToTwo, sharedColour, c(unsharedColour, sharedColour), edges),
                     InToBtwn(twoToOne, sharedColour, c(unsharedColour, sharedColour), edges)
                   ),
                   lapply(list(1:2, 2:3, c(1,3), 2:1, 3:2, c(3, 1)), function (choice)
                     BtwnToBtwn(oneToTwo, c(sharedColour, unsharedColours[choice[1]]), c(sharedColour, unsharedColours[choice[2]]), edges)
                   ), {
                     if (length(tile1Colours) > length(tile2Colours)) {
                       btwnColour <- tile2Colours[!(tile2Colours %in% sharedColour)]
                       inToDummy <- twoToOne
                       dummyColours <- tile1Colours
                     } else {
                       btwnColour <- tile1Colours[!(tile1Colours %in% sharedColour)]
                       inToDummy <- oneToTwo
                       dummyColours <- tile2Colours
                     }
                     if (length(dummyColours) == 3) {
                       # otherColours <- dummyColours[dummyColours!=sharedColour]
                       list(
                         InToTdm  (inToDummy, sharedColour, dummyColours[1], dummyColours[-1], edges),
                         InToTdm  (inToDummy, sharedColour, dummyColours[2], dummyColours[-2], edges),
                         InToTdm  (inToDummy, sharedColour, dummyColours[3], dummyColours[-3], edges),
                         BtwnToTdm(inToDummy, sharedColour, btwnColour, dummyColours[c(1, (1:3)[-1])], edges),
                         BtwnToTdm(inToDummy, sharedColour, btwnColour, dummyColours[c(2, (1:3)[-2])], edges),
                         BtwnToTdm(inToDummy, sharedColour, btwnColour, dummyColours[c(3, (1:3)[-3])], edges)
                       )
                     }
                   }, {
                     # tdmTo.tdm impossible
                     sisterRange <- ColNo(DDm(sharedColour, min(unsharedColours)))
                     sisterToShared <- unsharedColours[as.logical(colSums(edges))[sisterRange + 0:2]]
                     if (length(sisterToShared)) {
                       list(
                         ResolveDdm(oneToTwo, sharedColour, sisterToShared, edges),
                         ResolveDdm(twoToOne, sharedColour, sisterToShared, edges),
                         InToDdm(oneToTwo, sharedColour, sisterToShared, edges),
                         InToDdm(twoToOne, sharedColour, sisterToShared, edges)
                       )
                     }
                   },
                   # ddm to anything else impossible
                   list(InToDdmDdm(oneToTwo, sharedColour, edges)),
                   list(InToDdmDdm(twoToOne, sharedColour, edges))
      )
      RemoveNullsFromList(newList)
    } else {
      # No colour in common
      switch(length(tile1Colours), { # 1
        switch(length(tile2Colours), { # 1, 1
          list(In1ToIn2(oneToTwo, tile1Colours, tile2Colours, edges))
        }, { # 1, 2
          c(
            lapply(1:2, function (i) In1ToIn2(oneToTwo, tile1Colours, tile2Colours[i], edges)),
            list(InToBtwn(oneToTwo, tile1Colours, tile2Colours, edges))
          )
        }, { # 1, 3
          unlist(lapply(1:3, function (i)
            list(
              In1ToIn2(oneToTwo, tile1Colours, tile2Colours[i], edges),
              InToBtwn(oneToTwo, tile1Colours, tile2Colours[-i], edges),
              InToTdm(oneToTwo, tile1Colours, tile2Colours[i], tile2Colours[-i], edges)
            )
          ), recursive=FALSE)
        })
      }, { # 2
        switch(length(tile2Colours),
               { # 2, 1
                 c(
                   lapply(1:2, function (i) In1ToIn2(oneToTwo, tile1Colours[i], tile2Colours, edges)),
                   list(InToBtwn(twoToOne, tile2Colours, tile1Colours, edges))
                 )
               }, { # 2, 2
                 c(
                   lapply(list(c(1,1), c(1,2), c(2,1), c(2,2)), function (choice)
                     In1ToIn2(oneToTwo, tile1Colours[choice[1]], tile2Colours[choice[2]], edges)
                   ),
                   list(
                     InToBtwn(oneToTwo, tile1Colours[1], tile2Colours, edges),
                     InToBtwn(twoToOne, tile2Colours[1], tile1Colours, edges),
                     InToBtwn(oneToTwo, tile1Colours[2], tile2Colours, edges),
                     InToBtwn(twoToOne, tile2Colours[2], tile1Colours, edges)
                   )
                   ,
                   list(BtwnToBtwn(oneToTwo, tile1Colours, tile2Colours, edges)) # Previously marked as containing a bug, which I think is now resolved?
                 )
               }
        )
      }, { # 3
        unlist(lapply(1:3, function (i)
          list(
            In1ToIn2(twoToOne, tile2Colours, tile1Colours[i], edges),
            InToBtwn(twoToOne, tile2Colours, tile1Colours[-i], edges),
            InToTdm (twoToOne, tile2Colours, tile1Colours[i], tile1Colours[-i], edges)
          )
        ), recursive=FALSE)
      })
    }
  }), recursive=FALSE)))
  output
}

LinkTwoByColour <- function (tile1, tile2, colour, edgesList) {
  stopifnot(tile1 <= 4) # Designed for use with Hub configurations, thus will only join one tile to another.
  stopifnot(tile2 <= 4) # Designed for use with Hub configurations, thus will only join one tile to another.
  stopifnot(length(colour) == 1 && colour <= 4)
  if (is.null(edgesList)) return (NULL)
  if (!is.list(edgesList)) edgesList <- list(edgesList)
  oneToTwo <- ColNo(c(In(tile1), In(tile2), In(tile1), In(tile1), In(tile2), In(tile2), Btwn(c(tile1, tile2))))
  twoToOne <- oneToTwo[c(2,1,5,6,3,4,7)]
  output <- unique(RemoveNullsFromList(unlist(lapply(edgesList, function (edges) {
    if (is.null(edges)) return(NULL)
    tile1Colours <- ColoursOnTile(edges[tile1, ], tile1)
    tile2Colours <- ColoursOnTile(edges[tile2, ], tile2)
    stopifnot(length(tile1Colours) > 0 && length(tile2Colours) > 0)
    stopifnot((colourOn.1 <- colour %in% tile1Colours) || (colourOn.2 <- colour %in% tile2Colours))
    sharedColour <- tile1Colours[tile1Colours %in% tile2Colours]
    if (length(sharedColour) > 1 || (length(sharedColour) && (sharedColour != colour))) return (NULL)
    unsharedColours <- (1:4)[-colour]
    if (length(sharedColour) > 0) {
      # Colour in common; specifically, this colour must be 'colour'
      newList <- c(list(InToIn(oneToTwo, colour, edges)),
                   unlist(lapply(unsharedColours, function (unshared)
                     list(
                       InToBtwn(oneToTwo, colour, c(unshared, colour), edges),
                       InToBtwn(twoToOne, colour, c(unshared, colour), edges)
                     )
                   ), recursive=FALSE),
                   {
                     if (length(tile1Colours) > length(tile2Colours)) {
                       btwnColour   <- tile2Colours[!(tile2Colours %in% colour)]
                       inToDummy   <- twoToOne
                       dummyColours <- tile1Colours
                     } else {
                       btwnColour   <- tile1Colours[!(tile1Colours %in% colour)]
                       inToDummy   <- oneToTwo
                       dummyColours <- tile2Colours
                     }
                     if (length(dummyColours) == 3) {
                       # otherColours <- dummyColours[dummyColours!=sharedColour]
                       list(
                         InToTdm  (inToDummy, colour, dummyColours[1], dummyColours[-1], edges),
                         InToTdm  (inToDummy, colour, dummyColours[2], dummyColours[-2], edges),
                         InToTdm  (inToDummy, colour, dummyColours[3], dummyColours[-3], edges),
                         BtwnToTdm(inToDummy, colour, btwnColour, dummyColours[c(1, (1:3)[-1])], edges),
                         BtwnToTdm(inToDummy, colour, btwnColour, dummyColours[c(2, (1:3)[-2])], edges),
                         BtwnToTdm(inToDummy, colour, btwnColour, dummyColours[c(3, (1:3)[-3])], edges)
                       )
                     }
                   }, {
                     # tdmTo.tdm impossible
                     sisterRange <- ColNo(DDm(sharedColour, min(unsharedColours)))
                     sisterToShared <- unsharedColours[as.logical(colSums(edges))[sisterRange + 0:2]]
                     stopifnot(sisterToShared < 2)
                     if (length(sisterToShared)) {
                       list(
                         ResolveDdm(oneToTwo, sharedColour, sisterToShared, edges),
                         ResolveDdm(twoToOne, sharedColour, sisterToShared, edges),
                         InToDdm(oneToTwo, sharedColour, sisterToShared, edges),
                         InToDdm(twoToOne, sharedColour, sisterToShared, edges)
                       )
                     }
                   },
                   # ddm to anything else impossible
                   list(InToDdmDdm(oneToTwo, sharedColour, edges)),
                   list(InToDdmDdm(twoToOne, sharedColour, edges))
      )
      message('shared colour in TinkTwoByColour')
      RemoveNullsFromList(newList)
    } else {
      otherTile1Colours <- tile1Colours[!tile1Colours %in% colour]
      otherTile2Colours <- tile2Colours[!tile2Colours %in% colour]
      sisterRange     <- ColNo(DDm(colour, unsharedColours[1]))
      sisterToShared <- unsharedColours[as.logical(colSums(edges))[sisterRange + 0:2]]
      stopifnot(length(sisterToShared) <= 1)
      newList <- switch(length(tile1Colours), { # 1 colour on tile 1
        if (colourOn.1 && length(tile2Colours) > 1) {
          list(
            InToBtwn  (oneToTwo, colour, tile2Colours, edges),
            InToTdm   (oneToTwo, colour, colour, tile2Colours, edges),
            InToTdm   (oneToTwo, colour, tile2Colours[1], c(colour, tile2Colours[2]), edges),
            InToTdm   (oneToTwo, colour, tile2Colours[2], c(colour, tile2Colours[1]), edges),
            if (length(sisterToShared)) InToDdm   (oneToTwo, colour, sisterToShared, edges),
            if (length(sisterToShared)) ResolveDdm(oneToTwo, colour, 1, edges),
            if (length(sisterToShared)) ResolveDdm(oneToTwo, colour, 2, edges),
            if (length(sisterToShared)) ResolveDdm(oneToTwo, colour, 3, edges),
            if (length(sisterToShared)) ResolveDdm(oneToTwo, colour, 4, edges),
            InToDdmDdm(oneToTwo, colour, edges)
          )
        } else {
          # Will never work: unless you link to a boundary in the tile that lacks 'colour', you'll
          # end up making a 'Y' configuration.
        }
      }, { # 2 colours on tile 1
        switch(length(tile2Colours),
               { # 1 colour on tile 2
                 if (colourOn.1) {
                   # Link colour is one of two colours on tile one.
                   # No suitable boundary
                 } else {
                   # Link colour is only colour on tile two.
                   list(
                     InToBtwn  (twoToOne, colour, tile1Colours, edges),
                     InToTdm   (twoToOne, colour, colour, tile1Colours, edges),
                     InToTdm   (twoToOne, colour, tile1Colours[1], c(colour, tile1Colours[2]), edges),
                     InToTdm   (twoToOne, colour, tile1Colours[2], c(colour, tile1Colours[1]), edges),
                     if (length(sisterToShared)) InToDdm   (twoToOne, colour, sisterToShared, edges),
                     if (length(sisterToShared)) ResolveDdm(twoToOne, colour, 1, edges),
                     if (length(sisterToShared)) ResolveDdm(twoToOne, colour, 2, edges),
                     if (length(sisterToShared)) ResolveDdm(twoToOne, colour, 3, edges),
                     if (length(sisterToShared)) ResolveDdm(twoToOne, colour, 4, edges),
                     InToDdmDdm(twoToOne, colour, edges)
                   )
                 }
               }, { # 2 colours on tile 2
                 if (colourOn.1) {
                   list(
                     InToBtwn  (oneToTwo, colour, tile2Colours, edges),
                     BtwnToBtwn(oneToTwo, tile1Colours, tile2Colours, edges),
                     InToTdm   (oneToTwo, colour, colour, tile1Colours, edges),
                     InToTdm   (oneToTwo, colour, tile2Colours[1], c(colour, tile2Colours[2]), edges),
                     InToTdm   (oneToTwo, colour, tile2Colours[2], c(colour, tile2Colours[1]), edges),
                     BtwnToTdm (oneToTwo, colour, tile1Colours[!(tile1Colours %in% colour)], c(colour, tile2Colours), edges),
                     BtwnToTdm (oneToTwo, colour, tile1Colours[!(tile1Colours %in% colour)], c(tile2Colours[1], colour, tile2Colours[2]), edges),
                     BtwnToTdm (oneToTwo, colour, tile1Colours[!(tile1Colours %in% colour)], c(tile2Colours[2], colour, tile2Colours[1]), edges)
                   )
                 } else {
                   list(
                     InToBtwn  (twoToOne, colour, tile1Colours, edges),
                     BtwnToBtwn(twoToOne, tile2Colours, tile1Colours, edges),
                     InToTdm   (twoToOne, colour, colour, tile2Colours, edges),
                     InToTdm   (twoToOne, colour, tile1Colours[1], c(colour, tile1Colours[2]), edges),
                     InToTdm   (twoToOne, colour, tile1Colours[2], c(colour, tile1Colours[1]), edges),
                     BtwnToTdm (twoToOne, colour, tile2Colours[!(tile2Colours %in% colour)], c(colour, tile1Colours), edges),
                     BtwnToTdm (twoToOne, colour, tile2Colours[!(tile2Colours %in% colour)], c(tile1Colours[1], colour, tile1Colours[2]), edges),
                     BtwnToTdm (twoToOne, colour, tile2Colours[!(tile2Colours %in% colour)], c(tile1Colours[2], colour, tile1Colours[1]), edges)
                   )
                 }
               }
        )
      }, { # 3 colours on tile 1
        if (colourOn.1) {
          # No boundaries on tile 2
        } else {
          unlist(lapply(1:3, function (i)
            list(
              InToBtwn(twoToOne, colour, tile1Colours[-i], edges),
              InToTdm (twoToOne, colour, tile1Colours[i], tile1Colours[-i], edges)
            )
          ), recursive=FALSE)
        }
      })
      RemoveNullsFromList(newList)
    }
  }), recursive=FALSE)))
  output
}
TripleThreeColours <- function (tiles, colours, edges) {
  # Tiles must have no colours in common
  EditEdges(c(-1, tiles[1], colours[1],
              -1, tiles[2], colours[2],
              -1, tiles[3], colours[3],
              +2, tiles[1], colours[1],
              +2, tiles[2], colours[2],
              +2, tiles[3], colours[3],
              +1, TDm(tiles[1], tiles[-1]), TDm(colours[1], colours[-1]),
              +1, TDm(tiles[2], tiles[-2]), TDm(colours[2], colours[-2]),
              +1, TDm(tiles[3], tiles[-3]), TDm(colours[3], colours[-3])
  ), edges)
}
TripleFourColours <- function (twoColourTile, oneColourTiles, coloursTwo, coloursOne, edges) {
  list(EditEdges(c(-1, twoColourTile, coloursTwo[1],
                   -1, oneColourTiles[1], coloursOne[1],
                   -1, oneColourTiles[2], coloursOne[2],
                   +2, twoColourTile, coloursTwo[1],
                   +2, oneColourTiles[1], coloursOne[1],
                   +2, oneColourTiles[2], coloursOne[2],
                   
                   +1, TDm(twoColourTile, oneColourTiles), TDm(coloursTwo[1], coloursOne),
                   +1, TDm(oneColourTiles[1], c(twoColourTile, oneColourTiles[2])),
                   TDm(coloursOne[1], c(coloursTwo[1], coloursOne[2])),
                   +1, TDm(oneColourTiles[2], c(twoColourTile, oneColourTiles[1])),
                   TDm(coloursOne[2], c(coloursTwo[1], coloursOne[1]))),
                 edges),
       EditEdges(c(-1, twoColourTile, coloursTwo[2],
                   -1, oneColourTiles[1], coloursOne[1],
                   -1, oneColourTiles[2], coloursOne[2],
                   +2, twoColourTile, coloursTwo[2],
                   +2, oneColourTiles[1], coloursOne[1],
                   +2, oneColourTiles[2], coloursOne[2],
                   
                   +1, TDm(twoColourTile, oneColourTiles), TDm(coloursTwo[2], coloursOne),
                   +1, TDm(oneColourTiles[1], c(twoColourTile, oneColourTiles[2])),
                   TDm(coloursOne[1], c(coloursTwo[2], coloursOne[2])),
                   +1, TDm(oneColourTiles[2], c(twoColourTile, oneColourTiles[1])),
                   TDm(coloursOne[2], c(coloursTwo[2], coloursOne[1]))),
                 edges),
       EditEdges(c(-1, twoColourTile, Btwn(coloursTwo),
                   -1, oneColourTiles[1], coloursOne[1],
                   -1, oneColourTiles[2], coloursOne[2],
                   
                   +1, twoColourTile, DDm(coloursTwo[1], coloursTwo[2]),
                   +1, twoColourTile, DDm(coloursTwo[2], coloursTwo[1]),
                   +1, TDm(twoColourTile, oneColourTiles), 'DDmDDm',
                   
                   +2, oneColourTiles[1], coloursOne[1],
                   +1, TDm(oneColourTiles[1], c(twoColourTile, oneColourTiles[2])), DDm(coloursOne[1], coloursOne[2]),
                   
                   +2, oneColourTiles[2], coloursOne[2],
                   +1, TDm(oneColourTiles[2], c(twoColourTile, oneColourTiles[1])), DDm(coloursOne[2], coloursOne[1])),
                 edges))
}
LinkByTriple <- function (tiles, edgesList) {
  if (is.null(edges)) return (NULL)
  if (!is.list(edgesList)) edgesList <- list(edgesList)
  unlist(lapply(edgesList, function (edges) {
    coloursOnTile <- list(ColoursOnTile(edges[tiles[1], ], tiles[1]),
                          ColoursOnTile(edges[tiles[2], ], tiles[2]),
                          ColoursOnTile(edges[tiles[3], ], tiles[3]))
    nTilesWithColours <- table(unlist(coloursOnTile))
    if (sum(nTilesWithColours > 1) > 1) return (NULL) # We can't have a triple junction
    switch(max(nTilesWithColours), { # 1: No colours in common; triple junction
      nColours <- vapply(coloursOnTile, length, integer(1))
      if (max(nColours) == 1) {
        list(TripleThreeColours(tiles, unlist(coloursOnTile), edges))
      } else {
        twoColours <- (nColours == max(nColours))
        TripleFourColours(tiles[twoColours], tiles[!twoColours],
                          unlist(coloursOnTile[twoColours]),
                          unlist(coloursOnTile[!twoColours]), edges)
      }
    }, { # 2: Colour in common between two tiles; join these tiles then add a boundary to their junction
      colour.inCommon <- as.integer(names(nTilesWithColours[nTilesWithColours == max(nTilesWithColours)]))
      commonColour.present <- vapply(coloursOnTile,
                                     function (col) colour.inCommon %in% col,
                                     logical(1))
      tilesWithCommon <- tiles[commonColour.present]
      oddTile.out <- tiles[!commonColour.present]
      LinkTwoTiles(oddTile.out, ColNo(Btwn(tilesWithCommon)),
                   LinkTwoTiles(tilesWithCommon[1], tilesWithCommon[2], edges))
    }, { # 3: Colour in common to all tiles
      LinkTwoTiles(tiles[1], ColNo(Btwn(tiles[2:3])),
                   LinkTwoTiles(tiles[2], tiles[3], edges))
    })
  }), recursive = FALSE)
}
DoubleDummyLink <- function (edgesList) {
  # Returns all valid double dummies
  if (is.null(edgesList)) return (NULL)
  if (!is.list(edgesList)) edgesList <- list(edgesList)
  unlist(lapply(edgesList, function (edges) {
    coloursOnTile <- vapply(1:4, function (tile) {
      ret <- logical(4)
      ret[ColoursOnTile(edges[tile, ], tile)] <- TRUE
      ret
    }, logical(4))
    nTilesWithColour <- rowSums(coloursOnTile)
    nColoursOnTiles  <- colSums(coloursOnTile)
    commonColours <- which(nTilesWithColour == max(nTilesWithColour))
    
    valid.pairings <- switch(max(nTilesWithColour),
                             2:4, # 1: Each tile is exclusively a different colour
                             {# 2:
                               switch(length(commonColours), 2:4,
                                      which(coloursOnTile[commonColours[ifelse(coloursOnTile[commonColours[1], 1], 1, 2)], ])[-1],
                                      NULL, NULL)
                             }, { # 3: Three tiles share a colour
                               switch(sum(nTilesWithColour > 0), 2:4, return(NULL))
                             },
                             2:4) # 4: All tiles share a colour
    possibleLinks <- lapply(1:4, function (i) which(edges[i, ] + c(rep(2, 4), rep(0, 31)) > 0))
    linkOptions   <- vapply(possibleLinks, length, integer(1))
    linkConfigs   <- matrix(0, prod(linkOptions), ncol=4)
    linkConfigs[, 1] <- rep(possibleLinks[[1]], each=prod(linkOptions[2:4]))
    linkConfigs[, 2] <- rep(possibleLinks[[2]], each=prod(linkOptions[3:4]))
    linkConfigs[, 3] <- rep(possibleLinks[[3]], each=prod(linkOptions[4:4]))
    linkConfigs[, 4] <- possibleLinks[[4]]
    DoubleDummy <- function (edgesToBreak, sisterToOne) {
      # edgesToBreak is a list of the type of edge (column number in edgeNames) to be broken on each of tiles 1:4
      # The function will be sent the rows of linkConfigs as edgesToBreak
      # sisterToOne must be a valid sister tile;
      #   i.e. if 1&3 and 2&4 each share a distinct colour, 2 cannot be sisterToOne
      near.pair <- c(1, sisterToOne)
      far.pair  <- (1:4)[-near.pair]
      tilesInOrder <- c(near.pair, far.pair)
      tile.position <- match(1:4, tilesInOrder)
      sisterOf <- vapply(1:4, function (i) as.integer(tilesInOrder[switch(tile.position[i], 2, 1, 4, 3)]), integer(1))
      validConstruction <- vapply(commonColours, function (colour) {
        all(colourOn[colour, edgesToBreak[coloursOnTile[colour, ]]]) # is the colour 'colour' on the edges to break?
      }, logical(1))
      if (!all(validConstruction)) return (NULL)
      farNode   <- MeetEdges(edgesToBreak[ far.pair])
      nearNode  <- MeetEdges(edgesToBreak[near.pair])
      nearNearFarFar <- c(nearNode, nearNode, farNode, farNode)
      adjacentNode <- nearNearFarFar[tile.position]
      dummyEdge <- MeetEdges(c(farNode, nearNode))
      newEdges <- EditEdges(c(+1, 'DDmDDm', dummyEdge), edges)
      nearNode.easy <- nearNode <= 4
      outside.endOfMarginEdge <- switch(edge.types[dummyEdge], { # In(dummyEdge)
        rep(dummyEdge, 4) # To reconcile in ReconcileTile
      }, { # Btwn(farNode, nearNode)
        nearNearFarFar
      }, { # Tdm (easyNode, hardNode)
        if (nearNode.easy) {
          c(nearNode, nearNode, edgesToBreak[far.pair[1]], edgesToBreak[far.pair[2]])
        } else {
          c(edgesToBreak[near.pair[1]], edgesToBreak[near.pair[2]], farNode, farNode)
        }
      }, { # DDm (easyNode, one.hardTile)
        if (nearNode.easy) {
          if (edgesToBreak[far.pair[1]] <= 4) {
            c(nearNode, nearNode, edgesToBreak[far.pair[1]], 'DDmDDm')
          } else {
            c(nearNode, nearNode, 'DDmDDm', edgesToBreak[far.pair[2]])
          }
        } else {
          if (edgesToBreak[near.pair[1]] <= 4) {
            c(edgesToBreak[near.pair[1]], 'DDmDDm', farNode, farNode)
          } else {
            c('DDmDDm', edgesToBreak[near.pair[2]], farNode, farNode)
          }
        }
      }, { # DDmDDm
        c(edgesToBreak[near.pair[1]], edgesToBreak[near.pair[2]],
          edgesToBreak[ far.pair[1]], edgesToBreak[ far.pair[2]])
      })
      
      ReconcileTile <- function (i, tmpEdges) {
        # Replace the broken edge with two new edges, reclassifying TDms as DDms if necessary
        # tile:   the number of the tile we're working on at the moment
        tile <- tilesInOrder[i]
        # broken: the colouration of the edge upon that tile to be broken and replaced by new ones
        broken <- edgesToBreak[tile]
        # joinTo:  the colouration of the edge to be added, which will break the broken edge
        joinTo <- outside.endOfMarginEdge[i]
        # tmpEdges:  tile/colouration matrix.
        switch(edge.types[broken], { # 1: In()
          EditEdges(c(#-1, tile, edgeNames[broken],  ## TODO: This, and the below, is commented out in a crude attempt to handle Singletons.
            +1, DDm(tile, sisterOf[tile]), Btwn(c(broken, joinTo)),
            #+1, tile, edgeNames[broken]
            +1, tile, edgeNames[broken]), tmpEdges)
        }, { # 2: Btwn()
          brokenEnds <- endsOf[, broken]
          if (joinTo %in% brokenEnds) {
            EditEdges(c(-1, tile, edgeNames[broken],
                        +1, DDm(tile, sisterOf[tile]), joinTo,
                        +1, tile, Btwn(c(brokenEnds[1], joinTo)),
                        +1, tile, Btwn(c(brokenEnds[2], joinTo))), tmpEdges)
          } else { switch(EdgeTypes(1:5)[joinTo], { # In
            EditEdges(c(-1, tile, edgeNames[broken]
                        +1, DDm(tile, sisterOf[tile]), TDm(joinTo, brokenEnds),
                        +1, tile, TDm(brokenEnds[1], c(brokenEnds[2], joinTo)),
                        +1, tile, TDm(brokenEnds[2], c(brokenEnds[1], joinTo))), tmpEdges)
          }, { #Btwn
            stop("Program error: tried to join Btwn to Btwn")
          }, { # TDm
            EditEdges(c(-1, tile, edgeNames[broken]
                        +1, DDm(tile, sisterOf[tile]),  endsOf[1, joinTo],
                        +1, tile, Btwn(c(brokenEnds[1], endsOf[1, joinTo])),
                        +1, tile, Btwn(c(brokenEnds[2], endsOf[1, joinTo]))), tmpEdges)
          }, { # DDm
            stop("Program error: tried to join DDm to Btwn")
          }, { # DDmDDm
            # If we've added a double dummy, it must have intersected the only boundary on the tile!
            EditEdges(c(-1, tile, edgeNames[broken]
                        +1, DDm(tile, sisterOf[tile]), 'DDmDDm',
                        +1, tile, DDm(brokenEnds[1], brokenEnds[2]),
                        +1, tile, DDm(brokenEnds[2], brokenEnds[1])), tmpEdges)
          })
          }
        }, { # 3: TDm()
          tdCols <- tripleDummyEnds[, broken]
          if (any(which(colourOn[, joinTo]) %in% tdCols)) {
            colour.inCommon <- which(colourOn[, joinTo])[which(colourOn[, joinTo]) %in% tdCols]
            if (colour.inCommon == tdCols[1]) {
              EditEdges(c(+1, tile, colour.inCommon,
                          +1, DDm(tile, sisterOf[tile]), colour.inCommon
              ), tmpEdges)
            } else {
              EditEdges(c(-1, tile, TDm(tdCols[1], tdCols[1:3][-1]),
                          +2, tile, In(tdCols[1]), # Replacing the above, and the new inner edge
                          +1, DDm(tile, sisterOf[tile]), In(tdCols[1]),
                          -1, tile, TDm(tdCols[2], tdCols[1:3][-2]),
                          +1, tile, Btwn(tdCols[c(1, 2)]),
                          -1, tile, TDm(tdCols[3], tdCols[1:3][-3]),
                          +1, tile, Btwn(tdCols[c(1, 3)])), tmpEdges)
            }
          } else {
            EditEdges(c(-1, tile, TDm(tdCols[1], tdCols[1:3][-1]),
                        -1, tile, TDm(tdCols[2], tdCols[1:3][-2]),
                        -1, tile, TDm(tdCols[3], tdCols[1:3][-3]),
                        +1, tile, DDm(tdCols[1], joinTo), # Broken edge is sister to the Other (joined to) Colour
                        +1, DDm(tile, sisterOf[tile]), DDm(joinTo, tdCols[1]),
                        +1, tile, DDm(tdCols[2], tdCols[3])
                        +1, tile, DDm(tdCols[3], tdCols[2])
                        +1, tile, 'DDmDDm'), tmpEdges)
          }
        }, { # 4: DDm() # Note: comma added after In(joinTo) and removed on subsequent line
          # This may have introduced an error?
          EditEdges(c(+1, tile, In(joinTo),
                      +1, DDm(tile, sisterOf[tile]) ##### Edge to double dummy ... #####
          ), tmpEdges)
        }, { # 5: DDmDDm
          # Assumption: joinTo is a single colour, as all four colours are on this tile.
          # We therefore resolve the double dummy, and create a new triple dummy.
          newTripleColours <- (1:4)[-joinTo]
          sisterToJoinTo <- newTripleColours[which(tmpEdges[tile, ColNo(DDm(joinTo, newTripleColours[1])):ColNo(DDm(joinTo, newTripleColours[3]))] > 0)]
          these <- c(joinTo, sisterToJoinTo)
          others <- (1:4)[-these]
          EditEdges(c(-1, tile, edgeNames[broken],
                      -1, tile, DDm(these[1], these[2]),
                      -1, tile, DDm(these[2], these[1]),
                      -1, tile, DDm(others[1], others[2]),
                      -1, tile, DDm(others[2], others[1]),
                      +1, DDm(tile, sisterOf[tile]), joinTo,
                      +2, tile, joinTo, # one of which replaces the double dummy
                      +1, tile, Btwn(these),
                      +1, tile, TDm(joinTo, others),
                      +1, tile, TDm(others[1], c(joinTo, others[2])),
                      +1, tile, TDm(others[2], c(joinTo, others[1]))), tmpEdges)
        })
      }
      for (i in 1:4) newEdges <- ReconcileTile(i, newEdges)
      newEdges
    }
    unlist(apply(linkConfigs, 1, function(config) lapply(valid.pairings,
                                                         function (sisterToOne) DoubleDummy(config, sisterToOne))), recursive=FALSE)
  }), recursive=FALSE)
}


NumberOfTrees <- function (this.e) {
  inTileRoots <- btwnTileRoots <- tilesLinked.byCol <- subtreeEdges <- subtreeTips <- rootsTo.add <- this.e[1:4, 1:4]
  for (tile in 1:4) for (col in 1:4) {
    inTileRoots      [tile, col] <- sum(this.e[tile, edgeOf[, col]])
    tilesLinked.byCol[tile, col] <- sum(this.e[edgeOf[, tile], col])
    btwnTileRoots    [tile, col] <- sum(this.e[edgeOf[, tile], col],
                                        this.e[edgeOf[, tile], edgeOf[, col]])
  }
  subtreeRoots <- inTileRoots + btwnTileRoots
  subtreeRoots[subtreeEdges == -3] <- 0 # Resolve ambiguous edges such as tile=Btwn23, col=Btwn23; col 2 is either on tile 2 or 3, not both
  subtreeTips[] <- pmax(0, ((subtreeEdges - subtreeRoots) + 3) / 2)
  
  # Start assembling trees.  First build and label each individual subtree (i.e. tile-colour combination)
  edges.beforeRooting <- 2 * subtreeTips - 3
  waysToBuildUnrootedSubtrees <- DoubleFactorial(2 * subtreeTips - 5)
  
  # Now connect the subtrees on each individual tile
  rootsTo.add[] <- t(vapply(1:4, function (myTile) {
    tile <- this.e[myTile, -(1:4)]
    vapply(1:4, function (myCol)
      as.integer(sum(tile[edgeOf[, myCol]])),
      integer(1))
  }, integer(4)))
  internalEdgesAfterRootingTile <- edges.beforeRooting + rootsTo.add # Edges where both ends are the same colour
  externalEdgesAfterRootingTile <- rootsTo.add # Edges where only one end is the same colour
  additions <- rootsTo.add + tilesLinked.byCol
  all.waysToLinkSubtrees <- DoubleFactorial(edges.beforeRooting + (2 * (additions - 1))) / DoubleFactorial(edges.beforeRooting - 2)
  
  
  invalidBecause.ddm <- DoubleFactorial(edges.beforeRooting + pmax(0, (2 * (additions - 3)))) / DoubleFactorial(edges.beforeRooting - 2)
  invalidBecauseRoot.tdm <- DoubleFactorial(edges.beforeRooting + pmax(0, (2 * (additions - 2)))) / DoubleFactorial(edges.beforeRooting - 2) - ifelse(rootsTo.add == 3, invalidBecause.ddm, 0)
  invalidBecauseTileLink.tdm <- DoubleFactorial(edges.beforeRooting + pmax(0, (2 * (additions - 2)))) / DoubleFactorial(edges.beforeRooting - 2) - ifelse(tilesLinked.byCol == 3, invalidBecause.ddm, 0)
  invalidBecause.tdmAnd.tdm <- DoubleFactorial(edges.beforeRooting + pmax(0, (2 * (additions - 3)))) / DoubleFactorial(edges.beforeRooting - 2)
  invalidBecause.tdmAnd.ddm <- DoubleFactorial(edges.beforeRooting + pmax(0, (2 * (additions - 4)))) / DoubleFactorial(edges.beforeRooting - 2)
  invalidBecause.ddmAnd.tdm <- DoubleFactorial(edges.beforeRooting + pmax(0, (2 * (additions - 4)))) / DoubleFactorial(edges.beforeRooting - 2)
  invalidBecause.ddmAnd.ddm <- DoubleFactorial(edges.beforeRooting + pmax(0, (2 * (additions - 5)))) / DoubleFactorial(edges.beforeRooting - 2)
  
  
  possibleRoot.ddms <- ifelse(rootsTo.add >= 3, choose(rootsTo.add, 2), 0)
  possibleTileLink.ddms <- ifelse(tilesLinked.byCol >= 3, choose(tilesLinked.byCol, 2), 0)
  possibleRoot.tdms <- choose(rootsTo.add, 2)
  possibleTileLink.tdms <- choose(tilesLinked.byCol, 2)
  
  possibleRoot.tdmTileLink.tdm <- possibleRoot.tdms * possibleTileLink.tdms
  possibleRoot.tdmTileLink.ddm <- possibleRoot.tdms * possibleTileLink.ddms
  possibleRoot.ddmTileLink.tdm <- possibleRoot.ddms * possibleTileLink.tdms
  possibleRoot.ddmTileLink.ddm <- possibleRoot.ddms * possibleTileLink.ddms
  
  countedTwice <-
    possibleRoot.tdms * invalidBecauseRoot.tdm +
    possibleRoot.ddms * invalidBecause.ddm +
    possibleTileLink.tdms * invalidBecauseTileLink.tdm +
    possibleTileLink.ddms * invalidBecause.ddm -
    possibleRoot.tdmTileLink.tdm * invalidBecause.tdmAnd.tdm -
    possibleRoot.tdmTileLink.ddm * invalidBecause.tdmAnd.ddm -
    possibleRoot.ddmTileLink.tdm * invalidBecause.ddmAnd.tdm -
    possibleRoot.ddmTileLink.ddm * invalidBecause.ddmAnd.ddm
  
  waysToCompleteTiles <- all.waysToLinkSubtrees - countedTwice
  
  tileMarginEdges <- matrix(0, 4, 4)
  # Now connect tiles together
  waysToJoinTiles <- 1
  
  links <- which(t(this.e[-(1:4), ]) > 0)
  linkTiles  <- links %/% ncol(this.e) + 5 # tile
  linkColours <- links %% ncol(this.e)      # colour
  
  # Link tiles that share a colour first
  links <- c(links[linkColours <= 4], links[linkColours > 4])
  for (link in links) {
    tile   <- ((link - 1) %/% ncol(this.e)) + 5
    colour <- link %% ncol(this.e)
    tiles   <- which(memberOf[, tile])
    colours <- which(memberOf[, colour])
    
    for (myTile in tiles) {
      edgesOnTile <- internalEdgesAfterRootingTile[myTile, ] + externalEdgesAfterRootingTile[myTile, ]
      edgesTouchingTile <- edgesOnTile + tileMarginEdges[myTile, ]
      myColour <- colours[edgesOnTile[colours] > -3]
      tipsOfColour  <- subtreeTips[myTile, myColour]
      if (!length(tipsOfColour)) tipsOfColour <- 0
      
      if (colour <= 4) {
        waysToJoinTiles <- c(waysToJoinTiles, 1)        # we've already counted this, above
        internalEdgesAfterRootingTile[myTile, myColour] <- internalEdgesAfterRootingTile[myTile, myColour] + 1
        tileMarginEdges[myTile, myColour] <- tileMarginEdges[myTile, myColour] + 1
      } else if (FALSE) { #  Previous code
        # Single colour; easy
        if (tipsOfColour == 1) {
          waysToJoinThisTile <- switch (sum(this.e[myTile, c(Btwn(c(colour, 1)), Btwn(c(colour, 2)), Btwn(c(colour, 3)), Btwn(c(colour, 4)))[-colour]]) + 1,
                                        1,              # 0: perhaps a single tip in a triple dummy
                                        max(1, edgesOnTile[myColour]), # 1  ## EDITED! Could cause probs!
                                        2, # 2
                                        {
                                          message ('not sure about 3!')
                                          2 # 3
                                        })
          waysToJoinTiles <- c(waysToJoinTiles, waysToJoinThisTile)
          tileMarginEdges[myTile, myColour] <- tileMarginEdges[myTile, myColour] + 1
          internalEdgesAfterRootingTile   [myTile, myColour] <- internalEdgesAfterRootingTile   [myTile, myColour] + 1
        } else {
          # The colour IS already on the tile or 'colour' would be > 4
          #waysToJoinTiles <- c(waysToJoinTiles, edgesOnTile[myColour])
          #message('Check: Should this be internalEdgesAfterRootingTile instead of edgesOnTile?')
          waysToJoinTiles <- c(waysToJoinTiles, edgesOnTile[myColour])
          internalEdgesAfterRootingTile   [myTile, myColour] <- internalEdgesAfterRootingTile   [myTile, myColour] + 1
          tileMarginEdges[myTile, myColour] <- tileMarginEdges[myTile, myColour] + 1
        }
      } else if (colour <= 10) {
        waysToJoinTiles <- c(waysToJoinTiles, pmax(internalEdgesAfterRootingTile[myTile, myColour], 1))
        internalEdgesAfterRootingTile   [myTile, myColour] <- internalEdgesAfterRootingTile   [myTile, myColour] + 1
        tileMarginEdges[myTile, myColour] <- tileMarginEdges[myTile, myColour] + 1
      } else if (colour <= 22) {
        # Find out if colour is on this tile;
        if (tipsOfColour) {
          # If colour is on this tile: easy
          waysToJoinTiles <- c(waysToJoinTiles, pmax(internalEdgesAfterRootingTile[myTile, myColour], 1))
          internalEdgesAfterRootingTile   [myTile, myColour] <- internalEdgesAfterRootingTile   [myTile, myColour] + 1
          tileMarginEdges[myTile, myColour] <- tileMarginEdges[myTile, myColour] + 1
        } else {
          # There are already two edges leading to the dummy node that is the only place for this to go
          ###waysToJoinTiles <- c(waysToJoinTiles, 1) # TODO when debugging complete: DELETE
          tileMarginEdges[myTile, myColour] <- tileMarginEdges[myTile, myColour] + 1
        }
      } else if (colour <= 34) {
        # Find out if colour is on this tile;
        if (tipsOfColour) {
          # If colour is on this tile: easy
          waysToJoinTiles <- c(waysToJoinTiles, pmax(internalEdgesAfterRootingTile[myTile, myColour], 1))
          internalEdgesAfterRootingTile   [myTile, myColour] <- internalEdgesAfterRootingTile   [myTile, myColour] + 1
          tileMarginEdges[myTile, myColour] <- tileMarginEdges[myTile, myColour] + 1
        } else {
          # If not: there's already an edge leading to the dummy node
          ###waysToJoinTiles <- c(waysToJoinTiles, 1) # TODO when debugging complete: DELETE
          tileMarginEdges[myTile, myColour] <- tileMarginEdges[myTile, myColour] + 1
        }
      } else {
        # All necessary nodes already exist.
        ###waysToJoinTiles <- c(waysToJoinTiles, 1) # TODO when debugging complete: DELETE
        tileMarginEdges[myTile, myColour] <- tileMarginEdges[myTile, myColour] + 1
      }
      #message ('Tile ', myTile, ' (one of ', length(tiles), '): ', paste0(waysToJoinTiles, collapse = ' '))
      waysToJoinTiles
    }
  }
  prod(waysToBuildUnrootedSubtrees, waysToCompleteTiles, waysToJoinTiles)
}

# Housekeeping functions
RemoveNullsFromList <- function (L) L[!vapply(L, is.null, logical(1))]

#' Number of trees consistent with two splits
#'
#' @inheritParams MutualInformation
#' @author Martin R. Smith
#' @concepts Split information
#' @export
TreesConsistentWithTwoSplits <- function (n, A1, A2=A1) {
  smallSplit <- min(A1, A2)
  bigSplit <- max(A1, A2)
  
  zones <- matrix(c(smallSplit, bigSplit - smallSplit, 0L, n - bigSplit), 2L, 2L)
  
  if (all(zones %in% c(0L, n))) { # Splits are uninformative
    return(NUnrooted(n))
  }
  if (smallSplit == 0 || bigSplit == n) {
    return(prod(vapply(zones[zones > 0], NRooted, double(1))))
  }
  
  zoneEdges <- 2L * zones - 3L # Beware the -3!
  tipsOnTile <- rowSums(zones)
  edgesOnTile <- 2L * tipsOnTile - 3L
  coloursOnTile <- zones > 0
  nColOnTile <- rowSums(coloursOnTile)
  nTiles <- nrow(zones)
  twoSplitEdgeNames <- c("In 1", "In 2", "Btwn 12", "DDm 1.2")
  startEdges <- matrix(0, length(twoSplitEdgeNames), length(twoSplitEdgeNames),
                       dimnames=list(twoSplitEdgeNames, twoSplitEdgeNames))
  
  startEdges[1:2, 1:2] <- zoneEdges
  
  # Local funcs
  ColourInCommon  <- function(tile1, tile2) zones[tile1, ] > 0 & zones[tile2, ] > 0
  ColourOnTile    <- function(colour, tile) zones[tile, colour] > 0
  
  # Count assemblies from within tiles
  edgesAfterTileInternalsAssembled <- lapply(seq_len(2L), function (i) {
    # edgeOptions is a list of all possible edge configurations.
    # LEMMA: We can work out the number of consistent trees from edgeOptions
    if (nColOnTile[i] == 1) {
      list(startEdges[i, ]) # One colour
    } else {
      colsOnThisTile <- which(coloursOnTile[i, ])
      list(AddBoundaryToTile(colsOnThisTile[1:2], startEdges[i, ])) # Two colours
    }
  })
  validAssemblies <- lapply(edgesAfterTileInternalsAssembled,
                            function (x) which(!vapply(x, is.null, logical(1))))
  nAssemblies <- vapply(validAssemblies, length, integer(1))
  nCombinations <- prod(nAssemblies[nAssemblies > 0])
  combinations <- t(matrix(vapply(validAssemblies[nAssemblies > 0],
                                  rep_len, integer(nCombinations),
                                  nCombinations),
                           nrow=nCombinations))
  assemblies <- unique(lapply(as.list(data.frame(combinations)), function (x) {
    # assemblies is a list containing an item for every permutation of tiles' internal assemblies
    edges <- startEdges
    for (tile in seq_len(nTiles)) {
      edges[tile, ] <- edgesAfterTileInternalsAssembled[[tile]][[x[tile]]]
    }
    edges
  }))
  
  # Now for each possible assembly within each tile, connect the tiles together:
  edgesAfterTilesConnected <- unique(RemoveNullsFromList(
    if (nTiles == 1) {
      assemblies
    } else {
      unlist(lapply(assemblies, function (ass) LinkTwoTiles(ass)),
             recursive = FALSE)
    }
  )) # unique is neccessary; see dat.27 for example case
  
  validTrees <- vapply(edgesAfterTilesConnected, NumberOfTrees, double(1))
  
  # Return:
  round(sum(validTrees))
}

#' Number of unrooted trees consistent with splits
#' @template splitsParam
#' 
#' @references 
#' \insertRef{Carter1990}{TreeSearch}, Theorem 2.
#'
UnrootedTreesMatchingSplit <- function (splits) {
  splits <- splits[splits > 0L]
  totalTips <- sum(splits)
  round(prod(DoubleFactorial(2L * totalTips - 5L),
             DoubleFactorial(2L * splits - 3L)) /
          DoubleFactorial(2L * (totalTips - length(splits)) - 1L))
}

#' Double Factorial
#' 
#' @param n Vector of integers.
#' 
#' @return Calculates the double factorial, n x (n - 2) x (n - 4) x (n - 6) x ...
#' 
#' @examples {
#' DoubleFactorial (0) # Return 1 if n < 2
#' DoubleFactorial (2) # 2
#' DoubleFactorial (5) # 1 x 3 x 5
#' DoubleFactorial (8) # 2 x 4 x 6 x 8
#' }
#' 
#' @author Martin R. Smith
#' @export
DoubleFactorial <- function (x) exp(LogDoubleFactorial(x))

# Memoize this function at your peril...
#' @importFrom phangorn ldfactorial
#' @describeIn DoubleFactorial Returns the logarithm of the double factorial.
LogDoubleFactorial <- (function (x) {
  ifelse(x < 2, 0,
         ifelse(x %% 2,
                ldfactorial(x),
                lfactorial(x) - ldfactorial(x - 1L)
         )
  )
})
