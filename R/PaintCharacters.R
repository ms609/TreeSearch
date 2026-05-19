#' Colour characters by tree concordance
#'
#' `PaintCharacters()` assigns a colour to each character in `dataset` by
#' computing a perceptually weighted mean of the colours assigned to tree edges
#' by [TreeTools::PaintTree()], using the mutual information between each
#' character and each edge as the weight.
#'
#' For each character, the colour is the weighted mean (in CIELAB space, which
#' is perceptually uniform) of the colours of all tree edges that the character
#' concordantly supports.  The weight for each edge is the product of its
#' normalized mutual information (concordance quality) and its relative
#' information amount; discordant edges (quality \eqn{\le 0}) are excluded.
#' Characters with no concordant signal on the tree are coloured grey
#' (`"#888888"`).
#'
#' The intended use is to colour a character-space MDS plot so that characters
#' supporting the same clade share a similar hue.  If the returned colours look
#' desaturated ("murky"), try raising `threshold` to exclude low-information
#' edges, or inspect `ConcordanceTable()` directly to understand the
#' character–edge signal.
#'
#' @param dataset A `phyDat` object containing morphological character data,
#'   whose `names` match the tip labels of `tree`.
#' @param tree A `phylo` object whose tip labels match `names(dataset)`.
#' @param threshold Numeric scalar; edges whose information value (the
#'   `"hBest"` × `"n"` product from [ClusteringConcordance()]) is below this
#'   threshold are excluded from the weighted average regardless of their
#'   concordance.  Default `0` retains all concordant edges.  Raising the
#'   threshold suppresses low-information edges that would otherwise dilute the
#'   colour signal.
#' @param palette Palette specification passed to [TreeTools::PaintTree()].
#'   Either a character string (`"default"`, `"protanopia"`, `"tritanopia"`)
#'   or a function `function(h, s)` mapping hue (0–360°) and saturation (0–1)
#'   to hex colours.
#'
#' @return A character vector of hex colour strings, one entry per character in
#'   `dataset`, named by character index.  Grey (`"#888888"`) indicates
#'   characters with no concordant signal on the tree.
#'
#' @examples
#' data("congreveLamsdellMatrices", package = "TreeSearch")
#' dataset <- congreveLamsdellMatrices[[1]][, 1:12]
#' tree <- referenceTree
#' library("TreeTools", quietly = TRUE)
#'
#' cols <- PaintCharacters(dataset, tree)
#' conc <- ConcordanceTable(tree, dataset)
#' # Plot the tree alongside to interpret the colours:
#' paint <- PaintTree(tree)
#' plot(tree, edge.color = paint$edgeCol, edge.width = 2)
#'
#' @seealso [TreeTools::PaintTree()], [ConcordanceTable()]
#' @family split support functions
#' @importFrom grDevices col2rgb convertColor rgb
#' @importFrom TreeTools PaintTree
#' @export
PaintCharacters <- function(dataset, tree, threshold = 0,
                             palette = "default") {
  paint <- TreeTools::PaintTree(tree, palette)
  cc    <- ClusteringConcordance(tree, dataset, return = "all")

  # Replicate the ConcordanceTable extraction (without triggering its plot).
  # matrix() guards against dimension collapse when nChar == 1 or nEdge == 1.
  nEdge   <- dim(cc)[[2L]]
  nChar   <- dim(cc)[[3L]]
  info    <- matrix(cc["hBest", , ] * cc["n", , ], nEdge, nChar,
                    dimnames = dimnames(cc)[2:3])
  relInfo <- info / max(info, na.rm = TRUE)
  relInfo[is.na(relInfo)] <- 0
  quality <- matrix(cc["normalized", , ], nEdge, nChar)
  relInfo[is.na(quality)] <- 0
  quality[is.na(quality)] <- 0

  # Align PaintTree edge colours to ClusteringConcordance edge order.
  # Row names are child node IDs (non-trivial splits only).
  ctNodes  <- as.integer(rownames(info))
  edgeIdx  <- match(ctNodes, tree[["edge"]][, 2L])
  edgeCols <- paint$edgeCol[edgeIdx]

  # Convert edge colours to CIELAB (perceptually uniform; a*/b* are Cartesian
  # so weighted averages avoid the circular-mean issue of hue).
  labMat <- matrix(
    convertColor(t(col2rgb(edgeCols)) / 255, from = "sRGB", to = "Lab"),
    ncol = 3L
  )  # nEdges × 3

  # Weight matrix: concordant edges only, scaled by relative information.
  wMat <- pmax(quality, 0) * relInfo   # nEdges × nChars
  wMat[info < threshold] <- 0

  wSum     <- colSums(wMat)            # nChars
  noInfo   <- wSum == 0
  wSumSafe <- ifelse(noInfo, 1, wSum)

  # Weighted Lab mean: t(3×nEdges %*% nEdges×nChars) → nChars×3, then / wSum.
  labAvg <- t(t(labMat) %*% wMat) / wSumSafe

  # Convert back to sRGB; clamp out-of-gamut values; encode as hex.
  # matrix() guards against convertColor() dropping to a vector for 1 row.
  rgbAvg   <- matrix(
    pmax(0, pmin(1, convertColor(labAvg, from = "Lab", to = "sRGB"))),
    ncol = 3L
  )
  charCols <- rgb(rgbAvg[, 1L], rgbAvg[, 2L], rgbAvg[, 3L])
  charCols[noInfo] <- "#888888"
  names(charCols) <- colnames(info)

  # Return:
  charCols
}
