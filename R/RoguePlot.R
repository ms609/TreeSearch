#' Another approach: Calculate similarity of neighbouring clades
#' https://www.biorxiv.org/content/10.1101/2021.06.24.449756v1

#' Visualise position of rogue taxa
#' 
#' @param trees List or `multiPhylo` object containing phylogenetic trees
#' of class `phylo` to be summarized.
#' @param tip Numeric or character identifying rogue leaf, in format accepted
#' by `drop.tip()`.
#' @param p A numeric value between 0.5 and 1 giving the proportion for a clade
#' to be represented in the consensus tree (see `consensus()`).
#' @param \dots Additional parameters to `plot.phylo()`.
#' @param plot Logical specifying whether to plot the tree.
#' @param Palette Function that takes a parameter `n` and generates a colour 
#' palette with `n` entries.
#' @param thinEdge,fatEdge Numeric specifying width to plot edges if the rogue
#' tip does not / does attach to them.
#' @return `RoguePlot()` returns a vector of integers specifying the number of
#' trees in `trees` in which the rogue leaf is attached to each edge in turn
#' of the consensus tree.
#' @references 10.1371/journal.pone.0212942
#' @examples
#' trees <- list(read.tree(text = '(a1, (b1, (c1, (rogue, (d1, e1)))));'),
#'               read.tree(text = '(a1, (b1, (c1, (rogue, (d1, e1)))));'),
#'               read.tree(text = '(a1, (b1, (c1, (rogue, (d1, e1)))));'),
#'               read.tree(text = '(a1, (b1, (c1, (d1, (rogue, e1)))));'))
#' RoguePlot(trees, 'rogue')
#' @template MRS
#' @importFrom ape consensus drop.tip
#' @importFrom TreeTools AddTip as.Splits
#' @importFrom grDevices colorRampPalette
#' @export
RoguePlot <- function (trees, tip, p = 1, plot = TRUE,
                       Palette = colorRampPalette(c('black', 'red'), space = 'Lab'),
                       thinEdge = par('lwd'), fatEdge = thinEdge + 1L,
                       ...) {
  noRogue <- trees
  splits <- as.Splits(trees)
  noRogue[] <- lapply(trees, drop.tip, tip)
  cons <- consensus(noRogue, p = p)
  additionPoints <- lapply(cons$edge[, 2],
                           function (where) AddTip(cons, where, tip))
  nOnEdge <- vapply(
    as.Splits(additionPoints, splits),
    function (sp) sum(vapply(splits, function (st) all(sp %in% st), logical(1))),
    integer(1))
  
  if (plot) {
    plot(cons, edge.col = Palette(sum(nOnEdge) + 1L)[nOnEdge + 1L],
         edge.width = ifelse(nOnEdge > 0, fatEdge, thinEdge), ...)
  }
  
  # Return:
  nOnEdge
}