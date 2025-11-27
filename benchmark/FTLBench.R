library("TreeTools")
library("TreeSearch")

tree <- ape::read.tree(
  text = ("((a1, a2), (((b1, b2), (c, d)), ((e1, e2), (f, g))));"))
char <- "0102220333"
charDat <- StringToPhyDat(char, TipLabels(tree))
if (interactive()) {
  PlotCharacter(tree, charDat)
}
nullLengths <- c(mean = 5.8207, median = 6)
expect_equal(
  Consistency(StringToPhyDat(char, TipLabels(tree)), tree, nRelabel = 100),
  c(ci = m / s, ri = r, rc = r * m / s,
    rhi = h / (nullLengths[["median"]] - m),
    rhiBar = h / (nullLengths[["mean"]] - m)),
  tolerance = 0.01
)

expect_equal(
  Consistency(StringToPhyDat(char, TipLabels(tree)), tree, nRelabel = NULL),
  c(ci = m / s, ri = r, rc = r * m / s,
    rhi = h / (nullLengths[["median"]] - m),
    rhiBar = h / (nullLengths[["mean"]] - m)),
  tolerance = 0.01
)
pd <- StringToPhyDat(char, TipLabels(tree))

bench::mark(Consistency(pd, tree, nRelabel = 1000),
            Consistency(pd, tree, nRelabel = NULL), check = FALSE)
