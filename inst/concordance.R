devtools::load_all()
use("ape")
use("TreeTools")

labels <- letters[1:5]
tr <- PectinateTree(labels)
plot(tr)
nodelabels(adj = c(4, 0.5))
 char <- MatrixToPhyDat(cbind(
  c(a = 0, b = 0, c = 0, d = 1, e = 1), # 100% concordant with both
  c(0, 0, 1, 0, 1), # discordant
  c(0, 1, 0, 0, 1)  
  ))
QuartetConcordance(tr, char, weight = FALSE)

Conc <- lapply(labels, function(lab) {
  quartet <- DropTip(tr, lab)
  for (i in 1:3) {
    message(setdiff(labels, lab), " char ", i, ": ", QuartetConcordance(quartet, char[, i]))
  }
})
  
