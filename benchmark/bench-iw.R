library("TreeTools")

tmp_lib <- tempfile(pattern = "lib")
dir.create(tmp_lib)
devtools::install(args = paste0("--library=", tmp_lib))
library("TreeSearch", lib.loc = tmp_lib)

data("congreveLamsdellMatrices", package = "TreeSearch")
dataset <- congreveLamsdellMatrices[[42]]
tree <- AdditionTree(dataset) |> RenumberTips(names(dataset)) |> RootTree(1)

.IWScore <- function (edge, morphyObjs, weight, charSeq, concavity, 
                      minLength, target = Inf) {
  morphy_iw(edge, morphyObjs, weight, minLength, charSeq,
            concavity, target + epsilon)
}
nTip <- NTip(tree)
edge <- tree[["edge"]]

morphyObj <- PhyDat2Morphy(dataset)
on.exit(morphyObj <- UnloadMorphy(morphyObj), add = TRUE)

at <- attributes(dataset)
characters <- PhyToString(dataset, ps = "", useIndex = FALSE,
                          byTaxon = FALSE, concatenate = FALSE)
startWeights <- at[["weight"]]
minLength <- MinimumLength(dataset, compress = TRUE)
morphyObjects <- lapply(characters, SingleCharMorphy)
on.exit(morphyObjects <- vapply(morphyObjects, UnloadMorphy, integer(1)),
        add = TRUE)
  
nLevel <- length(at[["level"]])
nChar <- at[["nr"]]
nTip <- length(dataset)
cont <- at[["contrast"]]
colnames(cont) <- as.character(at[["levels"]])
simpleCont <- ifelse(rowSums(cont) == 1,
                     apply(cont != 0, 1, function (x) colnames(cont)[x][1]),
                     "?")

unlisted <- unlist(dataset, use.names = FALSE)
tokenMatrix <- matrix(simpleCont[unlisted], nChar, nTip)
charInfo <- apply(tokenMatrix, 1, CharacterInformation)

# Crude estimate of score added per unit processing time
rawPriority <- charInfo
priority <- startWeights * rawPriority
informative <- charInfo > 0
# Will work from end of sequence to start.
charSeq <- seq_along(charInfo)[informative][order(priority[informative])] - 1L
concavity <- 10

identical(
  morphy_iw(edge, morphyObjects, startWeights, minLength, charSeq, concavity,
            target = Inf),
  morphy_iw_charwise(edge, morphyObj, startWeights, minLength, charSeq, concavity,
                     target = Inf)
)

ub(
  morphy_iw(edge, morphyObjects, startWeights, minLength, charSeq, concavity,
            target = Inf),
  morphy_iw_charwise(edge, morphyObj, startWeights, minLength, charSeq,
                     concavity, target = Inf)
)
