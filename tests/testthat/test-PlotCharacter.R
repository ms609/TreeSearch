test_that("PlotCharacter.phylo()", {
  
  dataset <- TreeTools::StringToPhyDat("1111 1111 0000", tips = 12)
  expect_error(PlotCharacter(TreeTools::BalancedTree(14), dataset),
               "Taxa in tree missing from dataset:\\s*t13, t14$")
  expect_error(PlotCharacter(TreeTools::StarTree(12), dataset),
               "bifurcating")
  
  Character <- function (str, plot = FALSE, edges = FALSE, ...) {
    tree <- ape::read.tree(text = 
     "((((((a, b), c), d), e), f), (g, (h, (i, (j, (k, l))))));")
    if (edges) {
      tree$edge.length <- rep(2, 22)
    }
    dataset <- TreeTools::StringToPhyDat(str, tips = tree)
    PlotCharacter(tree, dataset,
                  edge.width = 3, plot = plot, ...)
  }
  
  expect_equal(
    Character("23--1??--032", updateTips = TRUE),
    structure(
      c(FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE,
        TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE,
        TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
        FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE,
        FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE,
        FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE,
        FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
        FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE,
        FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE,
        FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE,
        TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
        FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
        TRUE, FALSE, FALSE, FALSE, TRUE, TRUE), .Dim = c(23L, 5L),
      .Dimnames = list(NULL, c("-", "0", "1", "2", "3"))))
  
  skip_if_not_installed("vdiffr")

  Test <- if (interactive()) {
    function (str, edges = FALSE, ...) {
      invisible(Character(str, plot = TRUE, edges = edges, ...))
    }
  } else {
    function (str, edges = FALSE, ...) {
      vdiffr::expect_doppelganger(
        paste0("PlotChar_",
               gsub("?", "Q",
                    gsub("(", "d",
                         gsub(")", "b",
                              gsub("-", "I", str,
                                   fixed = TRUE), fixed = TRUE), fixed = TRUE), fixed = TRUE)),
        function() Character(str, plot = TRUE, edges = edges, ...))
    }
  }
  
  Test("23--1??--032", edges = TRUE)
  Test("23--1??(-0)-(01)32")
  Test("23??1????032")
  Test("11--????--11", unitEdge = TRUE)
  Test("000011????00")
  Test("????????????")
  Test("-------?????")
  Test("------------")
  Test("1234(45)AACGTTT")
  
  # From TGuillerme testing suite:
  Test("11-------100")
  Test("1100----1100")
  Test("000011110000")
  Test("1---1111---1")
  Test("----1111---1")
  Test("01----010101")
  Test("01---1010101")
  Test("1??--??--100")
  Test("21--3??--032")
  Test("11--1??--111")
  Test("11--1000001-")
  Test("01------0101")
  Test("110--?---100")
  Test("210--100--21")
  Test("????----1???")
  Test("23--1----032")
  Test("1----1----1-")
  Test("-1-1-1--1-1-")
  
  Test("--------0101")
  Test("10101-----01")
  Test("011--?--0011")
  Test("110--??--100")
  Test("21--1----012")
  Test("11----111111")
  Test("210210------")
  Test("----1111----")
  Test("230--??1--32")
  Test("023--??1--32")
  Test("023-???1--32")
  Test("23--1?1--023")
  Test("----1010----")
  Test("------11---1")
  Test("10----11---1")
  Test("320--??3--21")
})

test_that("Edge cases work", {
  tree <- ape::read.tree(text = "(a, (b, ((c, d), (e, f))));")
  dataset <- TreeTools::StringToPhyDat("-01100", tips = tree)
  if (interactive()) {
    PlotCharacter(tree, dataset)
  } else {
    expect_equal(c("-" = FALSE, "0" = TRUE, "1" = FALSE),
                 PlotCharacter(tree, dataset, plot = FALSE)[9, ])
  }
  
  tree <- ape::read.tree(text = "(a, (b, (c, (d, (e, f)))));")
  dataset <- TreeTools::StringToPhyDat("--0101", tips = tree)
  if (interactive()) {
    PlotCharacter(tree, dataset)
  } else {
    expect_equal(cbind("-" = c(1, 1, 0, 0, 0),
                       "0" = c(0, 0, 1, 1, 1),
                       "1" = c(0, 0, 1, 1, 1)),
                 1 * PlotCharacter(tree, dataset, plot = FALSE)[7:11, ])
  }
})

test_that("Out-of-sequence works", {
  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger("PlotChar_out-of-sequence", function () {
    PlotCharacter(ape::read.tree(text = "(a, (b, (c, d)));"),
                  TreeTools::StringToPhyDat("1342",
                                            tips = c("a", "c", "d", "b"))
                  )}
  )
})

test_that("PlotCharacter.multi()", {
  Bal <- TreeTools::BalancedTree
  a..h <- letters[1:8]
  expect_error(PlotCharacter(list(Bal(8), 9), "dataset"), "class `phylo`")
  expect_error(PlotCharacter(list(Bal(8), Bal(9)), "dataset"),
               "same tip labels")
  expect_error(PlotCharacter(list(Bal(8), Bal(a..h)), "dataset"),
               "same tip labels")
  
  trees <- ape::read.tree(text = c("(a, (b, (c, (d, ((g, h), (e, f))))));",
                                   "(a, (b, (c, ((d, e), (f, (g, h))))));"))
  
  dat <-  TreeTools::StringToPhyDat("00011011", tips = a..h)
  expect_equal(PlotCharacter(trees[1], dat, plot = FALSE),
               PlotCharacter(trees[[1]], dat, plot = FALSE))
                             
                             
  state1 <- PlotCharacter(trees[[1]], dat, plot = FALSE)
  state2 <- PlotCharacter(trees[[2]], dat, plot = FALSE)
  stateCons <- PlotCharacter(trees, dat, plot = FALSE)
  expect_equal(stateCons, state1[-c(13, 15), ] | 
                 state2[c(match(TipLabels(trees[[1]]), TipLabels(trees[[2]])),
                          9:12, 15), ])
  
  
  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger("PlotChar_consensus", function() {
    PlotCharacter(trees, dat)
  })
  vdiffr::expect_doppelganger("PlotChar_invariant", function() {
    inv <- TreeTools::StringToPhyDat("00000000", tips = a..h)
    PlotCharacter(trees, inv)
  })
  vdiffr::expect_doppelganger("PlotChar_invar_ambig", function() {
    invq <- TreeTools::StringToPhyDat("000?00?{01}", tips = a..h)
    PlotCharacter(trees, invq)
  })
})
