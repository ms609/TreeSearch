context('Tree differences')

test_that("Split combatibility is correctly established", {
  expect_true(SplitsCompatible(as.logical (c(0,0,1,1,0)), as.logical(c(0,0,1,1,0))))
  expect_true(SplitsCompatible(as.logical (c(0,0,1,1,0)), !as.logical(c(0,0,1,1,0))))
  expect_true(SplitsCompatible(as.logical (c(0,0,1,1,0)), as.logical(c(1,0,1,1,0))))
  expect_true(SplitsCompatible(!as.logical(c(0,0,1,1,0)), as.logical(c(1,0,1,1,0))))
  expect_false(SplitsCompatible(as.logical(c(0,0,1,1,0)), as.logical(c(1,1,0,1,0))))
})

test_that('Tree differences are correctly calculated', {
  # Labels in different order to confound Tree2Splits
  treeSym8 <- ape::read.tree(text='((e, (f, (g, h))), (((a, b), c), d));')
  treeBal8 <- ape::read.tree(text='(((e, f), (g, h)), ((a, b), (c, d)));')
  treeAb.Cdefgh <- ape::read.tree(text='((a, b), (c, d, e, f, g, h));')
  treeAbc.Defgh <- ape::read.tree(text='((a, b, c), (d, e, f, g, h));')
  treeAbcd.Efgh <- ape::read.tree(text='((a, b, c, d), (e, f, g, h));')
  treeTwoSplits <- ape::read.tree(text="(((a, b), c, d), (e, f, g, h));")

  # Labels differ
  expect_error(InfoTreeDist(treeSym8, ape::read.tree(text='((a, b, c, D), (e, f, g, h));')))
  expect_equal(22.53747, round(InfoTreeDist(treeSym8, treeSym8), 5))
  expect_equal(13.75284, round(InfoTreeDist(treeSym8, treeBal8), 5))
  expect_equal(-log2(945/10395), InfoTreeDist(treeSym8, treeAb.Cdefgh))
  expect_equal(-log2(315/10395), InfoTreeDist(treeSym8, treeAbc.Defgh))
  # Test symmetry of small vs large splits
  expect_equal(InfoTreeDist(treeSym8, treeAbc.Defgh), InfoTreeDist(treeAbc.Defgh, treeSym8))
  expect_equal(-log2(225/10395), InfoTreeDist(treeSym8, treeAbcd.Efgh))
  expect_equal(-log2(225/10395) - log2(945/10395),
               InfoTreeDist(treeSym8, treeTwoSplits))
  expect_equal(MutualInformation(8, 4, 3),
               InfoTreeDist(treeTwoSplits, treeAbc.Defgh))
  
  expect_equal(InfoTreeDist(treeSym8, list(treeSym8, treeBal8)), InfoTreeDist(list(treeSym8, treeBal8), treeSym8))
  expect_equal(matrix(c(InfoTreeDist(treeSym8, treeSym8), InfoTreeDist(treeBal8, treeSym8),
                      InfoTreeDist(treeSym8, treeAbc.Defgh), InfoTreeDist(treeBal8, treeAbc.Defgh),
                      InfoTreeDist(treeSym8, treeAbcd.Efgh), InfoTreeDist(treeBal8, treeAbcd.Efgh)),
                      3L, 2L, byrow=TRUE,
                      dimnames=list(c('sym', 'abc', 'abcd'), c('sym', 'bal'))), 
    InfoTreeDist(list(sym=treeSym8, bal=treeBal8), 
               list(sym=treeSym8, abc=treeAbc.Defgh, abcd=treeAbcd.Efgh)))
  
})
