# Extracted from test-ts-iw.R:119

# prequel ----------------------------------------------------------------------
skip_on_cran()
ts_iw <- function(tree, ds, min_steps, k) {
  TreeSearch:::ts_fitch_score(tree$edge, ds$contrast, ds$tip_data, ds$weight,
                              ds$levels, min_steps = min_steps, concavity = k)
}
result_phylo <- function(result, ref_tree) {
  structure(
    list(edge = result$edge, tip.label = ref_tree$tip.label,
         Nnode = ref_tree$Nnode),
    class = "phylo"
  )
}
iw_ref <- list(
  Vinther2008 = list(
    ew_pect = 140,
    pect = c(`3` = 16.0214285714, `10` = 6.5712620713, `100` = 0.7738979605),
    ew_rand = 206,
    rand = c(`3` = 22.6902597403, `10` = 10.5984709735, `100` = 1.3952233575)
  ),
  Agnarsson2004 = list(
    ew_pect = 1124,
    pect = c(`3` = 105.6437092319, `10` = 53.9811403206, `100` = 7.8903937162),
    ew_rand = 2025,
    rand = c(`3` = 140.6803269787, `10` = 84.6971261078, `100` = 15.5259369359)
  ),
  Wills2012 = list(
    ew_pect = 501,
    pect = c(`3` = 40.4743589744, `10` = 21.5671299686, `100` = 3.3844536694),
    ew_rand = 752,
    rand = c(`3` = 51.2034153662, `10` = 30.4867399246, `100` = 5.5136222612)
  ),
  Aria2015 = list(
    ew_pect = 184,
    pect = c(`3` = 18.8750000000, `10` = 8.7590840532, `100` = 1.1607895426),
    ew_rand = 311,
    rand = c(`3` = 28.6465091926, `10` = 15.4228339170, `100` = 2.3271752178)
  ),
  Zhu2013 = list(
    ew_pect = 2274,
    pect = c(`3` = 165.9651353803, `10` = 100.6971122772, `100` = 18.1011399968),
    ew_rand = 2219,
    rand = c(`3` = 166.2086376593, `10` = 100.1915980265, `100` = 17.6983351390)
  ),
  Loconte1991 = list(
    ew_pect = 1099,
    pect = c(`3` = 67.3159008309, `10` = 42.6637628190, `100` = 8.3203496488),
    ew_rand = 1121,
    rand = c(`3` = 67.2907477021, `10` = 42.7966001749, `100` = 8.4684795670)
  )
)
steps_ref <- list(
  # Updated 2026-03-19 after T-097 NA ambiguity fix (ts_na_char_steps)
  Vinther2008 = as.integer(c(0, 2, 1, 2, 1, 1, 1, 2, 1, 2, 3, 2, 3, 2, 2,
                  4, 4, 3, 3, 5, 2, 2, 2, 0, 3, 3, 3, 5, 3, 2, 2, 4, 2,
                  4, 3, 2, 2, 4, 3, 1, 0, 3, 0, 6, 2, 2, 2, 4, 4, 2)),
  Aria2015 = as.integer(c(2, 7, 2, 2, 9, 2, 3, 3, 6, 2, 4, 3, 2, 5, 2, 2,
               3, 2, 1, 3, 4, 5, 6, 4, 2, 3, 17, 8, 5, 2, 1, 2, 2, 2, 3,
               2, 6, 2, 4, 3, 2, 3, 5, 2, 1, 5, 5, 8, 3, 2))
)

# test -------------------------------------------------------------------------
skip_on_cran()
data("inapplicable.phyData", package = "TreeSearch")
for (nm in names(iw_ref)) {
    dataset <- inapplicable.phyData[[nm]]
    ds <- make_ts_data(dataset)
    minSteps <- MinimumLength(dataset, compress = TRUE)

    set.seed(5729)
    tree <- TreeTools::Preorder(TreeTools::RandomTree(dataset, root = TRUE))

    for (k_str in c("3", "10", "100")) {
      k <- as.numeric(k_str)
      score <- ts_iw(tree, ds, minSteps, k)
      expect_equal(score, iw_ref[[nm]]$rand[[k_str]], tolerance = 1e-8,
                   label = paste(nm, "rand k =", k))
    }
  }
