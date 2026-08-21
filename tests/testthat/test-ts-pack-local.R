# Tests for TS_PACK_LOCAL per-block-local-alphabet packing (default-ON).
# Packing is score-EXACT per tree (a bijective per-block state relabel plus the
# "?"->char-alphabet substitution) but changes plane indexing, so:
#   (1) the scorer must return identical lengths ON vs OFF, and
#   (2) reconstruction output must map planes back to GLOBAL state labels
#       (CharBlock::plane_state[]), else a packed block prints the plane index
#       instead of the state it represents.
# The scorer is index-agnostic, so only (2) exercises plane_state[] -- hence the
# hand-built non-contiguous-alphabet fixture below.
skip_on_cran()

# --- (1) score-exactness: TreeLength identical ON vs OFF on real NA data ------

test_that("TS_PACK_LOCAL is score-exact per tree (NA dataset)", {
  skip_if_not_installed("TreeTools")
  data(inapplicable.phyData)
  phy <- inapplicable.phyData[["Sansom2010"]]
  trees <- lapply(c(1L, 7L, 99L), FixedTree, dataset = phy)

  lenOff <- withr::with_envvar(c(TS_PACK_LOCAL = "0"),
    vapply(trees, TreeLength, double(1), dataset = phy))
  lenOn <- withr::with_envvar(c(TS_PACK_LOCAL = "1"),
    vapply(trees, TreeLength, double(1), dataset = phy))

  expect_identical(lenOn, lenOff)
})

# --- (2) plane_state[] relabel: reconstruction prints GLOBAL labels -----------
# Fixture forcing a NON-contiguous local alphabet in one block:
#   levels c("-","0","1","2"); applicable labels "0","1","2".
#   char A (has "-") uses "0","1","2" -> all three become global levels and A
#     sits in its own has_inapplicable block.
#   char B (no "-") uses ONLY "1","2" -> its block's local alphabet SKIPS global
#     "0". Packed planes are [plane 0 -> "1", plane 1 -> "2"], so a tip resolved
#     to "2" occupies plane 1. Correct label "2"; the pre-fix code (plane index)
#     printed "1".

test_that("packed reconstruction maps planes back to global state labels", {
  skip_if_not_installed("TreeTools")
  levels <- c("-", "0", "1", "2")
  contrast <- rbind("-" = c(1, 0, 0, 0), "0" = c(0, 1, 0, 0),
                    "1" = c(0, 0, 1, 0), "2" = c(0, 0, 0, 1))
  colnames(contrast) <- levels
  # contrast rows: 1="-", 2="0", 3="1", 4="2".  Columns = char A, char B.
  tipData <- matrix(c(
    1L, 3L,   # t1  A="-"  B="1"
    2L, 3L,   # t2  A="0"  B="1"
    3L, 4L,   # t3  A="1"  B="2"
    3L, 4L,   # t4  A="1"  B="2"
    4L, 3L,   # t5  A="2"  B="1"
    4L, 4L    # t6  A="2"  B="2"
  ), ncol = 2, byrow = TRUE)
  weight <- c(1L, 1L)
  tree <- TreeTools::RenumberTips(
    ape::read.tree(text = "((t1,t2),((t3,t4),(t5,t6)));"), paste0("t", 1:6))
  edge <- TreeTools::Preorder(tree)[["edge"]]

  debugChar <- function(pack, target) {
    withr::with_envvar(c(TS_PACK_LOCAL = if (pack) "1" else "0"),
      TreeSearch:::ts_na_debug_char(edge, contrast, tipData, weight, levels,
                                    target))
  }

  # char B is pattern 2; tips t1,t2,t5 -> "1" and t3,t4,t6 -> "2".
  rOff <- debugChar(FALSE, 2L)
  rOn <- debugChar(TRUE, 2L)
  tipLabOff <- as.character(rOff$prelim)[seq_len(rOff$n_tip)]
  tipLabOn <- as.character(rOn$prelim)[seq_len(rOn$n_tip)]

  # Packing MUST be active for this to test anything: block B collapses from the
  # 3-state global alphabet to its 2-state local alphabet.
  expect_identical(rOff$n_states, 3L)
  expect_identical(rOn$n_states, 2L)

  # The discriminating assertion: the "2"-tips (t3,t4,t6) read "2", NOT the plane
  # index "1". Pre-fix code prints "1" here.
  expect_identical(tipLabOn[c(3L, 4L, 6L)], c("2", "2", "2"))
  expect_identical(tipLabOn[c(1L, 2L, 5L)], c("1", "1", "1"))

  # Tips are resolved singletons, so packed labels equal the unpacked truth.
  expect_identical(tipLabOn, tipLabOff)
})
