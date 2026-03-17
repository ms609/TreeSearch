# Extracted from test-length_range.R:2

# test -------------------------------------------------------------------------
expect_equal(MinimumLength(1:3), expect_warning(MinimumSteps(1:3)))
