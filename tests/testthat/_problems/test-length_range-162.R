# Extracted from test-length_range.R:162

# test -------------------------------------------------------------------------
manyStates <- c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512)
expect_silent(result <- MaximumLength.numeric(manyStates))
