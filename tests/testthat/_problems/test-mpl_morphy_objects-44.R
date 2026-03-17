# Extracted from test-mpl_morphy_objects.R:44

# test -------------------------------------------------------------------------
morphyObj <- SingleCharMorphy("1")
on.exit(UnloadMorphy(morphyObj))
expect_error(morphy_profile(matrix(NA, 10, 2), list(morphyObj),
                              1, 1L, matrix(1), 1),
               "Number of edges does not match Morphy object dimensions")
