## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----brazeau, eval = FALSE----------------------------------------------------
# MaximizeParsimony(dataset)

## ----hsj, eval = FALSE--------------------------------------------------------
# hierarchy <- CharacterHierarchy("1" = 2:5, "6" = 7:8)
# MaximizeParsimony(dataset, hierarchy = hierarchy,
#                   inapplicable = "hsj", hsj_alpha = 1.0)

## ----xform, eval = FALSE------------------------------------------------------
# hierarchy <- CharacterHierarchy("1" = 2:5)
# MaximizeParsimony(dataset, hierarchy = hierarchy,
#                   inapplicable = "xform")

## ----hierarchy-manual---------------------------------------------------------
library(TreeSearch)

# Character 1 controls characters 2-5
h <- CharacterHierarchy("1" = 2:5)

# Multiple controlling primaries
h <- CharacterHierarchy("1" = 2:5, "6" = 7:8)

## ----hierarchy-auto-----------------------------------------------------------
names <- c("sup_tail", "sub_tail_colour", "sub_tail_shape",
           "sup_wing", "sub_wing_venation", "other_char")
HierarchyFromNames(names)

