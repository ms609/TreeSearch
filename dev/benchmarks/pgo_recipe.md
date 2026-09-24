# PGO (Profile-Guided Optimization) Build Recipe

## Overview

PGO lets GCC optimize branch prediction, function layout, and inlining
decisions based on actual runtime behavior. Requires two compilation
passes: one instrumented build to gather profile data, then a second
build that uses that data for optimization.

## Results (2026-03-16, GCC 13 / rtools45, Windows x86_64)

| Benchmark | Baseline (s) | PGO (s) | Speedup |
|-----------|-------------|---------|---------|
| Vinther EW (23 tips) | 0.240 | 0.240 | 0% |
| Vinther IW (23 tips) | 0.170 | 0.190 | -12% |
| Zhu EW (75 tips) | 4.010 | 3.790 | 5% |
| Zhu IW (75 tips) | 5.340 | 4.990 | 7% |
| Agnarsson EW (62 tips) | 2.200 | 2.080 | 5% |

PGO provides a modest ~5-7% speedup on medium-sized datasets where the
C++ hot path dominates. On small datasets, R overhead and startup time
swamp any C++ improvement. Scores are identical (correctness verified:
53/53 driven search tests pass).

## Build Steps

All steps run from the package root directory.

### Step 1: Baseline build (no PGO)

Ensure no `src/Makevars.win` exists:

```bash
rm -f src/Makevars.win src/*.o src/*.dll
R CMD INSTALL --library=.agent-pgo .
```

### Step 2: Instrumented build

Create `src/Makevars.win`:

```makefile
PROFILE_DIR = C:/Users/pjjg18/GitHub/TreeSearch/.pgo-data
PKG_CXXFLAGS = -fprofile-generate=$(PROFILE_DIR)
PKG_CFLAGS = -fprofile-generate=$(PROFILE_DIR)
PKG_LIBS = -fprofile-generate
```

Build and install:

```bash
rm -rf .pgo-data && mkdir .pgo-data
rm -f src/*.o src/*.dll
R CMD INSTALL --library=.agent-pgo-gen .
```

### Step 3: Training workload

Load the instrumented build and exercise all major code paths:

```r
library(TreeSearch, lib.loc = ".agent-pgo-gen")
data(inapplicable.phyData, package = "TreeSearch")

# EW + IW on small and medium datasets
MaximizeParsimony(inapplicable.phyData[["Vinther2008"]],
                  maxReplicates = 5L, targetHits = 3L, verbosity = 0L)
MaximizeParsimony(inapplicable.phyData[["Vinther2008"]], concavity = 10,
                  maxReplicates = 5L, targetHits = 3L, verbosity = 0L)
MaximizeParsimony(inapplicable.phyData[["Zhu2013"]],
                  maxReplicates = 3L, targetHits = 2L, verbosity = 0L)
MaximizeParsimony(inapplicable.phyData[["Zhu2013"]], concavity = 10,
                  maxReplicates = 3L, targetHits = 2L, verbosity = 0L)
MaximizeParsimony(inapplicable.phyData[["Agnarsson2004"]],
                  maxReplicates = 3L, targetHits = 2L, verbosity = 0L)
```

The `.gcda` files appear under `.pgo-data/C~/Users/.../src/`.

### Step 4: PGO-use build

Replace `src/Makevars.win`:

```makefile
PROFILE_DIR = C:/Users/pjjg18/GitHub/TreeSearch/.pgo-data
PKG_CXXFLAGS = -fprofile-use=$(PROFILE_DIR) -fprofile-correction
PKG_CFLAGS = -fprofile-use=$(PROFILE_DIR) -fprofile-correction
PKG_LIBS = -fprofile-use
```

Build (**note: takes 3-5 minutes**, much longer than normal):

```bash
rm -f src/*.o src/*.dll
R CMD INSTALL --library=.agent-pgo-use .
```

### Step 5: Clean up

**Always remove `src/Makevars.win` after PGO builds** — leaving PGO
flags in place will cause segfaults (instrumented build) or broken
builds (PGO-use without matching `.gcda` files):

```bash
rm -f src/Makevars.win src/*.o src/*.dll
```

## Notes

- `-fprofile-correction` is needed because some source files may have
  changed since profile generation. It tells GCC to accept mismatched
  profiles gracefully rather than erroring.
- The `.pgo-data/` directory contains machine-specific binary data.
  Do not commit to version control.
- PGO-use compilation is 2-5× slower than normal. Allow 5 minutes for
  a full rebuild (30+ source files).
- GCC on Windows (rtools45) nests `.gcda` files under a path encoding
  like `C~/Users/...`. This is expected behavior.
