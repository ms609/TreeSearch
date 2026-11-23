# R Package Conventions

Load this when: adding `.R` files, writing roxygen docs, updating vignettes,
or running pre-commit documentation checks.

---

## R source file ordering

`DESCRIPTION` has an explicit `Collate:` field. When adding a new `.R` file,
**update the Collate field** — otherwise R sources alphabetically, which can
break if one file's top-level code depends on a later file.

---

## Documentation checks (mandatory)

After any change to a function signature or roxygen block, run:

```r
devtools::check_man()
```

After writing or updating documentation prose, also run:

```r
spelling::spell_check_package()
```

Both should be clean before committing. `check_man` catches Rd parse errors,
cross-ref failures, `\usage` mismatches; `spell_check_package` catches typos
in `@description`/`@details`/`@param` text.

References are added using Rdpack's `\insertCite{}`, with
`\insertAllCited{}` in the references section.

---

## Algorithm vignette (mandatory updates)

`vignettes/search-algorithm.Rmd` documents the search algorithm for
publication. **Any change that modifies search behaviour** — new heuristics,
parameter tuning, scoring methods, stopping criteria, pool management, or
rearrangement operators — **must be accompanied by an update to this vignette.**

- Published techniques: add a short summary and `@Key` citation.
- Novel contributions: describe the algorithm in enough detail for a reader
  to understand the design and rationale. Include empirical results where
  available (e.g. benchmark deltas).
- New references: add `@article{Key, ...}` to `inst/REFERENCES.bib`.

The vignette uses pandoc-style `@Key` citations (same as the other
vignettes), not Rdpack `\insertCite{}`.
