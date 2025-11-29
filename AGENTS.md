# Repository Guidelines

## Project Structure & Module Organization
- `R/` – exported R functions for spatial/feature adjacency, kernels, and helpers (e.g., `neighbor_graph.R`, `knn_weights.R`).
- `src/` – Rcpp/Armadillo implementations (e.g., `weight_funs.cpp`); compiled objects may appear locally but should not be committed.
- `tests/testthat/` – unit tests grouped by feature (`test-spatial-weights.R`, `test-commute-time.R`, etc.).
- `data-raw/` – scripts that generate sample datasets; outputs land in `test_data/` or `data/`.
- `docs/` + `_pkgdown.yml` – pkgdown site config; `man/` is roxygen output; `README.Rmd` is the source for `README.md`.

## Build, Test, and Development Commands
- `R -q -e "devtools::document()"` – regenerate roxygen docs and `NAMESPACE` after API changes.
- `R -q -e "devtools::test()"` – run the `testthat` suite.
- `R CMD build .` – create a source tarball.
- `R CMD check neighborweights_*.tar.gz --as-cran` – full CRAN-style checks (slow but authoritative).
- `R -q -e "pkgdown::build_site()"` – rebuild docs site locally.
- `R -q -e "devtools::install()"` – install the package for interactive use.

## Coding Style & Naming Conventions
- Base R style with 2-space indents, no tabs; keep lines ≲100 characters.
- Functions use `snake_case`; S3 methods documented with `@method` and `@export`.
- Space around `=` in arguments (`sigma = 1`); prefer explicit defaults.
- Keep roxygen examples minimal and runnable; tag helpers with `@keywords internal`.
- C++: follow existing Rcpp style, prefer `const` references, and add brief `//` comments for non-obvious logic.

## Testing Guidelines
- Place tests in `tests/testthat/` as `test-<topic>.R`; mirror file names in `R/` when feasible.
- Use `testthat` v3 expectations (`expect_equal(..., tolerance = 1e-8)` for numeric results).
- Reuse fixtures in `test_data/`; keep random seeds fixed when generating data.
- Run `devtools::test()` before pushing; run `R CMD check` before tagging a release.

## Commit & Pull Request Guidelines
- Commit messages: short, imperative summaries (e.g., “Fix example errors in R CMD check”).
- PRs should describe intent, key changes, and verification commands; link issues when available.
- Note any API changes and regenerated artifacts (`man/`, `NAMESPACE`, pkgdown); avoid committing local `.o`/`.so` builds.

## Security & Configuration Tips
- Keep `DESCRIPTION`/`NAMESPACE` in sync with code changes; update `Imports`/`LinkingTo` when adding new dependencies.
- If touching Rcpp, run `Rcpp::compileAttributes()` before documenting; clean stray binaries with `devtools::clean_dll()` when needed.
