# DILS (Day 1)

Hands-on ABC random forests for two-population DILS. Simulations are already done. You will not run `ms`, Snakemake, or the original DILS pipeline.

Slides are the main teaching material. This folder is the working copy for the R scripts.

## How to work

1. Change to this directory (`DILS/`).
2. Start with `scripts/00_check_setup.R`.
3. Open and complete `scripts/01_...` through `scripts/10_...` in order.
4. Do not source or run all scripts at once.
5. Outputs are written to a local `results/` folder (created as needed).
6. Work in your own writable copy of the repository. Do not share a single results directory.

In RStudio you can open `DILS_workshop.Rproj` so the working directory is this folder.

Assignments in R use `=`, not `<-`.

## Packages

Install:

```r
install.packages(c("tidyverse", "abcrf"))
```

`ranger` is pulled in by `abcrf`.

See the repository root `README.md` for tested versus required versions. Exact random-forest results can change with R, `abcrf`, `ranger`, the seed, and parallelism. Workshop forests use 500 trees (`n_trees` in `R/helpers.R`).
