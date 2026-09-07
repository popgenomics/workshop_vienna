# Server checks

Non-destructive checks for the Linux workshop server.

From the repository root:

```bash
bash server/check_prerequisites.sh
```

The shell script:

- does not write into the repository
- compiles Aphid in a temporary directory (`mktemp -d`) and deletes it on exit
- does not run IQ-TREE, ASTRAL, or a full DILS analysis

`server/check_R_packages.R` uses base R only. It reports R and package versions, compares them to the versions tested on the preparation machine, and checks that the five DILS teaching tables exist.

## What is blocking vs a warning

Blocking (non-zero exit) if missing: `bash`, `git`, `R`, `Rscript`, `gcc`, DILS packages (`tidyverse`, `abcrf`, `ranger`), DILS teaching tables, a writable temp directory, or a failed Aphid test compile.

Warnings (do not by themselves mean the DILS day is broken): R or package version differences other than `abcrf`; missing `iqtree3` or `astral`; unconfirmed ASTRAL implementation. **Exact ASTRAL implementation and version: pending confirmation from Arthur Boddaert.** Presence of a command named `astral` is not a scientific validation of the Day 2 pipeline.

`abcrf` **1.9** is required. Do not silently replace it with a newer version: random-forest results can change.

## Student copies

Each student needs a writable clone or copy. `DILS/scripts/` write under `DILS/results/`. Aphid will write under `Aphid/trees_topo/`, `Aphid/trees/`, and `Aphid/outputs/`.
