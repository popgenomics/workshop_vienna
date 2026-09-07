# Vienna workshop

This repository contains the material for a two-day workshop:

- `DILS/`: demographic inference with DILS;
- `Aphid/`: gene-flow analyses with Aphid.

## Requirements

The workshop will run on a Linux server with:

- R;
- the R packages `tidyverse` and `abcrf`;
- GCC;
- IQ-TREE 3, available as `iqtree3`;
- ASTRAL, available as `astral`.

The Aphid R scripts use base R and do not require additional R packages.

## Compiling Aphid

Compile Aphid directly on the teaching server:

```bash
gcc -O2 -std=c11 Aphid/software/aphid.0.11.c \
  -lm \
  -o Aphid/software/aphid
```
