# Vienna workshop

This repository contains the material for a two-day workshop:

- `DILS/`: demographic inference with DILS;
- `Aphid/`: gene-flow analyses with Aphid.
- `How_To_Use_LiSC.md`: Instructions for connecting to and using the LiSC server.

## Requirements

The workshop will run on a Linux server with:

- R;
- the R packages `tidyverse` and `abcrf`;
- GCC;
- IQ-TREE 3, available as `iqtree3`;
- ASTRAL, available as `astral`;
- Python 3 with `pandas<=2.1.0`.

Aphid post-processing requires Python 3 and `pandas<=2.1.0`. Any R scripts used for Aphid rely on base R.

## Compiling Aphid

Compile Aphid directly on the teaching server:

```bash
gcc -O2 -std=c11 Aphid/software/aphid.0.11.c \
  -lm \
  -o Aphid/software/aphid
```
