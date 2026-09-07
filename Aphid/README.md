# Aphid (Day 2)

Day 2 of the Vienna workshop, led by Arthur Boddaert.

This folder currently contains the Aphid C program and empty destinations for data, trees, inputs, outputs, and R scripts. **The Day 2 pipeline is not reproducible from this repository yet.** Alignments, trees, Aphid input files, and workshop R scripts are still expected from Arthur.

## Software

| File | Role |
| --- | --- |
| `software/aphid.0.11.c` | Upstream source (do not edit to silence compiler warnings) |
| `software/LICENSE` | GPL-3.0 from the upstream project |
| `software/README_upstream.md` | Original Aphid README |
| `software/containerize.bash` | Upstream container helper |
| `software/aphid_precompiled_linux_x86_64_glibc_2.34` | ELF Linux x86-64 binary linked against **GLIBC_2.34**. Not a universal binary. |

SHA256 (verified on copy):

```
0511e5adea40a1d81ccc913c1a9c39122d80ee84aa255b186b934a0c9eca51f5  software/aphid_precompiled_linux_x86_64_glibc_2.34
702dbafec1b97d617b2e2cded08cf4285250dc7a32fbfff64c9afca63a4a8b02  software/aphid.0.11.c
```

Upstream clone (kept locally, not in this public repository): https://gitlab.com/iago-lito/aphid.git (tag `v0.11`).

## Compile on the server

Do **not** treat the precompiled binary as portable. Compile from source:

```bash
gcc -O2 -std=c11 -Wall -Wextra \
  Aphid/software/aphid.0.11.c -lm \
  -o Aphid/software/aphid
```

`-lm` must come **after** the source (or object) file. Compilation emits several warnings; they are expected. Do not modify Arthur's scientific C code to remove them.

The original precompiled file must stay untouched. Write the server build to `Aphid/software/aphid` (a different name).

With no arguments the program prints:

```text
usage: aphid tree_file taxon_file option_file outfile
```

It expects four arguments.

## Planned pipeline (not runnable yet)

Commands below are the intended Day 2 workflow. They will fail until Arthur's files are added.

### Gene-tree topology

```bash
cd Aphid/data

iqtree3 \
  -s ./alignments/<locus>.fas \
  -m MFP \
  -o Labidesthes_sicculus \
  -pre ../trees_topo/<locus> \
  -T AUTO \
  --threads-max 2
```

### Branch lengths on third codon positions

```bash
iqtree3 \
  -s ./third_posi_codon/<locus>.fas \
  -te ../trees_topo/<locus>.treefile \
  -m MFP \
  -o Labidesthes_sicculus \
  -pre ../final_trees/<locus> \
  -T AUTO \
  --threads-max 2
```

IQ-TREE may use up to two threads with these commands.

### Concatenation

```bash
cat ./final_trees/*.treefile > cichlids_concated.treefile
```

### ASTRAL

```bash
astral \
  -i cichlids_concated.treefile \
  -o cichlids_speciestree_astral.treefile
```

**Exact ASTRAL implementation and version: pending confirmation from Arthur Boddaert.** Do not assume ASTRAL-III, ASTRAL-IV, or ASTER. If the chosen tool needs Java, that dependency will be added after confirmation.

### Aphid

```bash
Aphid/software/aphid \
  ./inputs_aphid/cichlids.in \
  ./inputs_aphid/cichlids.tax \
  ./inputs_aphid/config_v.opt \
  ./outputs/gene_cichlids.cichlids
```

## Pending material from Arthur Boddaert

- Alignments (`data/alignments/`)
- Third-position alignments (`data/third_posi_codon/`)
- Optional precomputed trees
- `inputs_aphid/` files `.in`, `.tax`, and `.opt`
- Workshop R scripts under `scripts/`
- Exact ASTRAL implementation and version
- Final Day 2 instructions

## R

Future Aphid R scripts are assumed to use **base R only**. No extra R packages are required for Day 2 at this time.
