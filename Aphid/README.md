# Aphid

Practical analyses use an Isoptera (termite) data set [(Bucek et al., 2019)](https://doi.org/10.1016/j.cub.2019.08.076), with `Empusa_pennata` as the outgroup.

## Layout

```text
Aphid/
├── requirements.txt
├── software/
│   ├── aphid.0.11.c
│   └── aphid                  # created by the compile command below
├── data/
│   ├── alignments/
│   └── third_posi_codon/
├── trees_topo/
├── trees/
├── isoptera_genetrees.treefile
├── inputs_aphid/
├── outputs/
└── scripts/
    └── process_aphid.py
    └── plot_contributions.py
```

All commands below assume the working directory is `Aphid/`.

## Compile Aphid

```bash
gcc -O2 -std=c11 software/aphid.0.11.c \
  -lm \
  -o software/aphid
```

## Gene-tree topology

From `data/`. Replace `<locus>` by the locus name without the extension.

```bash
cd data

iqtree3 \
  -s ./alignments/<locus>.fas \
  -m MFP \
  -o Empusa_pennata \
  -pre ../trees_topo/<locus> \
  -T AUTO \
  --threads-max 2
```

## Branch lengths on third codon positions

Still from `data/`:

```bash
iqtree3 \
  -s ./third_posi_codon/<locus>.fas \
  -te ../trees_topo/<locus>.treefile \
  -m MFP \
  -o Empusa_pennata \
  -pre ../trees/<locus> \
  -T AUTO \
  --threads-max 2
```
The scripts `1_run_iqtree_topology.sh` and `2_run_iqtree_branch_length.sh` can be used to automatically generate the gene trees:

```bash
bash scripts/1_run_iqtree_topology.sh # Used to get the topology
bash scripts/2_run_iqtree_branch_length.sh # Used to recalculate the branch lengths
```
Because generating a tree for each gene is time-consuming, precomputed trees are provided:
```bash
tar -xzf ./trees.tar.gz
```

## Concatenation

From `Aphid/`:

```bash
cat ./trees/*.fas.tree > isoptera_genetrees.treefile
```

## ASTRAL

```bash
astral \
  -i isoptera_genetrees.treefile \
  -o astral_isoptera_speciestree.tree
```

## Create Aphid input file
```bash
python ./scripts/aphid_infile.py --trees ./trees --alignments ./data/third_posi_codon --extension "fas" --output ./inputs_aphid/isoptera.in
```

## Aphid (verbose)

First triplet:

```bash
./software/aphid \
  inputs_aphid/isoptera.in \
  inputs_aphid/isoptera_1.tax \
  inputs_aphid/config_v.opt \
  outputs/isoptera_1.csv
```

Second triplet:

```bash
./software/aphid \
  inputs_aphid/isoptera.in \
  inputs_aphid/isoptera_2.tax \
  inputs_aphid/config_v.opt \
  outputs/isoptera_2.csv
```

## Standard (non-verbose) tables

`config_v.opt`(`verbose = 1`) outputs a human-readable version of the Aphid results. However, this output is difficult to parse automatically. To address this, use `config.opt`(`verbose = 0`), which prints a single summary line to stdout. This line can be redirected to a file for downstream analyses, such as calculating the median timing of GF and the contribution to phylogenetic conflict for each species pair.

First triplet:

```bash
echo "dataset,option" > prov1
echo "nb_gene,ntopo0,ntopo1,ntopo2,ntopo3,mean_lg,tau1,tau1_l,tau1_h,tau2,tau2_l,tau2_h,theta,theta_l,theta_h,pab,pab_l,pab_h,pac,pac_l,pac_h,pbc,pbc_l,pbc_h,pa,pa_l,pa_h,noevent,noconflict_ILS,noconflict_GF,conflict_ILS,conflict_ILS_l,conflict_ILS_h,conflict_GF,conflict_GF_l,conflict_GF_h,imbalance_ILS,dominant_ILS,imbalance_GF,dominant_GF,lnL" > prov2

echo ./inputs_aphid/isoptera_1.tax,basic.opt >> prov1
./software/aphid inputs_aphid/isoptera.in inputs_aphid/isoptera_1.tax inputs_aphid/config.opt outputs/isoptera_1.csv >> prov2

paste -d"," prov1 prov2 > ./isoptera_1.csv

rm prov1
rm prov2
```

Second triplet:

```bash
echo "dataset,option" > prov1
echo "nb_gene,ntopo0,ntopo1,ntopo2,ntopo3,mean_lg,tau1,tau1_l,tau1_h,tau2,tau2_l,tau2_h,theta,theta_l,theta_h,pab,pab_l,pab_h,pac,pac_l,pac_h,pbc,pbc_l,pbc_h,pa,pa_l,pa_h,noevent,noconflict_ILS,noconflict_GF,conflict_ILS,conflict_ILS_l,conflict_ILS_h,conflict_GF,conflict_GF_l,conflict_GF_h,imbalance_ILS,dominant_ILS,imbalance_GF,dominant_GF,lnL" > prov2

echo ./inputs_aphid/isoptera_2.tax,basic.opt >> prov1
./software/aphid inputs_aphid/isoptera.in inputs_aphid/isoptera_2.tax inputs_aphid/config.opt outputs/isoptera_2.csv >> prov2

paste -d"," prov1 prov2 > ./isoptera_2.csv

rm prov1
rm prov2
```

## Python post-processing

```bash
python3 -m pip install -r requirements.txt
```

First triplet:

```bash
python3 scripts/process_aphid.py \
  -t Isoptera \
  --aphid_standard ./isoptera_1.csv \
  --aphid_output ./outputs/isoptera_1.csv \
  --output ./isoptera_processed_1.csv \
  --times 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1
```

Second triplet:

```bash
python3 scripts/process_aphid.py \
  -t Isoptera \
  --aphid_standard ./isoptera_2.csv \
  --aphid_output ./outputs/isoptera_2.csv \
  --output ./isoptera_processed_2.csv \
  --times 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1
```

## Plotting the contribution

First triplet:

```bash
Rscript scripts/plot_results.R -r ./isoptera_processed_1.csv -p ./isoptera_1_plot.pdf
```

Second triplet:

```bash
Rscript scripts/plot_results.R -r ./isoptera_processed_2.csv -p ./isoptera_2_plot.pdf
```