# Aphid

Practical analyses use an Isoptera (termite) data set, with `Empusa_pennata` as the outgroup.

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
├── inputs_aphid/
├── outputs/
└── scripts/
    └── process_aphid.py
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

## Concatenation

From `Aphid/`:

```bash
cat ./trees/*.treefile > isoptera_genetrees.treefile
```

## ASTRAL

```bash
astral \
  -i isoptera_genetrees.treefile \
  -o astral_isoptera_speciestree.tree
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

`config.opt` (`verbose = 0`) prints one summary line on stdout. The `option` field `basic.opt` is only a label in that table; `process_aphid.py` does not read it. The file actually used for the run is `inputs_aphid/config.opt`.

First triplet:

```bash
echo "dataset,option,nb_gene,ntopo0,ntopo1,ntopo2,ntopo3,av_lg,tau1,tau1_low,tau1_high,tau2,tau2_low,tau2_high,theta,theta_low,theta_high,pab,pab_low,pab_high,pac,pac_low,pac_high,pbc,pbc_low,pbc_high,pa,pa_low,pa_high,no_event,noconflict_ILS,noconflict_GF,conflict_ILS,conflict_ILS_low,conflict_ILS_high,conflict_GF,conflict_GF_low,conflict_GF_high,imbalance_ILS,dominant_ILS,imbalance_GF,dominant_GF,max_lnL" > isoptera_1.csv

./software/aphid \
  inputs_aphid/isoptera.in \
  inputs_aphid/isoptera_1.tax \
  inputs_aphid/config.opt \
  prov1 \
  | awk '{print "inputs_aphid/isoptera_1.tax,basic.opt," $0}' >> isoptera_1.csv

rm -f prov1
```

Second triplet:

```bash
echo "dataset,option,nb_gene,ntopo0,ntopo1,ntopo2,ntopo3,av_lg,tau1,tau1_low,tau1_high,tau2,tau2_low,tau2_high,theta,theta_low,theta_high,pab,pab_low,pab_high,pac,pac_low,pac_high,pbc,pbc_low,pbc_high,pa,pa_low,pa_high,no_event,noconflict_ILS,noconflict_GF,conflict_ILS,conflict_ILS_low,conflict_ILS_high,conflict_GF,conflict_GF_low,conflict_GF_high,imbalance_ILS,dominant_ILS,imbalance_GF,dominant_GF,max_lnL" > isoptera_2.csv

./software/aphid \
  inputs_aphid/isoptera.in \
  inputs_aphid/isoptera_2.tax \
  inputs_aphid/config.opt \
  prov2 \
  | awk '{print "inputs_aphid/isoptera_2.tax,basic.opt," $0}' >> isoptera_2.csv

rm -f prov2
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
