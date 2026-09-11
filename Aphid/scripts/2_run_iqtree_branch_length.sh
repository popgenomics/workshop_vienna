#!/bin/bash

for fasta in ./data/alignments/*.fas; do
	name=$(basename "$fasta" .fas)
	iqtree3 -s ./data/third_posi_codon/"$name".fas -te ./trees_topo/"$name".fas.treefile -m MFP -o Empusa_pennata -pre ./trees/"$name" -T AUTO --threads-max 2
done
