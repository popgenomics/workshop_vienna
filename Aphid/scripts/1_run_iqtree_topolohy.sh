#!/bin/bash

for fasta in ./data/alignments/*.fas; do
	name=$(basename "$fasta" .cds)
	iqtree3 -s "$fasta" -m MFP -o Empusa_pennata -pre ./trees_topo/"$name" -T AUTO --threads-max 2
done
