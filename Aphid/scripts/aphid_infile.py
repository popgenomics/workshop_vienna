import os
import argparse
from Bio import SeqIO

def read_fasta(alignment):
    with open(alignment) as filein:
        filein.readline()
        seq = filein.readline()
        seq = seq.rstrip()
    return len(seq)

def tree2aphid(trees, alignments, extension, output):
    input_string = ''
    for file in os.listdir(trees):
        if not file.startswith('.'):
            name = file.split('.')[0]
            alignment = os.path.join(alignments, f'{name}.{extension}')
            tree = os.path.join(trees, file)
            length = read_fasta(alignment=alignment)
            with open(tree) as filein:
                newick_tree = filein.readline()
                newick_tree = newick_tree.rstrip()
            input_string += f'{newick_tree}\t{length}\t{name}\n'
        with open (output, 'w') as fileout:
            fileout.write(input_string)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--trees')
    parser.add_argument('-a', '--alignments')
    parser.add_argument('-e', '--extension')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()

    tree2aphid(trees=args.trees, alignments=args.alignments, extension=args.extension, output=args.output)