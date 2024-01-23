#!/usr/bin/env python

import sys, os, csv
from Bio import SeqIO
from collections import defaultdict
import pandas as pd

# Define inputs and outputs from command line arguments
core_gene_ref = open(sys.argv[1])  # Reference genome fasta file
gene_pres_abs = open(sys.argv[2])  # Gene presence-absence file (Roary output)
genomes_for_dnds = os.listdir(sys.argv[3])  # Directory containing .ffn files from Prokka
output_folder = sys.argv[4]  # Directory for output fasta files
os.makedirs(output_folder, exist_ok=True)  # Create output folder if it doesn't exist
ref = sys.argv[1].replace('.ffn', '')  # Reference genome name

# Import gene presence-absence data into a pandas DataFrame
df = pd.read_csv(gene_pres_abs)

# Import reference genome and extract locus tags
cgd = SeqIO.to_dict(SeqIO.parse(core_gene_ref, "fasta"))
locus_tags = list(cgd)

# Initialize containers for genome data
genome_locus_tags = []
genomes = {}
genomes[ref] = SeqIO.to_dict(SeqIO.parse(sys.argv[1], 'fasta'))

# Read and store locus tags and sequences for each genome
for sequence in genomes_for_dnds:
    # Filter and add locus tags for each genome based on reference
    genome_locus_tags.append(df.loc[df[ref].isin(locus_tags), sequence.replace('.ffn', '')].values.tolist())
    sequence_path = os.path.join(sys.argv[3], sequence)
    genomes[sequence.replace('.ffn', '')] = SeqIO.to_dict(SeqIO.parse(sequence_path, 'fasta'))

print("Locus tags/isolates with split genes:")

# Process each locus tag
for i, locus_tag in enumerate(locus_tags):
    outputfile = os.path.join(output_folder, locus_tag + '.fa')
    locus_tag_seqs = []

    # Gather sequences corresponding to each locus tag from different genomes
    for j, sequence in enumerate(genomes_for_dnds):
        if ';' in genome_locus_tags[j][i]:
            print(f"{locus_tag},{sequence.replace('.ffn', '')},({genome_locus_tags[j][i]})")
            continue
        locus_tag_seqs.append(genomes[sequence.replace('.ffn', '')][genome_locus_tags[j][i]])

    # Write the gathered sequences to a fasta file
    SeqIO.write(locus_tag_seqs, outputfile, 'fasta')
