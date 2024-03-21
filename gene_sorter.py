#!/usr/bin/env python
 
import sys, os, csv
from Bio import SeqIO
from collections import defaultdict
import pandas as pd
import numpy as np
 
if len(sys.argv) != 5:
    print("Usage: script.py <core_gene_ref> <gene_pres_abs> <genomes_for_dnds_dir> <output_folder>")
    sys.exit(1)
 
# Define inputs and outputs from command line arguments
core_gene_ref = open(sys.argv[1])  # Reference genome fasta file
gene_pres_abs = open(sys.argv[2])  # Gene presence-absence file (Roary output)
genomes_for_dnds = os.listdir(sys.argv[3])  # Directory containing .ffn files from Prokka
output_folder = sys.argv[4]  # Directory for output fasta files
 
            # Testing
            # core_gene_ref="/home/wordenp/projects/MGT_project/analyses/Output/AFB_50_Representative_Genomes_for_MGT/Prokka_FFN_Only_AFB_50_Reps/prokka_GCF_002951875.1_ASM295187v1_ERIC-I_genomic_subset.ffn"
            # gene_pres_abs="/home/wordenp/projects/MGT_project/analyses/Output/AFB_50_Representative_Genomes_for_MGT/Roary_AFB_50_Reps/subsetForMGT_gene_presAbs.csv"
            # genomes_for_dnds="/home/wordenp/projects/MGT_project/analyses/Output/AFB_50_Representative_Genomes_for_MGT/Prokka_FFN_Only_AFB_50_Reps"
            # output_folder="/home/wordenp/projects/MGT_project/analyses/Output/AFB_50_Representative_Genomes_for_MGT/1_eachGeneInAllGenomes_FASTA"
 
os.makedirs(output_folder, exist_ok=True)  # Create output folder if it doesn't exist
ref = os.path.basename(sys.argv[1].replace('.ffn', ''))  # Reference genome name
 
# Import gene presence-absence data into a pandas DataFrame
df = pd.read_csv(gene_pres_abs, low_memory=False)
 
# Import reference genome and extract locus tags
cgd = SeqIO.to_dict(SeqIO.parse(core_gene_ref, "fasta"))
locus_tags = list(cgd)
 
# Initialize containers for genome data
genome_locus_tags = defaultdict(dict)
genomes = {}
genomes[ref] = SeqIO.to_dict(SeqIO.parse(sys.argv[1], 'fasta'))
 
# Read and store locus tags and sequences for each genome
for sequence in genomes_for_dnds:
    # Filter and add locus tags for each genome based on reference
    sequence_path = os.path.join(sys.argv[3], sequence)
    genomes[sequence.replace('.ffn', '')] = SeqIO.to_dict(SeqIO.parse(sequence_path, 'fasta'))
    for locus_tag in locus_tags:
        genome_locus_tags[locus_tag][sequence] = df.query(f"{ref}=={locus_tag}")[sequence.replace('.ffn', '')]
    #genome_locus_tags.append(df.loc[df[ref].isin(locus_tags), sequence.replace('.ffn', '')].values.tolist())
 
 
print("Locus tags/isolates with split genes:")
 
# Transpose the list of lists to align locus tags properly with their respective genomes as columns
#genome_locus_tags_transposed = list(map(list, zip(*genome_locus_tags)))
 
# Create a DataFrame from the transposed data
#genome_locus_tags_df = pd.DataFrame(genome_locus_tags_transposed)
#genome_locus_tags_df.columns = genomes_for_dnds
 
# Process each locus tag
locus_tags_not_found = []
for locus_tag in locus_tags:
    outputfile = os.path.join(output_folder, locus_tag + '.fa')
    locus_tag_seqs = []
    # Gather sequences corresponding to each locus tag from different genomes
    for sequence in genomes_for_dnds:
        # print(genomes[sequence.replace('.ffn', '')])
        try:
            if ';' in genome_locus_tags[locus_tag][sequence] or "\t" in genome_locus_tags[locus_tag][sequence]:
                print(f"{locus_tag},{sequence.replace('.ffn', '')},({genome_locus_tags[locus_tag][sequence]})")
                continue
        #except TypeError:
            #continue
        except KeyError:
            locus_tags_not_found.append(locus_tag)
        locus_tag_seqs.append(genomes[sequence.replace('.ffn', '')][genome_locus_tags[locus_tag][sequence]])
    # Write the gathered sequences to a fasta file
    SeqIO.write(locus_tag_seqs, outputfile, 'fasta')
 
outputfile = os.path.join(output_folder, 'locus_tags_not_found.txt')
with open(outputfile, 'w') as outh:
    outh.write("\n".join(locus_tags_not_found) + '\n')
