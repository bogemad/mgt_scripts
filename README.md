# Gene Sequence Sorter for Global Alignment
This Python script is designed for bioinformatics applications, specifically to sort homologous genes from bacterial isolates into separate FASTA files. These files can then be used for global alignment using tools like MAFFT or MUSCLE.

## Requirements
Python 3
Biopython
Pandas
Installation
Ensure that Python 3 is installed on your system.
Install the required Python packages if not already installed:

```pip install biopython pandas```

## Usage
The script takes the following command-line arguments:

A FASTA file of genes from a reference genome.
A gene presence-absence file formatted like the output of the Roary pipeline.
A folder containing gene sequence FASTA files (.ffn files from Prokka).
An output folder for the resulting FASTA files.
Run the script using the following command:

```python gene_sequence_sorter.py [reference.fasta] [gene_presence_absence.csv] [prokka_ffn_folder/] [output_folder/]```

## Description

The script begins by reading the reference genome and the gene presence-absence file.
It then iterates through each .ffn file in the specified Prokka folder, gathering and sorting the homologous gene sequences based on the reference.
For each gene, the corresponding sequences from different bacterial isolates are collected and written into separate FASTA files in the output folder.
These output files are ready for global alignment analysis using alignment tools like MAFFT or MUSCLE.

## Output

The script generates a series of FASTA files, each containing sequences of a homologous gene from different bacterial isolates.
The files are named after the reference.fasta gene's locus tag and are saved in the specified output folder.
