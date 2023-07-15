#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
In the cancer related genes annotation file, find the canonical transcript for
each genes given genes and create a BED file with all the genes' intron,
CDS and UTRs.
"""
import re
import pickle
import numpy as np

from utils.trie import Trie

# Input and output files
GENES = snakemake.input.genes   # List of genes of interest
INBED = snakemake.input.bed     # All genes bed file
OUTBED = snakemake.output.bed   # Genes of interest bed
UNFOUND = snakemake.params.unfound
FEATURES = snakemake.output.features    # Dict of transcript features for similar genes search

# # Test input and output files
# GENES = 'data/Cancer_Genes.txt'
# INBED = 'data/all_canonical_transcripts.bed'
# OUTBED = './data/Cancer_Genes.bed'
# UNFOUND = './data/unfound_genes.txt'
# FEATURES = './data/cancer_genes_FEATURES.pkl'

with open(GENES, 'r', encoding='utf-8') as f:
    genes = [line.rstrip() for line in f]

# Create Trie pattern to fastly retrieve lines with genes of interest
trie = Trie()
trie.add_list(genes)
pattern = re.compile(f'^({trie.pattern()})($|-[0-9])')

# Initiate dict to store the length of each regions of genes-of-interest
transcript_features = {gene: np.zeros(4, dtype=int) for gene in genes}
features_idx = {'introns': 0,
                'CDSs': 1,
                'five_prime_UTRs': 2,
                'three_prime_UTRs': 3,
                }

interest_lines = []    # List to store lines about cancer related genes
with open(INBED, 'r', encoding='utf-8') as bed:
    for line in bed:
        name, type_ = line.split('\t')[3].split('::')
        length = int(line.split('\t')[2]) - int(line.split('\t')[1])
        # If the line is about a gene in list
        if bool(pattern.search(name)):
            interest_lines.append(line)   # Add line to bed file
            # Update transcript features
            transcript_features[name][features_idx[f'{type_}s']] += length


# Write a list of genes with no transcript found in annotation
unfound_genes = [gene for gene, features in transcript_features.items()
                if np.sum(features) == 0]

if bool(unfound_genes):         # Write list to file if there is genes with no transcript
    with open(UNFOUND, 'w', encoding='utf-8') as f:
        f.writelines(gene + '\n' for gene in unfound_genes)

# Write bed file
with open(OUTBED, 'w', encoding='utf-8') as f:
    f.writelines(interest_lines)

# Dump transcript_features to pickle file for further use (similar genes search)
with open(FEATURES, 'wb') as f:
    pickle.dump({gene: features for gene, features in transcript_features.items()
                    if gene not in unfound_genes},
                f)
