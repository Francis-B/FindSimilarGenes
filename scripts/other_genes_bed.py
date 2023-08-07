#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
From the bed file containing all genes, extract lines about genes not in lists
of genes of interest.

"""
import pickle
import re
import copy
import numpy as np

from utils.trie import Trie

# Input and output files
GENES_LISTS = snakemake.input.genes_lists   # List of genes of interest
INBED = snakemake.input.bed     # All genes bed file
OUTBED = snakemake.output.bed
FEATURES = snakemake.output.features    # Dict of transcript features for similar genes search

# # Test input and output files
# GENES_LISTS = ['data/Cancer_Predisposition_Genes.txt', 'data/Cancer_Genes.txt']
# INBED = 'data/all_canonical_transcripts.bed'
# OUTBED = './data/test.bed'
# FEATURES = './data/test.pkl'

genes = []
for GENES in GENES_LISTS:
    with open(GENES, 'r', encoding='utf-8') as f:
        genes.extend([line.rstrip() for line in f])

# Create Trie pattern to fastly retrieve lines with genes of interest
trie = Trie()
trie.add_list(genes)
pattern = re.compile(f'^({trie.pattern()})($|-[0-9])')

# Initiate list to store the canonical transcript of each gene
transcript_features = {}
result_dict = {'length': np.zeros(4, dtype=int), 'number': np.zeros(4, dtype=int)}
features_idx = {'introns': 0,
                'CDSs': 1, 
                'five_prime_UTRs': 2,
                'three_prime_UTRs': 3,
                }

other_genes_lines = []
with open(INBED, 'r', encoding='utf-8') as bed:
    for line in bed:
        name, type_ = line.split('\t')[3].split('::')
        length = int(line.split('\t')[2])+1 - int(line.split('\t')[1])
        # If the line is not about a gene in list
        if not bool(pattern.search(name)):
            other_genes_lines.append(line)
            # Update transcript features
            transcript_features.setdefault(name, copy.deepcopy(result_dict))
            transcript_features[name]['length'][features_idx[f'{type_}s']] += length
            transcript_features[name]['number'][features_idx[f'{type_}s']] += 1


# Write bed file with all longest transcripts
with open(OUTBED, 'w', encoding='utf-8') as f:
    f.writelines(other_genes_lines)

# Dump transcript features to pickle file for further use (similar genes search)
with open(FEATURES, 'wb') as f:
    pickle.dump(transcript_features, f)
