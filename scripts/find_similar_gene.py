#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Find the most similar not cancer related gene for each cancer related gene. 
"""

import re
import pickle
import numpy as np

from utils.trie import Trie

# Input and output files
BED = snakemake.input.other_bed   # Other genes bed file
INTEREST_FEATURES = snakemake.input.interest_features
OTHER_FEATURES = snakemake.input.other_features
SIMILAR_BED = snakemake.output.similar
GENE_PAIRS = snakemake.output.gene_pairs

# # Test input and output files
# BED = 'data/other_genes.bed'
# INTEREST_FEATURES = 'data/Cancer_Genes_FEATURES.pkl'
# OTHER_FEATURES = 'data/other_genes_FEATURES.pkl'
# SIMILAR_BED = 'data/Cancer_Genes_similar.bed'
# GENE_PAIRS = 'data/Cancer_Genes_Best_Matches.pkl'

with open(INTEREST_FEATURES, 'rb') as f:
    cancer_features_dict = pickle.load(f)
with open(OTHER_FEATURES, 'rb') as f:
    not_cancer_features = pickle.load(f)

null_genes = np.array(list(not_cancer_features.keys()))   # genes name
null_features = np.array(list(not_cancer_features.values()))    # features

best_matches = dict.fromkeys(cancer_features_dict.keys())

print('Finding similar genes... This may take few minutes')
# Loop through cancer related genes to find similar not cancer related genes
for gene, features in cancer_features_dict.items():
    # Create mask to exclude not-cancer-related genes with regions of length 0 
    # which are not null in cancer related genes
    non_zero_regions = np.nonzero(features)   # Index of null length regions
    mask = [np.all(arr[non_zero_regions]) for arr in null_features]
    features_to_compare = iter(null_features[mask])  # Features of genes with non zero regions
    # Calculate the difference between genes
    diff = [np.sum(np.abs(features - arr)) for arr in features_to_compare]
    best_matches[gene] = null_genes[mask][np.argmin(diff)]  # Gene which yields the smallest difference

# Initialize Trie to speed up genes search
trie = Trie()
trie.add_list(list(best_matches.values()))
pattern = re.compile(f'^({trie.pattern()})($|-[0-9])')

# Loop through all not-cancer-related genes and extract all best matches
similar_genes = []
with open(BED, 'r', encoding='utf-8') as bed:
    for line in bed:
        if bool(pattern.search(line.split('\t')[3].split('::')[0])):
            similar_genes.append(line)

# Write bed file with all similar genes
with open(SIMILAR_BED, 'w', encoding='utf-8') as f:
    f.writelines(similar_genes)

# Save best_matches dict to pickle file
with open(GENE_PAIRS, 'wb') as f:
    pickle.dump(best_matches, f)
