#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
Create an Excel file to compare the length and number of each regions for all
the similar genes pairs.
"""

import os
import pickle
import json
import pandas as pd

INTEREST_FEATURES = snakemake.input.interest_features
OTHER_FEATURES = snakemake.input.other_features
GENE_PAIRS = snakemake.input.gene_pairs
EXCEL = snakemake.output.excel

features_idx = {'introns': 0,
                'CDSs': 1,
                'five_prime_UTRs': 2,
                'three_prime_UTRs': 3,
                }

with open(OTHER_FEATURES, 'rb') as f:
    other_features = pickle.load(f)

with pd.ExcelWriter(EXCEL) as writer: # pylint: disable=abstract-class-instantiated
    for interest_feature_file, gene_pairs_file in zip(INTEREST_FEATURES, GENE_PAIRS):
        with open(interest_feature_file, 'rb') as f:
            interest_features = pickle.load(f)
        with open(gene_pairs_file, 'rb') as f:
            gene_pairs = json.load(f)

        index = pd.MultiIndex.from_product([['Cancer_Genes', 'Similar_Genes'],
                                            ['intron', 'five_prime_UTR', 'CDS', 'three_prime_UTR'],
                                            ['length', 'number']])

        # Master dict containing counts for both cancer genes and similar genes
        master_dict = {'Cancer_Genes': interest_features, 'Similar_Genes': other_features}
        gene_lists = {'Cancer_Genes': gene_pairs.keys(), 'Similar_Genes': gene_pairs.values()}
        data = {}
        # For each multiindex tuple, collect corresponding data in master_dict
        for distribution, region, stat in index:
            # Append new list (column) to data
            dict_ = master_dict[distribution]   # Get corresponding distribution's dict
            data[(distribution, region, stat)] = [dict_[gene][stat][features_idx[f'{region}s']]
                                                for gene in gene_lists[distribution]]

        # Add gene names to data
        new_index = [('Cancer_Genes', '', 'Names')] + \
                    index.to_list() + \
                    [('Similar_Genes', '', 'Names')]

        for distribution in master_dict:
            data[(distribution, '', 'Names')] = gene_lists[distribution]

        df = pd.DataFrame(data, columns=pd.MultiIndex.from_tuples(new_index))

        sheet_name = os.path.basename(interest_feature_file).split('.')[0]
        df.to_excel(writer, sheet_name=sheet_name)
