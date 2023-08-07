""" 
From a given list of genes, an annotation file and a genome fasta file, this
pipeline will return the sequences of all genes of interest and the sequences
of all their most similar gene.

During find_similar_gene, a pickle file will be created ({list_name}_gene_pairs.pkl).
This file will contains a dictionary of all the gene pairs and thus can be used 
for further analysis. 

N.B. This pipeline use bedtools to retrieve the sequences. Make sure to 
install it before running it.
N.B.2. Because this pipeline uses Bedtools, wget and gzip commands, it must be 
run on a Unix terminal. 
"""

import os

configfile: "config.yaml"

lists = [os.path.basename(list_).split('.')[0] for list_ in config['lists']]
annotation_name = config['annotation_url'].rsplit('/', 1)[-1].rsplit('.', 1)[0]
fasta_name = config['fasta_url'].rsplit('/', 1)[-1].rsplit('.', 1)[0]

wildcard_constraints:
    list_name = '|'.join(lists),
    distribution = '|'.join(['', '_similar'])



rule all:
    input:
        expand('data/{genes_list}.fa', genes_list = lists),
        expand('data/{genes_list}_similar.fa', genes_list=lists),
        'data/Similar_Genes_Features.xlsx'

rule write_features_to_excel:
    """ Write an excel file to compare features of similar genes pairs. """
    input:
        gene_pairs = expand('data/{list_name}_gene_pairs.json', list_name=lists),
        interest_features = expand('data/{list_name}_FEATURES.pkl', list_name=lists),
        other_features = 'data/other_genes_FEATURES.pkl'
    output:
        excel = 'data/Similar_Genes_Features.xlsx'
    script:
        'scripts/features_to_excel.py'


rule get_regions_sequences:
    """ Use the genome fasta file and the precedently created bed files to retrieve
        the sequences of all genes of interest and similar genes. """
    input:
        bed = 'data/{list_name}{distribution}.bed',
        fasta = f'data/{fasta_name}'
    output:
        sequences = 'data/{list_name}{distribution}.fa'
    run:
        shell("bedtools getfasta -fi {input.fasta} -bed {input.bed} -name -s | fold -w 60 > {output.sequences}"),


rule find_similar_gene:
    """ Read the genes-of-interest bed file and, for each gene, find the most 
        similar "other" gene. These similar genes will then be used as the null
        distribution. """
    input:
        other_bed = 'data/other_genes.bed',
        other_features = 'data/other_genes_FEATURES.pkl',
        interest_features = 'data/{list_name}_FEATURES.pkl'
    output:
        similar = 'data/{list_name}_similar.bed',
        gene_pairs = 'data/{list_name}_gene_pairs.json'
    script:
        'scripts/find_similar_gene.py'


rule create_other_genes_bedfile:
    """ Create a bed file containing only the lines about genes which are not
        in the given genes list. If multiple lists are given, all list are merged
        together for this rule.  """
    input:
        genes_lists = expand('data/{list_name}.txt', list_name = lists),
        bed = 'data/all_canonical_transcripts.bed'
    output:
        bed = temp('data/other_genes.bed'),
        features = temp('data/other_genes_FEATURES.pkl')
    script:
        'scripts/other_genes_bed.py'


rule create_genes_of_interest_bedfiles:
    """ For each genes list, create a bed file containing only the lines about 
        genes in the list and and compile features of each transcript. """
    input:
        genes = 'data/{list_name}.txt',
        bed = 'data/all_canonical_transcripts.bed'
    output:
        bed = temp('data/{list_name}.bed'),
        features = temp('data/{list_name}_FEATURES.pkl')  # Transcript features

    params:
        # List of genes for which no transcript was found in the annotation file
        unfound = 'data/{list_name}_NO_TRANSCRIPT.txt'
    script:
        'scripts/genes_of_interest_bed.py'


rule convert_annotation_to_bed:
    """ Convert the annotation file to a bed file with only the longest valid 
        transcript for each gene. """
    input:
        gff = f'data/{annotation_name}'
    output:
        bed = 'data/all_canonical_transcripts.bed'
    script:
        'scripts/convert_annotation_to_bed.py'

# rule download_annotation:
#     """ Download and decompress the gencode comprehensive annotation. 
#         COMMENT THIS RULE IF ANNOTATION HAS ALREADY BEEN DOWNLOADED OR IF USING
#         CUSTOM FILE. """
#     output:
#         gff = f"data/{annotation_name}"
#     params:
#         annotation_url = annotation_url
#     shell:
#         """
#         wget -O {output.gff}.gz {params.annotation_url} &&\
#         gzip -d {output.gff}.gz
#         """

# rule download_genome:
#     """ Download the genome fasta file from the given url. 
#         COMMENT THIS RULE IF GENOME HAS ALREADY BEEN DOWNLOADED OR IF USING 
#         CUSTOM FILE. """
#     output:
#         fasta = f"data/{fasta_name}"
#     params:
#         fasta_url = fasta_url
#     shell:
#         """
#         wget -O {output.fasta}.gz {params.fasta_url} &&\
#         gzip -d {output.fasta}.gz
#         """
