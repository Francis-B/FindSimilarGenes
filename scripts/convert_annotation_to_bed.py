""" 
This script filters the annotation file to keep only the longest valid transcript
for each gene.
"""
from utils.gff_line import GFFLine
from utils.transcript import Transcript

# Input and output files
GFF = snakemake.input.gff
BED = snakemake.output.bed

# # Test input and output files
# GFF = 'data/gencode.v43.chr_patch_hapl_scaff.annotation.gff3'
# BED = 'data/all_canonical_transcripts.bed'

# Initiate dict to store the canonical transcript of each gene
canonical_transcripts = {}
# Type for which we want informations
features = ['five_prime_UTR', 'three_prime_UTR', 'CDS', 'start_codon', 'stop_codon']

with open(GFF, 'r', encoding='utf-8') as annotation:
    for line in annotation:
        if line.startswith('#'):
            continue

        gff_line = GFFLine(line)
        if gff_line.type == 'transcript':
            # Add last builded transcript to dict
            try:
                # Skip transcripts with no given start or end codon
                if not actual_transcript.is_valid():
                    continue
                # If there is no transcript for the actual gene in dict, add it
                if actual_transcript.gene_name not in canonical_transcripts:
                    canonical_transcripts[actual_transcript.gene_name] = actual_transcript
                # If the transcript stored in dict is shorter than the actual transcript,
                # replace it
                elif canonical_transcripts[actual_transcript.gene_name].length < actual_transcript.length:
                    canonical_transcripts[actual_transcript.gene_name] = actual_transcript
            # Will catch error for the first loop when no transcript exist
            except NameError:
                pass

            # Initialise a new transcript from actual gff line
            actual_transcript = Transcript(gff_line)

        elif gff_line.type in features:
            actual_transcript.add_region(gff_line)

# Write bed file with the longest transcripts for each gene. This automatically 
# infer introns localisation for each genes due to Transcript internal method
with open(BED, 'w', encoding='utf-8') as f:
    for transcript in canonical_transcripts.values():
        f.writelines(transcript.format_to_bed())
