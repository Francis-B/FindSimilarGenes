#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
This class make it easier to extract informations about transcript from annotation.
and infer introns' 
genomic localisation from gff3 files. It is used in create_bed.py 
"""


class Transcript():
    """ From a list of gff3 lines about a specific transcript, store genomic
      localisation of interest """

    def __init__(self, gff_line):
        self.regions = []   # List to store all regions of the transcript with add_regions()
        self.transcript_name = gff_line.attributes['transcript_name']
        self.gene_name = gff_line.attributes['gene_name']
        self.chr = gff_line.chr
        self.length = self._get_transcript_length(gff_line)
        self.strand = gff_line.strand
        self.cds_limits = dict.fromkeys(['start', 'end'], None)

    def _get_transcript_length(self, gff_line):
        """ Get the start and end localisation of the transcript from tuple 
            used to initialise the object. """
        if gff_line.type == 'transcript':
            return int(gff_line.end) - int(gff_line.start)
        raise TypeError(
            f'Unable to infer transcript bounds for {gff_line.attributes["transcript_name"]}')

    def add_region(self, gff_line):
        """ Add a new tuple to the instance. All added tuples are stored in 
            self.regions, used to infer introns' localisation and write it to bed. """
        if gff_line.type not in ['start_codon', 'stop_codon']:
            self.regions.append(tuple([gff_line.type, gff_line.start, gff_line.end]))
        elif gff_line.type == 'start_codon':
            self.cds_limits['start'] = gff_line.start
        elif gff_line.type == 'stop_codon':
            self.cds_limits['end'] = gff_line.end

    def is_valid(self):
        """ Return True if the transcript has a start and an end codon and has regions """
        if len([codon for codon in self.cds_limits.values() if codon is None]) == 0 and \
                len(self.regions) != 0:
            return True
        return False

    def _infer_introns_positive_strand(self):
        """ Infer introns localisation for positive strand transcript. Append
            results to self.localisation_data """
        # Sort by start
        sorted_regions = sorted(self.regions, key=lambda tup: tup[1], reverse=False)  

        for i in range(len(sorted_regions) - 1):
            last_region_end = int(sorted_regions[i][2])
            next_region_start = int(sorted_regions[i+1][1])
            # If there is a gap between two regions, it is an intron
            if next_region_start - last_region_end > 1:
                self.regions.append(('intron', last_region_end+1, next_region_start-1))

    def _infer_introns_negative_strand(self):
        # Sort by start
        sorted_regions = sorted(self.regions, key=lambda tup: tup[1], reverse=True)
        for i in range(len(sorted_regions) - 1):
            last_region_start = int(sorted_regions[i][1])
            next_region_end = int(sorted_regions[i+1][2])
            # If there is a gap between two regions, it is an intron
            if last_region_start - next_region_end > 1:
                self.regions.append(
                    ('intron', next_region_end+1, last_region_start-1))

    def _format_bed_line(self, line):
        """ Format a line of self.regions to a bed format. This function is used
            by self.write_to_bed()."""
        try:
            type_, start, end = line
            name = '::'.join([self.gene_name, type_])
            bed_line = '\t'.join(
                [self.chr, str(start), str(end), name, '0', self.strand])
            return f'{bed_line}\n'
        except ValueError:
            raise ValueError('Unable to format bed line')

    def format_to_bed(self):
        """ Infer introns localisation and return all lines of self.regions in 
            a list of bed formatted string lines. """
        if self.strand == '+':
            self._infer_introns_positive_strand()
            sorted_regions = sorted(self.regions, key=lambda tup: int(
                tup[1]), reverse=False)  # sort by start
        if self.strand == '-':
            self._infer_introns_negative_strand()
            sorted_regions = sorted(self.regions, key=lambda tup: int(
                tup[1]), reverse=True)  # sort by start

        return [self._format_bed_line(line) for line in sorted_regions]

    def write_to_bed(self, file):
        """ Write extracted information about transcript to a bed format. The 
            given file must be open in write mode. """
        if self.strand == '+':
            self._infer_introns_positive_strand()
            data = sorted(self.regions, key=lambda tup: int(
                tup[1]), reverse=False)  # sort by start
        if self.strand == '-':
            self._infer_introns_negative_strand()
            data = sorted(self.regions, key=lambda tup: int(
                tup[1]), reverse=True)  # sort by start

        for line in data:
            file.write(self._format_bed_line(line))
