#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
This script contains GFFLine class and read_annotation function which make it 
easier to work with annotation.
"""


class GFFLine():
    """ This class create an object from a gff line and make it easier to extract
        information."""

    def __init__(self, gff_line):
        self.line = gff_line.split('\t')
        self.chr = self.line[0]
        self.source = self.line[1]
        self.type = self.line[2]
        self.start = int(self.line[3])
        self.end = int(self.line[4])
        self.score = self.line[5]
        self.strand = self.line[6]
        self.phase = self.line[7]
        self.attributes = self._get_attributes(self.line[8])


    def _get_attributes(self, attributes):
        attributes_dict = {}
        _attributes = attributes.split(';')
        for attrib in _attributes:
            name, value = attrib.split('=')
            attributes_dict[name] = value
        return attributes_dict

    def get_line(self):
        " Return the original line"
        return '\t'.join(self.line)


def read_annotation(annotation_iterable):
    """ Make a GFF3Line object from an annotation's iterable.
    
        Return:
            - GFF3Line object if there is a line to read
            - None if there is no line to read 
    """
    try:
        line = next(annotation_iterable)
        while line.startswith('#'):
            line = next(annotation_iterable) # skip comments header
        return GFFLine(line)
    except StopIteration:
        return None
    