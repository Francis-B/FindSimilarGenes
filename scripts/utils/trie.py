#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
This class is used to create a trie structure from a list of strings and then
return a regex pattern. This allow a much faster search than a simple regex
pattern union. 

This script is a adapted version of Eric Duminil's script, which can be found at:
https://gist.github.com/EricDuminil/8faabc2f3de82b24e5a371b6dc0fd1e0
"""

import re


class Trie:
    """ Regexp::Trie in python. Creates a Trie out of a list of words. The trie 
        can be exported to a Regexp pattern. The corresponding Regexp should 
        match much faster than a simple Regexp union."""

    def __init__(self):
        self.data = {}

    def add_list(self, word_list):
        """ Add a list of strings to the trie."""
        for word in word_list:
            self.add(word)

    def add(self, word):
        """ Add one string to the trie"""
        ref = self.data
        for char in word:
            ref[char] = ref.get(char, {})
            ref = ref[char]
        ref[''] = {}

    def flush(self):
        """ Remove all words from trie. """
        self.data = {}

    def dump(self):
        """ Return the trie """
        return self.data

    def _pattern(self, data, is_root=False):
        # Last and only node in dict
        if len(data) == 1 and '' in data:
            return None

        deeper = []
        current = []
        end_reached = False

        for char in sorted(data):
            if char == '':
                end_reached = True
                continue

            # Call function on deeper nodes
            recurse = self._pattern(data[char])
            if recurse is None:     # Reached end of branch
                current.append(re.escape(char))
            else:
                deeper.append(re.escape(char) + recurse)

        final = list(deeper)

        if current:  # If this node is also a branch
            if len(current) == 1:
                final.append(current[0])
            else:
                final.append('[' + ''.join(current) + ']')

        if len(final) == 1:
            result = final[0]
        elif not is_root:
            result = '(?:' + '|'.join(sorted(final)) + ')'
        else:
            result = '(?:' + '|'.join(sorted(final)) + ')'

        if end_reached:
            if not deeper:  # If deeper is empty
                result += '?'
            else:
                result = '(?:' + result + ')?'

        return result

    def pattern(self):
        """ Return regex pattern. """
        return self._pattern(self.dump(), is_root=True)
