#!/usr/bin/env python
# encoding: utf-8

"""
Given: A protein string of length at most 1000 aa.

Return: The total number of different RNA strings from which the protein could
have been translated, modulo 1,000,000. (Don't neglect the importance of the
stop codon in protein translation.)
"""

from sys import argv
from collections import Counter
from operator import mul
import re

__author__ = "Ethan Hart"
__date__ = '2015-08-07'

# RNA codon to protein lookup table
rna_codon_table = {'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V', 'UUC': 'F',
                   'CUC': 'L', 'AUC': 'I', 'GUC': 'V', 'UUA': 'L', 'CUA': 'L',
                   'AUA': 'I', 'GUA': 'V', 'UUG': 'L', 'CUG': 'L', 'AUG': 'M',
                   'GUG': 'V', 'UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A',
                   'UCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A', 'UCA': 'S',
                   'CCA': 'P', 'ACA': 'T', 'GCA': 'A', 'UCG': 'S', 'CCG': 'P',
                   'ACG': 'T', 'GCG': 'A', 'UAU': 'Y', 'CAU': 'H', 'AAU': 'N',
                   'GAU': 'D', 'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
                   'CAA': 'Q', 'AAA': 'K', 'GAA': 'E', 'CAG': 'Q', 'AAG': 'K',
                   'GAG': 'E', 'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G',
                   'UGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G', 'CGA': 'R',
                   'AGA': 'R', 'GGA': 'G', 'UGG': 'W', 'CGG': 'R', 'AGG': 'R',
                   'GGG': 'G', 'UGA': 'Stop', 'UAG': 'Stop', 'UAA': 'Stop'}


def main():
    with open(argv[1], 'r') as inf:
        data = inf.read()

    r = []
    for l in data.strip():
        r.append(l)
    r.append('Stop')

    item_counts = Counter(val for val in rna_codon_table.values())

    counts = []
    for i in r:
        counts.append(item_counts[i])
    print reduce(mul, counts, 1) % 1000000

if __name__ == "__main__":
    main()
