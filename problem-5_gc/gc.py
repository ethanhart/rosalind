#!/usr/bin/env python
# encoding: utf-8

from sys import argv
import re

"""
Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).

Return: The ID of the string having the highest GC-content, followed by the
GC-content of that string.
"""


__author__ = "Ethan Hart"


def read_data(data_file):
    with open(data_file, 'r') as inf:
        data = inf.read()
    return data


def extract_sequence(data):
    """
    Extract DNA sequence from the following format
    >ID1
    AGCTAGCT
    >ID2
    CCATCTCAGCTAGT
    """

    data = re.sub('\s+', '', data)
    matches = re.findall('>Rosalind_[0-9]{4}[AGCT]+', data)
    sequences = {}
    for i in matches:
        seq_id = i[1:14]
        dna = i[14:]
        sequences[seq_id] = dna

    return sequences


def gc_content(data):
    """
    Count instances of DNA nucleotides in given string,
    return GC-content as percentage
    """

    strings = ["A", "C", "G", "T"]
    counts = {}
    for nuc in strings:
        counts[nuc] = data.count(nuc)

    GC = counts["G"] + counts["C"]
    return (GC / float(sum(counts.values()))) * 100


def main():
    data_file = argv[1]
    data = read_data(data_file)
    sequences = extract_sequence(data)

    gc_percentages = {}
    for k, v in sequences.items():
        gc_percentages[k] = gc_content(v)

    max_key = max(gc_percentages, key=gc_percentages.get)
    print max_key
    print gc_percentages[max_key]


if __name__ == "__main__":
    main()
