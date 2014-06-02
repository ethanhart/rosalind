#!/usr/bin/env python
# encoding: utf-8

from sys import argv

"""
Given: A DNA string s of length at most 1000 nt.

Return: Four integers (separated by spaces) counting the respective number of
        times that the symbols 'A', 'C', 'G', and 'T' occur in s.
"""

__author__ = "Ethan Hart"


def read_data(data_file):
    with open(data_file, 'r') as inf:
        lines = inf.read()
    return lines


def count_dna_strings(data):
    """
    Count instances of DNA nucleotides in given string
    Set: A C G T
    """

    strings = ["A", "C", "G", "T"]
    counts = map(lambda s: str(data.count(s)), strings)
    return ' '.join(counts)


def main():
    data = argv[1]
    datalines = read_data(data)
    print count_dna_strings(datalines)


if __name__ == "__main__":
    main()
