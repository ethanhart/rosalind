#!/usr/bin/env python
# encoding: utf-8

"""
Given: A protein string P of length at most 1000 aa.

Return: The total weight of P. Consult the monoisotopic mass table.
"""

from sys import argv


__author__ = "Ethan Hart"


# Monoisotopic Mass Table
mmt = {'A': 71.03711,
       'C': 103.00919,
       'D': 115.02694,
       'E': 129.04259,
       'F': 147.06841,
       'G': 57.02146,
       'H': 137.05891,
       'I': 113.08406,
       'K': 128.09496,
       'L': 113.08406,
       'M': 131.04049,
       'N': 114.04293,
       'P': 97.05276,
       'Q': 128.05858,
       'R': 156.10111,
       'S': 87.03203,
       'T': 101.04768,
       'V': 99.06841,
       'W': 186.07931,
       'Y': 163.06333}


def main():
    with open(argv[1], 'r') as inf:
        protein_string = inf.read().strip()

    total_weight = 0
    for char in protein_string:
        total_weight += mmt[char]

    print total_weight


if __name__ == "__main__":
    main()
