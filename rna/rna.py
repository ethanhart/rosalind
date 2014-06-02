#!/usr/bin/env python
# encoding: utf-8

"""
Given: A DNA string t having length at most 1000 nt.

Return: The transcribed RNA string of t.
"""

from sys import argv

__author__ = "Ethan Hart"


def dna2rna(dna):
    """
    Convert DNA string to RNA by replaceing
    all Ts with Us
    """

    rna = dna.replace('T', 'U')

    return rna.strip()


def main():
    with open(argv[1], 'r') as inf:
        dna = inf.read().strip()

    rna = dna2rna(dna)
    print rna

if __name__ == "__main__":
    main()
