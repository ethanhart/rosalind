#!/usr/bin/env python
# encoding: utf-8

from sys import argv

with open(argv[1], 'r') as inf:
    data = inf.read()

def dna2rna(dna):
    """
    Convert DNA string to RNA by replaceing
    all Ts with Us
    """

    rna = dna.replace('T', 'U')

    return rna.strip()

print dna2rna(data)
