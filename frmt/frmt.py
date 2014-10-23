#!/usr/bin/env python
# encoding: utf-8

"""
Given: A collection of n (n≤10) GenBank entry IDs.

Return: The shortest of the strings associated with the IDs in FASTA format.
"""

from sys import argv
from Bio import Entrez
from Bio import SeqIO


__author__ = "Ethan Hart"


def get_nuc_fasta(nucleotide_ids):
    """nucleotide_ids should be string of
    nucleotide ids separated by a comma
    """

    Entrez.email = "email@email.com"
    handle = Entrez.efetch(db="nucleotide", id=[nucleotide_ids], rettype="fasta")
    records = list(SeqIO.parse(handle, "fasta"))

    shortest = min(records, key=lambda x: len(x.seq))
    print shortest.format('fasta')


def main():
    with open(argv[1], 'r') as inf:
        data = inf.read()

    nucleotide_ids = ', '.join(data.strip().split())
    get_nuc_fasta(nucleotide_ids)


if __name__ == "__main__":
    main()
