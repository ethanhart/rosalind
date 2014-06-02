#!/usr/bin/env python
# encoding: utf-8

"""
Given: A collection of at most 10 DNA strings of equal length (at most 1 kbp)
in FASTA format.

Return: A consensus string and profile matrix for the collection. (If several
possible consensus strings exist, then you may return any one of them.)
"""

from sys import argv
import operator
import re

__author__ = "Ethan Hart"


def extract_sequence_from_fasta(data):
    """
    Extract DNA sequence from the FASTA format
    >ID1
    AGCTAGCT
    >ID2
    CCATCTCAGCTAGT
    """

    data = re.sub('\s+', '', data)
    matches = re.findall('>Rosalind_[0-9]{1,}[ACGT]+', data)
    sequences = {}
    for i in matches:
        seq_id = re.findall('Rosalind_[0-9]{1,}', i)[0]
        dna = re.findall('[ACGT]+', i)[0]
        sequences[seq_id] = dna

    return sequences


def create_empty_matrix(r, c):
    """
    Create an r,c matrix with all values set to 0
    """

    mat = []
    for i in range(r):
        row = []
        for i in range(c):
            row.append(0)
        mat.append(row)

    return mat


def consensus_string(sequences):
    """
    For each sequence (of identical lengths),
    tally nucleotide count at each index
    """

    cols = len(sequences[0])
    table = create_empty_matrix(4, cols)

    for seq in sequences:
        for index, char in enumerate(seq):
            if char == "A":
                table[0][index] += 1
            if char == "C":
                table[1][index] += 1
            if char == "G":
                table[2][index] += 1
            if char == "T":
                table[3][index] += 1

    string_table = [map(str, x) for x in table]

    consensus = []
    for item in range(cols):
        test = {"A": table[0][item],
                "C": table[1][item],
                "G": table[2][item],
                "T": table[3][item]}

        nuc = max(test.iteritems(), key=operator.itemgetter(1))[0]
        consensus.append(nuc)

    print ''.join(consensus)

    rows = ['A', 'C', 'G', 'T']
    for index, r in enumerate(rows):
        print '{0}: {1}'.format(r, ' '.join(string_table[index]))


def main():
    with open(argv[1], 'r') as inf:
        data = inf.read()

    sequences = extract_sequence_from_fasta(data)
    dna_seqs = sequences.values()
    consensus_string(dna_seqs)


if __name__ == "__main__":
    main()
