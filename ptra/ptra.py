#!/usr/bin/env python
# encoding: utf-8

"""
Given: A DNA string s of length at most 10 kbp, and a protein string translated
by s.

Return: The index of the genetic code variant that was used for translation.
(If multiple solutions exist, you may return any one.)
"""

from sys import argv
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


def rna_to_protein(rna_string):
    """Translate RNA sequence to protein string"""

    proteins = []
    for seq in range(0, len(rna_string) - 1, 3):
        rna_codon = rna_string[seq:seq + 3]
        protein = rna_codon_table[rna_codon]
        if protein == 'Stop':
            pass
            #return ''.join(proteins)
        else:
            proteins.append(protein)

    return ''.join(proteins)


def dna2rna(dna):
    """
    Convert DNA string to RNA by replaceing
    all Ts with Us
    """

    rna = dna.replace('T', 'U')

    return rna.strip()


def main():
    with open(argv[1], 'r') as inf:
        data = inf.readlines()

    dna = data[0].strip()
    protein = data[1].strip()

    rna_seq = dna2rna(dna)
    rna_as_prot = rna_to_protein(rna_seq)

    index = 1
    while True:
        if rna_as_prot.startswith(protein):
            protein_index = index
            rna_index = (index * 3) - 2
            print rna_index
            exit()
        else:
            rna_as_prot = rna_as_prot[1:]
            index += 1


if __name__ == "__main__":
    main()
