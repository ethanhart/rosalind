#!/usr/bin/env python
# encoding: utf-8

# Each codon consists of three nucleotides.
# Go though an RNA sequence three nucleotides at a time
# and replace the RNA nucleotides with the corresponding
# amino acid. This will produce a protein string.

from sys import argv

__author__ = "Ethan Hart"

# RNA codon to protein lookup table
rna_codon_table = {'UUU': 'F',
                    'CUU': 'L',
                    'AUU': 'I',
                    'GUU': 'V',
                    'UUC': 'F',
                    'CUC': 'L',
                    'AUC': 'I',
                    'GUC': 'V',
                    'UUA': 'L',
                    'CUA': 'L',
                    'AUA': 'I',
                    'GUA': 'V',
                    'UUG': 'L',
                    'CUG': 'L',
                    'AUG': 'M',
                    'GUG': 'V',
                    'UCU': 'S',
                    'CCU': 'P',
                    'ACU': 'T',
                    'GCU': 'A',
                    'UCC': 'S',
                    'CCC': 'P',
                    'ACC': 'T',
                    'GCC': 'A',
                    'UCA': 'S',
                    'CCA': 'P',
                    'ACA': 'T',
                    'GCA': 'A',
                    'UCG': 'S',
                    'CCG': 'P',
                    'ACG': 'T',
                    'GCG': 'A',
                    'UAU': 'Y',
                    'CAU': 'H',
                    'AAU': 'N',
                    'GAU': 'D',
                    'UAC': 'Y',
                    'CAC': 'H',
                    'AAC': 'N',
                    'GAC': 'D',
                    'UAA': 'Stop',
                    'CAA': 'Q',
                    'AAA': 'K',
                    'GAA': 'E',
                    'UAG': 'Stop',
                    'CAG': 'Q',
                    'AAG': 'K',
                    'GAG': 'E',
                    'UGU': 'C',
                    'CGU': 'R',
                    'AGU': 'S',
                    'GGU': 'G',
                    'UGC': 'C',
                    'CGC': 'R',
                    'AGC': 'S',
                    'GGC': 'G',
                    'UGA': 'Stop',
                    'CGA': 'R',
                    'AGA': 'R',
                    'GGA': 'G',
                    'UGG': 'W',
                    'CGG': 'R',
                    'AGG': 'R',
                    'GGG': 'G'}


def rna_to_protein(rna_string):
    """Translate RNA sequence to protein string"""

    proteins = []
    for seq in range(0, len(rna_string) - 1, 3):
        rna_codon = rna_string[seq:seq + 3]
        protein = rna_codon_table[rna_codon]
        if protein == 'Stop':
            return ''.join(proteins)
        else:
            proteins.append(protein)

    return ''.join(proteins)


def main():
    rna_doc = argv[1]
    with open(rna_doc, 'r') as inf:
        rna_string = inf.read().strip()

    protein_string = rna_to_protein(rna_string)
    print protein_string


if __name__ == "__main__":
    main()
