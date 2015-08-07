#!/usr/bin/env python
# encoding: utf-8

"""
Given: A DNA string s (of length at most 1 kbp) and a collection of substrings
of s acting as introns. All strings are given in FASTA format.

Return: A protein string resulting from transcribing and translating the exons
of s. (Note: Only one solution will exist for the dataset provided.)
"""

from sys import argv
import re

__author__ = "Ethan Hart"
__date__ = '2015-08-07'

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
    sequences = []
    for i in matches:
        seq_id = re.findall('Rosalind_[0-9]{1,}', i)[0]
        dna = re.findall('[ACGT]+', i)[0]
        #sequences[seq_id] = dna
        sequences.append(dna)  # Change to list as order matters, id not needed

    return sequences


def remove_intron(dna_seq, introns):
    """
    Remove intron substring from dna_seq
    """

    for intron in introns:
        dna_seq = re.sub(intron, '', dna_seq)

    return dna_seq


def dna2rna(dna):
    """
    Convert DNA string to RNA by replaceing
    all Ts with Us
    """

    rna = dna.replace('T', 'U')

    return rna.strip()


def main():
    with open(argv[1], 'r') as inf:
        data = inf.read()

    sequences = extract_sequence_from_fasta(data)
    dna_seq = sequences[0]
    introns = sequences[1:]

    dna_less_introns = remove_intron(dna_seq, introns)
    rna_seq = dna2rna(dna_less_introns)
    prot = rna_to_protein(rna_seq)

    print prot


if __name__ == "__main__":
    main()
