#!/usr/bin/env python
# encoding: utf-8


"""
Given: A collection of k (k≤100) DNA strings of length at most 1 kbp each in
FASTA format.

Return: A longest common substring of the collection. (If multiple solutions
exist, you may return any single solution.)
"""


from sys import argv
import re


__author__ = "Ethan Hart"


def read_data(data_file):
    with open(data_file, 'r') as inf:
        data = inf.read()
    return data


def extract_sequence_from_fasta(data):
    """
    Extract DNA sequence from the following format
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
        sequences.append(dna)

    return sequences


def get_substring(string, length):
    """Return all substrings of given length"""

    num_substrings = (len(string) - length) + 1
    substrings = []
    for s in range(num_substrings):
        substrings.append(string[s:length + s])

    return substrings


def find_longest_common_substring(sequences):
    sorted_seq = sorted(sequences, key=len)  # Shortest string is first

    # We know that item[0] will be the max length of substring
    short = sorted_seq[0]
    if short in sorted_seq[1]:
        return sorted_seq[0]

    # Return longest substring; try as long as substring is not 0
    dec = 1
    while dec <= len(short):
        substrings = get_substring(short, len(short) - dec)

        for s in substrings:
            if all(s in i for i in sorted_seq):
                return s
        dec += 1


def main():
    data_file = argv[1]
    data = read_data(data_file)
    sequences = extract_sequence_from_fasta(data)

    print find_longest_common_substring(sequences)


if __name__ == "__main__":
    main()
