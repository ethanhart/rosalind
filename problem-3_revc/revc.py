#!/usr/bin/env python
# encoding: utf-8

from sys import argv

"""
Given: A DNA string s of length at most 1000 bp.

Return: The reverse complement sc of s.
"""

__author__ = "Ethan Hart"


def reverse_compliment(data):
    """Returns the compliment string"""

    compliments = {"A": "T", "T": "A", "G": "C", "C": "G"}

    new_sequence = []

    for ind, char in enumerate(data):
        new_sequence.insert(0, compliments[char])

    return ''.join(new_sequence)


def main():
    with open(argv[1], 'r') as inf:
        dna = inf.read().strip()

    rc = reverse_compliment(dna)
    print rc


if __name__ == "__main__":
    main()
