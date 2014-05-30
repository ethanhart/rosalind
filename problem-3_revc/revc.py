#!/usr/bin/env python
# encoding: utf-8

from sys import argv

with open(argv[1], 'r') as inf:
    data = inf.read().strip()


def reverse_compliment(data):
    """Returns the compliment string"""

    compliments = {"A": "T", "T": "A", "G": "C", "C": "G"}

    new_sequence = []

    for ind, char in enumerate(data):
        new_sequence.insert(0, compliments[char])

    return ''.join(new_sequence)

print reverse_compliment(data)
