#!/usr/bin/env python
# encoding: utf-8

from sys import argv

"""
Given: Two DNA strings s and t of equal length (not exceeding 1 kbp).

Return: The Hamming distance dH(s,t).
"""

__author__ = "Ethan Hart"


def hamming_distance(str1, str2):

    hd = 0
    for x, y in zip(str1, str2):
        if x != y:
            hd += 1

    return hd


def main():
    with open(argv[1], 'r') as inf:
        lines = inf.readlines()

    a = lines[0].strip()
    b = lines[1].strip()
    hd = hamming_distance(a, b)
    print hd


if __name__ == "__main__":
    main()
