#!/usr/bin/env python
# encoding: utf-8

from sys import argv


with open(argv[1], 'r') as inf:
    lines = inf.readlines()

a = lines[0].strip()
b = lines[1].strip()


def hamming_distance(str1, str2):

    hd = 0
    for x, y in zip(str1, str2):
        if x != y:
            hd += 1

    return hd


print hamming_distance(a, b)
