#!/usr/bin/env python
# encoding: utf-8

from sys import argv

with open(argv[1], 'r') as inf:
    data = inf.read()

n = int(data.split()[0])  # Months
k = int(data.split()[1])  # Rabbit pairs produced per pair

def rabbits(n, k):
    fib_seq = []

    for i in range(n):
        if i < 2:
            fib_seq.append(1)
        else:
            adults = fib_seq[-1]     # Previous month's rabbits are all adults
            babys = fib_seq[-2] * k  # All rabbits from 2 months ago are having babys
            fib_seq.append(adults + babys)

    return fib_seq[-1]

print rabbits(n, k)
