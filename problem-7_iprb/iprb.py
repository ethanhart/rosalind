#!/usr/bin/env python
# encoding: utf-8

from __future__ import division
from sys import argv

"""
Given: a population containing k+m+n organisms:
        k individuals are homozygous dominant for a factor
        m are heterozygous
        n are homozygous recessive
Return: The probability that two randomly selected organisms will produce an
        individual possessing a dominant allele. Any two organisms can mate.
"""

with open(argv[1], 'r') as inf:
    data = inf.read().strip()

k, m, n = map(float, data.split())

t = k + m + n
prob = 1

# Find all parent combinations that will not allow
# for dominant allele and subtract them from 1

# Prob both parents being homozygous recessive
# Prob first parent homozygous recessive = (num of homozygous recessive (n) / total parents (t))
# Prob second parent homozygous recessive = (n-1/t-1) (since we removed 1 parent from pool)
p_homozygous_rec_x2 = (n / t) * ((n - 1) / (t - 1))
prob -= p_homozygous_rec_x2

# Prob parent being homozygous recessive = (n/t)
# Prob parent being heterozygous = (m/t-1) * prob of having recessive allele (0.5)
p_homozygous = (n / t)
p_heterozygous_rec = (m / (t - 1)) * 0.5
p_homozygous_heterozygous_rec = (p_homozygous * p_heterozygous_rec) * 2
prob -= p_homozygous_heterozygous_rec

# Prob first parent being heterozygous (m/t)
# Prob second parent being heterozygous (m-1/t-1)
# Probability of recessive allele 0.5 (for each parent)
p_heterozygous_rec_x2 = (m / t) * ((m - 1) / (t - 1)) * .25

prob -= p_heterozygous_rec_x2

print prob
