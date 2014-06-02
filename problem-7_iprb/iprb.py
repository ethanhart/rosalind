#!/usr/bin/env python
# encoding: utf-8

from __future__ import division
from sys import argv

"""
Given: Three positive integers k, m, and n, representing a population
containing k+m+n organisms: k individuals are homozygous dominant for a factor,
m are heterozygous, and n are homozygous recessive.

Return: The probability that two randomly selected mating organisms will
produce an individual possessing a dominant allele (and thus displaying the
dominant phenotype). Assume that any two organisms can mate."""

__author__ = "Ethan Hart"


def calcaulate_prob(k, m, n):
    t = k + m + n  # Total population
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

    return prob


def main():
    with open(argv[1], 'r') as inf:
        data = inf.read().strip()

    k, m, n = map(float, data.split())

    prob = calcaulate_prob(k, m, n)
    print prob


if __name__ == "__main__":
    main()
