#!/usr/bin/env python
# encoding: utf-8

"""
Given: Six positive integers, each of which does not exceed 20,000. The
integers correspond to the number of couples in a population possessing each
genotype pairing for a given factor. In order, the six given integers represent
the number of couples having the following genotypes (with the probability of
the offspring expressing the dominant trait):

    AA-AA = 1
    AA-Aa = 1
    AA-aa = 1
    Aa-Aa = .75
    Aa-aa = .5
    aa-aa = 0

Return: The expected number of offspring displaying the dominant phenotype in
the next generation, under the assumption that every couple has exactly two
offspring.
"""


__author__ = "Ethan Hart"


from sys import argv


def read_data(data_file):
    with open(data_file, 'r') as inf:
        data = inf.read().strip()
    return data


def main():
    data_file = argv[1]
    data = read_data(data_file)
    couples = map(int, data.split())

    # Prob of offspring expressing dominant trait
    dom_prob = [1.0, 1.0, 1.0, 0.75, 0.5, 0.0]

    expected_offspring = sum([p * c for p, c in zip(dom_prob, couples)])
    print expected_offspring * 2  # For both offspring

if __name__ == "__main__":
    main()
