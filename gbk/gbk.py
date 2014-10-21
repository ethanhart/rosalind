#!/usr/bin/env python
# encoding: utf-8

"""
Given: A genus name, followed by two dates in YYYY/M/D format.

Return: The number of Nucleotide GenBank entries for the given genus that were
published between the dates specified.
"""

from sys import argv
from Bio import Entrez


def search_entrez(gene, start_date, end_date):
    Entrez.email = "ethan.john.hart@gmail.com"
    organism = '"{0}"[Organism]'.format(gene)
    date_field = '"{0}"[Publication Date]'
    start_date_field = date_field.format(start_date)
    end_date_field = date_field.format(end_date)
    date_span = '{0} : {1}'.format(start_date_field, end_date_field)
    term = '{0} AND {1}'.format(organism, date_span)
    handle = Entrez.esearch(db="nucleotide", term=term)
    record = Entrez.read(handle)
    return record["Count"]


def main():
    with open(argv[1], 'r') as inf:
        data = inf.readlines()

    gene = data[0].strip()
    start_date = data[1].strip()
    end_date = data[2].strip()

    print search_entrez(gene, start_date, end_date)


if __name__ == "__main__":
    main()
