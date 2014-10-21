#!/usr/bin/env python
# encoding: utf-8

"""
Given: The UniProt ID of a protein.

Return: A list of biological processes in which the protein is involved
(biological processes are found in a subsection of the protein's
"Gene Ontology" (GO) section).
"""

# Requires biopython package
from Bio import ExPASy
from Bio import SwissProt
from sys import argv


__author__ = "Ethan Hart"


def get_processes(record):
    cross_refs = record.cross_references
    procs = []
    for cr in cross_refs:
        if cr[0] == 'GO':
            for i in cr:
                if i.startswith('P:'):
                    procs.append(i.replace('P:', ''))

    return procs


def main():
    with open(argv[1], 'r') as inf:
        data = inf.read()

    protein = data.strip()
    handle = ExPASy.get_sprot_raw(protein)
    record = SwissProt.read(handle)

    for i in get_processes(record):
        print i


if __name__ == "__main__":
    main()
