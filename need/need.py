#!/usr/bin/env python
# encoding: utf-8

"""
Given: Two GenBank IDs.

Return: The maximum global alignment score between the DNA strings associated with these IDs.
"""

from sys import argv
from Bio import SeqIO
from Bio import Entrez
import subprocess
import os


__author__ = "Ethan Hart"


def ids_to_fasta(ids):
    """Convert GenBank IDs to FASTA,
    write out FASTA to file with GenBank ID as filename
    """

    Entrez.email = "ethan.john.hart@gmail.com"
    handle = Entrez.efetch(db="nucleotide", id=ids, rettype="fasta")
    records = list(SeqIO.parse(handle, "fasta"))

    for i, r in enumerate(records):
        SeqIO.write(records, ids[i], 'fasta')


def calc_align_score(ids, outf):
    """Run EMBOSS's needle to align FASTA strings"""

    cmd = ["./needle", "-asequence", ids[0], "-bsequence", ids[1], outf,
            "-gapopen", "10", "-gapextend", "1", "-endweight", "-endopen",
            "10", "-endextend", "1"]

    subprocess.call(cmd, stderr=open(os.devnull, 'w'))


def read_score(outf):
    """Get score from needle output file"""

    with open(outf, 'r') as inf:
        data = inf.readlines()

    scores = []
    for line in data:
        if "Score" in line:
            scores.append(line.split()[-1])

    return scores[-1]


def cleanup(ids, outf):
    """Remove intermediate files"""

    os.remove(outf)
    for fname in ids:
        os.remove(fname)


def main():
    with open(argv[1], 'r') as inf:
        data = inf.read()

    ids = data.strip().split()
    outf = "tempout"

    ids_to_fasta(ids)
    calc_align_score(ids, outf)
    score = read_score(outf)
    cleanup(ids, outf)

    print score


if __name__ == "__main__":
    main()
