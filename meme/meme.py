#!/usr/bin/env python
# encoding: utf-8

"""
Given: A set of protein strings in FASTA format that share some motif with
minimum length 20.

Return: Regular expression for the best-scoring motif.

Note: MEME Suite (http://meme.nbcr.net/meme/) is required to run this code
"""

from sys import argv
import subprocess


__author__ = "Ethan Hart"


def run_motif(fasta_file):
    cmd = ["./meme", fasta_file, "-text", "-minw", "20"]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    out, err = process.communicate()

    return out


def extract_regex(motif):
    count = 0
    motif_lines = motif.split('\n')
    for line in motif_lines:
        if "regular expression" in line:
            return motif_lines[count + 2]
        count += 1


def main():
    inf = argv[1]

    motif = run_motif(inf)
    regex = extract_regex(motif)

    print regex


if __name__ == "__main__":
    main()
