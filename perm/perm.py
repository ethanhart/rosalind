#!/usr/bin/env python
# encoding: utf-8

"""
Given: A positive integer n≤7.

Return: The total number of permutations of length n, followed by a list of all
such permutations (in any order).
"""

from sys import argv


__author__ = "Ethan Hart"


def perm(seq):
    #print 'SEQ: ', seq
    if len(seq) <= 1:
        yield seq
    else:
        for p in perm(seq[1:]):
            for i in range(len(seq) + 1):
                yield p[:i] + seq[0:1] + p[i:]


def main():
    with open(argv[1], 'r') as inf:
        num = inf.read()
    num = int(num.strip())
    seq = []
    for i in range(num):
        seq.append(i + 1)

    all_seqs = []
    perms = list(perm(seq))
    for x in perms:
        if x not in all_seqs:
            all_seqs.append(x)

    print len(all_seqs)
    for y in all_seqs:
        print ' '.join(map(str, y))

if __name__ == "__main__":
    main()
