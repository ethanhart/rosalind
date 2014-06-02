#!/usr/bin/env python
# encoding: utf-8

"""
Given: Two DNA strings s and t (each of length at most 1 kbp).

Return: All locations of t as a substring of s.
"""

from sys import argv

__author__ = "Ethan Hart"


def find_substring(orig_string, substring):
    instances = []
    for i in range(len(orig_string)):
        test_string = orig_string[i:len(substring) + i]
        if substring == test_string:
            instances.append(str(i + 1))

    return ' '.join(instances)


def main():
    with open(argv[1], 'r') as inf:
        data = inf.readlines()

    orig_string = data[0].strip()
    substring = data[1].strip()
    print find_substring(orig_string, substring)


if __name__ == "__main__":
    main()
