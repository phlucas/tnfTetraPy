#!/usr/bin/python -tt

# Senthil/UNLV October 15, 2012
# Given a multifasta file calculate tetranucleotide frequency
# for each of the sequence and output a 257 column, that contains
# the name of the sequence and the 256 possible tetranucleotide frequency
# Requires python version > 2.7 (for collections)

# July 14 2015: Modified since python3 version gives error

# February 24 2017: Fixed bugs
#                   Changed print command to function for python3

# March 3 2017: Forked to tetraMap from tetra.py to run in hadoop
#
# April 10, 2017: Added SEQ to file in order to sort by both seq and tnf
#
# The MAP portion of the MapReduce framework is to form key/value pairs.
# In our case, the key will be a concatonation of the sequence name and
# tetranucleotide, and the value will be 1.  The Reduce program will match
# the same keys and add the values, effectively counting them.
#
from __future__ import print_function

import sys
from random import choice
from Bio import SeqIO
from itertools import product, islice
from collections import Counter
import numpy
import string


def main():
    origlist = make_tetra()
    for i in SeqIO.parse(sys.stdin, 'fasta'):
        fandrseq = i.seq.upper()
        if len(i.seq) >= 2000:
            rcseq = i.seq.reverse_complement()
            fandrseq = (i.seq + rcseq)
            fandrseq = fandrseq.upper()
            movwindow = ["".join(x) for x in window(fandrseq, 4)]
            for quad in movwindow:
                if quad in origlist:
                    #if quad == "AAAA":
                    print(i.id + '\t' + quad)


def make_tetra():
    a = []
    for x in product('ATGC', repeat=4):
        a.append(''.join(x))
    return a


def window(seq, n=4):
    """ Returns a sliding window (of width n) over data from the
        iterable s -> (s0, s1, ... s[n-1]), (s1, s2, ..., sn), ... """
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result

if __name__ == '__main__':
    main()
