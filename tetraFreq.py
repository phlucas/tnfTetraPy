#!/usr/bin/python -tt
#
# Senthil/UNLV October 15, 2012
# Given a multifasta file calculate tetranucleotide frequency
# for each of the sequence and output a 257 column, that contains
# the name of the sequence and the 256 possible tetranucleotide frequency
# Requires python version > 2.7 (for collections)
#
# July 14 2015: Modified since python3 version gives error
#
# February 24 2017: Fixed bugs
#                   Changed print command to function for python3
#
# March 3 2017: Forked to tetraReduce from tetra.py to run in hadoop
#
#
# April 10, 2017: Added SEQ to file in order to sort by both seq and tnf
#
# tetraFreq written to complete the actual frequency counts

from __future__ import print_function

import sys
from random import choice
from Bio import SeqIO
from itertools import product, islice
from collections import Counter
import numpy
import string

import fileinput

def main():
    if len(sys.argv) != 1:
        sys.exit('Simple TetraFreq Usage: python '
                 + sys.argv[0])
    
    origlist = make_tetra()

    cnt = Counter()
    print ('\t', end='')
    for x in origlist:
        print (x + '\t', end='')

# read values and accumulate them for "print_result"
current_seq = None
    total_count = 0
    
    for line in sys.stdin:
        #for line in fileinput.input():
        line = line.strip()
        seq, tetra, count = line.split('\t')
        
        if seq == current_seq:
            cnt[tetra] = int(count)
        else:
            if current_seq:
                print_result(current_seq, origlist, cnt)
                current_seq = seq
                cnt[tetra] = int(count)
            else:
                current_seq = seq
                cnt[tetra] = int(count)
                current_seq = seq

if current_seq == seq:
    print_result(current_seq, origlist, cnt)
    
    print('\n')
    fileinput.close()

def make_tetra():
    a = []
    for x in product('ATGC', repeat=4):
        a.append(''.join(x))
    return a

def print_result(cname, ol, cf):
    print ('\n' +  cname + '\t', end='')
    # totaltet means the sum of occurrences of all tetranucleotide in
    # the sequence
    totaltet = sum(cf.values())
    # Nobel et al 1998
    for x in ol:
        print ("%16.14f"% (cf[x]/float(totaltet)), end='')
        print ('\t', end='')

if __name__ == '__main__':
    main()
