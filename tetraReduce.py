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
# April 10, 2017: Added SEQ to file in order to sort by both seq and tnf
#
# The Reduce portion of the MapReduce framework is to read pre-sorted
# key1,key2,value pairs and calculate the frequencies


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
        sys.exit('Simple TetraReduce. Usage: python '
                 + sys.argv[0])
    
    # read values and accumulate them for "print_result"
    current_seq = None
    current_tetra = None
    count = 0
    
    for line in sys.stdin:
        #for line in fileinput.input():
        line = line.strip()
        seq, tetra = line.split('\t')
        
        if seq == current_seq and tetra == current_tetra:
            count += 1
        else:
            if current_seq == seq:
                print(current_seq + '\t' + current_tetra + '\t' + str(count))
                current_tetra = tetra
                count = 1
            else:
                count = 1
                current_seq = seq
                current_tetra = tetra

if current_seq == seq:
    print(current_seq + '\t' + current_tetra + '\t' + str(count))

if __name__ == '__main__':
    main()
