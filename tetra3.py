#!/usr/bin/python -tt

# Senthil/UNLV October 15, 2012
# Given a multifasta file calculate tetranucleotide frequency
# for each of the sequence and output a 257 column, that contains
# the name of the sequence and the 256 possible tetranucleotide frequency
# Requires python version > 2.7 (for collections)

# July 14 2015: Modified since python3 version gives error

# February 24 2017: Fixed bugs
#                   Changed print command to function for python3

from __future__ import print_function

import sys
from random import choice
from Bio import SeqIO
from itertools import product, islice
from collections import Counter
import numpy
import string

def main():
    if len(sys.argv) != 2:
        sys.exit('Simple Tetra. Usage: python '
                 + sys.argv[0] + ' <in_fasta>')
#    iupac_N = ('A', 'T', 'G', 'C')
    infile = sys.argv[1]
    f = open(infile, 'rU')
    origlist = make_tetra()
    # print tetramers once
    # old print 'Contigs'
    print ('\t', end='')
    for x in origlist:
        print (x + '\t', end='')
    for i in SeqIO.parse(f,'fasta'):
        # do something with
        # i.id and i.seq
        # rcseq reverse complement sequence
        # fandrseq forward+reverse complement seq
        fandrseq = i.seq.upper()
        if len(i.seq) >= 2000:
            rcseq = i.seq.reverse_complement()
            fandrseq = (i.seq + rcseq)
            fandrseq = fandrseq.upper()
        # print "\n", fandrseq
        # bases = list(fandrseq)
        # for base in range(len(bases)):
        #     if bases[base] == 'N':
        #         bases[base] = choice(iupac_N)
        # fandrseq = ''.join(bases)
        # print fandrseq
            movwindow = ["".join(x) for x in window(fandrseq, 4)]
            cfreq = calc_freq(movwindow, origlist)
            print_result(i.id, origlist, cfreq)
    print ('\n')
    f.close()

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
# This does not make sense so commented out April 7 2015
# (it handles the first nucleotide so removed comment Feb 24, 2017
    if len(result) == n:
        yield result    
    for elem in it:
        result = result[1:] + (elem,)
        yield result
        

def calc_freq(mw, ol):
    cnt = Counter()
    for quad in mw:
        if quad in ol:
            cnt[quad] += 1
    return cnt

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
