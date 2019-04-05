#!/usr/bin/env python3
"""
    Proportion of sites in an alignment that are segregating and are
    adjacent to another segregating site. Gaps and Ns are not
    counted as mismatches.

    count_adjacent_snps_alignment_quality.py <fasta input file(s)>
"""

import argparse
import os
import sys

from Bio import AlignIO

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('-d', action='store_true', default=False)
parser.add_argument('inputs', nargs='+')
args = parser.parse_args()

if args.d:
    input_files = [args.inputs[0] + '/' + f for f in os.listdir(args.inputs[0])]
else:
    input_files = args.inputs

for fname in input_files:
    aln = AlignIO.read(fname, 'fasta') 
    adj = 0
    if len(aln) > 1:
        segregating = set()
        for i in range(len(aln[0])):
            alleles = set()
            for j in range(len(aln)):
                a = aln[j, i].upper()
                if a != 'N' and a != '-':
                    alleles.add(a)
            if len(alleles) > 1:
                segregating.add(i)
        for i in segregating:
            if i+1 in segregating or i-1 in segregating:
                adj += 1
    print("%s\t%f" % (fname, adj/len(aln[0])))
