#!/usr/bin/env python3
"""
    Produce an aligned fasta file from the output of
    parse_blat_psl_to_pairwise_coords.py.

    aligned_fasta_from_blat.py [--spacing (3)] --output <output file>
        <alignment> <reference genome 1> <reference genome 2>

    Reference genome 1 should be the reference given to blat -- listed
    first in the blat command and in the 3rd and 4th columns of the
    alignment file.
"""

import argparse
import csv

from Bio import SeqIO, Seq
import pysam

def strand_match(b, strand):
    if strand == '+':
        return b
    else:
        return Seq.complement(b)

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--spacing', default=3, type=int)
parser.add_argument('aln')
parser.add_argument('ref1')
parser.add_argument('ref2')
args = parser.parse_args()

## Read reference genomes
ref1 = {}
for rec in SeqIO.parse(args.ref1, 'fasta'):
    ref1[rec.id] = str(rec.seq)
ref2 = {}
for rec in SeqIO.parse(args.ref2, 'fasta'):
    ref2[rec.id] = str(rec.seq)

## Read alignment file
s1 = []
s2 = []
with open(args.aln, 'rt') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    c = None
    p = None
    first = True
    for i, row in enumerate(rdr):
        c1 = row[2]
        p1 = int(row[3]) - 1
        c2 = row[0]
        p2 = int(row[1]) - 1
        st = row[4]
        if (c1 != c or (p1 != p-1 and p1 != p+1)) and not first:
            s1.append('-' * args.spacing)
            s2.append('-' * args.spacing)
        first = False
        p = p1
        c = c1
        if p1 < 0:
            s1.append('-')
            s2.append(strand_match(ref2[c2][p2], st))
        elif p2 < 0:
            s1.append(ref1[c1][p1])
            s2.append('-')
        else:
            s1.append(ref1[c1][p1])
            s2.append(strand_match(ref2[c2][p2], st))

with open(args.output, 'w') as oh:
    oh.write('>Ref\n')
    oh.writelines(s1)
    oh.write('\n')
    oh.write('>Qry\n')
    oh.writelines(s2)
    oh.write('\n')
