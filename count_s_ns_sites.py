#!/usr/bin/env python2
"""
    Count the number of potential synonymous and non-synonymous
    mutations using a very simple method.

    count_s_ns_sites.py --output <output file> <reference> <gff>
"""

import argparse
import csv
import itertools
import re

from Bio import SeqIO

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.izip_longest(*args, fillvalue=fillvalue)

# codon :proportion of mutations that are synonymous assuming
# equal base frequencies (rounded to 3 digits)
codons = {
          'ACC':(0.0, 0.0, 1.0),
          'ATG':(0.0, 0.0, 0.0),
          'AAG':(0.0, 0.0, 0.333),
          'AAA':(0.0, 0.0, 0.333),
          'ATC':(0.0, 0.0, 0.667),
          'AAC':(0.0, 0.0, 0.333),
          'ATA':(0.0, 0.0, 0.667),
          'AGG':(0.333, 0.0, 0.333),
          'CCT':(0.0, 0.0, 1.0),
          'CTC':(0.0, 0.0, 1.0),
          'AGC':(0.0, 0.0, 0.333),
          'ACA':(0.0, 0.0, 1.0),
          'AGA':(0.333, 0.0, 0.333),
          'CAT':(0.0, 0.0, 0.333),
          'AAT':(0.0, 0.0, 0.333),
          'ATT':(0.0, 0.0, 0.667),
          'CTG':(0.333, 0.0, 1.0),
          'CTA':(0.333, 0.0, 1.0),
          'ACT':(0.0, 0.0, 1.0),
          'CAC':(0.0, 0.0, 0.333),
          'ACG':(0.0, 0.0, 1.0),
          'CAA':(0.0, 0.0, 0.333),
          'AGT':(0.0, 0.0, 0.333),
          'CAG':(0.0, 0.0, 0.333),
          'CCG':(0.0, 0.0, 1.0),
          'CCC':(0.0, 0.0, 1.0),
          'TAT':(0.0, 0.0, 0.333),
          'GGT':(0.0, 0.0, 1.0),
          'TGT':(0.0, 0.0, 0.333),
          'CGA':(0.333, 0.0, 1.0),
          'CCA':(0.0, 0.0, 1.0),
          'TCT':(0.0, 0.0, 1.0),
          'GAT':(0.0, 0.0, 0.333),
          'CGG':(0.333, 0.0, 1.0),
          'CTT':(0.0, 0.0, 1.0),
          'TGC':(0.0, 0.0, 0.333),
          'GGG':(0.0, 0.0, 1.0),
          'TAG':(0.0, 0.0, 0.333),
          'GGA':(0.0, 0.0, 1.0),
          'TAA':(0.0, 0.333, 0.333),
          'GGC':(0.0, 0.0, 1.0),
          'TAC':(0.0, 0.0, 0.333),
          'TTC':(0.0, 0.0, 0.333),
          'TCG':(0.0, 0.0, 1.0),
          'TTT':(0.0, 0.0, 0.333),
          'TTG':(0.333, 0.0, 0.333),
          'TCC':(0.0, 0.0, 1.0),
          'GAA':(0.0, 0.0, 0.333),
          'TCA':(0.0, 0.0, 1.0),
          'GCA':(0.0, 0.0, 1.0),
          'GTA':(0.0, 0.0, 1.0),
          'GCC':(0.0, 0.0, 1.0),
          'GTC':(0.0, 0.0, 1.0),
          'TGA':(0.0, 0.333, 0.0),
          'GCG':(0.0, 0.0, 1.0),
          'GTG':(0.0, 0.0, 1.0),
          'GAG':(0.0, 0.0, 0.333),
          'GTT':(0.0, 0.0, 1.0),
          'GCT':(0.0, 0.0, 1.0),
          'TTA':(0.333, 0.0, 0.333),
          'GAC':(0.0, 0.0, 0.333),
          'CGT':(0.0, 0.0, 1.0),
          'TGG':(0.0, 0.0, 0.0),
          'CGC':(0.0, 0.0, 1.0)
         } 

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('reference')
parser.add_argument('gff')
args = parser.parse_args()

refidx = SeqIO.index(args.reference, 'fasta')

genes = []
with open(args.output, 'wb') as out:
    out.write('replicon\tpos\tsyn\tnonsyn\n')
    with open(args.gff, 'rb') as handle:
        rdr = csv.reader(handle, delimiter='\t')
        for j, row in enumerate(rdr):
            replicon = row[0]
            start = int(row[3]) - 1
            end = int(row[4])
            strand = row[6]
            name = re.sub('ID=(.+?);.+', '\\1', row[8])
            if strand == '+':
                seq = str(refidx[replicon][start:end].seq).upper()
            else:
                seq = str(refidx[replicon][start:end].seq.reverse_complement()).upper()
            for i, codon in enumerate(grouper(seq, 3)):
                cd = ''.join(codon)
                try:
                    degen = codons[cd]
                except KeyError:
                    degen = (0.0, 0.0, 0.0)
                if strand == '+':
                    p1 = (start + 1) + (i * 3)
                    p2 = p1 + 1; p3 = p1 + 2
                else:
                    p1 = end - (i * 3)
                    p2 = p1 - 1; p3 = p1 - 2
                out.write(replicon + '\t' +
                          str(p1) + '\t' +
                          str(degen[0]) + '\t' + str(1-degen[0]) + '\n')
                out.write(replicon + '\t' +
                          str(p2) + '\t' +
                          str(degen[1]) + '\t' + str(1-degen[1]) + '\n')
                out.write(replicon + '\t' +
                          str(p3) + '\t' +
                          str(degen[2]) + '\t' + str(1-degen[2]) + '\n')
