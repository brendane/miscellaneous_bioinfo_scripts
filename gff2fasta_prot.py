#!/usr/bin/env python2.7
"""
    Given a gff file and a fasta file, output a file with protein
    translations of genes.

    gff2fasta_prot.py [--prefix <prefix>] <gff> <fasta>
"""

import argparse
import csv
import sys

from Bio import SeqIO, Seq

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--prefix')
parser.add_argument('gff')
parser.add_argument('fasta')
args = parser.parse_args()

if args.prefix is not None:
    prefix = args.prefix
else:
    prefix = ''

idx = SeqIO.index(args.fasta, 'fasta')
lengths = {i:len(idx[i]) for i in idx}

with open(args.gff, 'rb') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        chrom = row[0]
        start = int(row[3]) - 1
        end = int(row[4])
        strand = row[6]
        lt = {x.split('=')[0]:x.split('=')[1]
              for x in row[8].split(';')}['locus_tag']
        if end > lengths[chrom]-1:
            #seq = idx[chrom][start:] + \
            #        idx[chrom][0:(end-lengths[chrom])]
            seq = idx[chrom][start:]
        else:
                seq = idx[chrom][start:end]
        overhang = 3 - (len(seq) % 3)
        if overhang != 3:
            seq.seq += 'N' * overhang
        if strand == '-':
            seq = seq.reverse_complement()
        sys.stdout.write('>' + prefix + lt + '\n')
        sys.stdout.write(str(seq.seq.translate()) + '\n')
