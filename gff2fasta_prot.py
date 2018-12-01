#!/usr/bin/env python2.7
"""
    Given a gff file and a fasta file, output a file with protein
    translations of genes.

    gff2fasta_prot.py [--id <id field (locus_tag)>] [--prefix <prefix>] <gff> <fasta>
"""

import argparse
import csv
import sys

from Bio import SeqIO, Seq

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--prefix')
parser.add_argument('--id', default='locus_tag')
parser.add_argument('gff')
parser.add_argument('fasta')
args = parser.parse_args()

id_field = args.id
if args.prefix is not None:
    prefix = args.prefix
else:
    prefix = ''

seqs = {rec.id:rec for rec in SeqIO.parse(args.fasta, 'fasta')}
lengths = {i:len(seqs[i]) for i in seqs}

with open(args.gff, 'rb') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        chrom = row[0]
        start = int(row[3]) - 1
        end = int(row[4])
        strand = row[6]
        lt = {x.split('=')[0]:x.split('=')[1]
              for x in row[8].split(';')}[id_field]
        if end > lengths[chrom]-1:
            #seq = seqs[chrom][start:] + \
            #        seqs[chrom][0:(end-lengths[chrom])]
            seq = seqs[chrom][start:]
        else:
                seq = seqs[chrom][start:end]
        overhang = 3 - (len(seq) % 3)
        if overhang != 3:
            seq.seq += 'N' * overhang
        if strand == '-':
            seq = seq.reverse_complement()
        sys.stdout.write('>' + prefix + lt + '\n')
        sys.stdout.write(str(seq.seq.translate()) + '\n')
