#!/usr/bin/env python2.7
"""
    Given a gff file and a fasta file, output a file with protein
    translations of genes.

    gff2fasta_prot.py [--no-dups] [--id <id field (locus_tag)>]
        [--table <trans table (Bacterial)>]
        [--prefix <prefix>] <gff> <fasta>

    Uses Bacterial translation table, and checks if the sequences is a
    true coding sequence. If it is, alternative start codons are translated
    as "M". Ns are added to the end to make a sequence
    with complete codons after reverse-complementing and before
    translation.

    If the id field is not present, the record is skipped. If the
    no-dups flag is set, sequences with the same name and length
    are skipped after they are first seen.
"""

import argparse
import csv
import sys

from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--prefix')
parser.add_argument('--table', default='Bacterial')
parser.add_argument('--id', default='locus_tag')
parser.add_argument('--no-dups', default=False, action='store_true')
parser.add_argument('gff')
parser.add_argument('fasta')
args = parser.parse_args()

table = args.table
id_field = args.id
nd = args.no_dups
if args.prefix is not None:
    prefix = args.prefix
else:
    prefix = ''

seqs = {rec.id:rec for rec in SeqIO.parse(args.fasta, 'fasta')}
lengths = {i:len(seqs[i]) for i in seqs}
seen = {}

with open(args.gff, 'rb') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        if row[0].startswith('#'):
            continue
        chrom = row[0]
        start = int(row[3]) - 1
        end = int(row[4])
        strand = row[6]
        try:
            lt = {x.split('=')[0]:x.split('=')[1]
                  for x in row[8].split(';')}[id_field]
        except KeyError:
            continue
        if end > lengths[chrom]-1:
            #seq = seqs[chrom][start:] + \
            #        seqs[chrom][0:(end-lengths[chrom])]
            seq = seqs[chrom][start:]
        else:
                seq = seqs[chrom][start:end]
        if lt in seen and seen[lt] == len(seq):
            continue
        seen[lt] = len(seq)
        if strand == '-':
            seq = seq.reverse_complement()
        overhang = 3 - (len(seq) % 3)
        if overhang != 3:
            seq.seq += 'N' * overhang
        sys.stdout.write('>' + prefix + lt + '\n')
        try:
            sys.stdout.write(str(seq.seq.translate(table=table, cds=True)) + '\n')
        except TranslationError:
            sys.stdout.write(str(seq.seq.translate(table=table)) + '\n')
