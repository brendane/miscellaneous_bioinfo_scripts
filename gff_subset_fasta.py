#!/usr/bin/env python
"""
    Subset a FASTA alignment using a gff file

    gff_subset_fasta.py --output <output prefix> <gff> <fasta>
        <replicon> <feature types>
"""

#==============================================================================#

import argparse
import csv
import gzip
import re

from Bio import AlignIO, SeqIO

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('gff')
parser.add_argument('fasta')
parser.add_argument('replicon')
parser.add_argument('features')
args = parser.parse_args()

replicon = args.replicon
features = set(args.features.split(','))

if args.fasta.endswith('.gz'):
    fhandle = gzip.open(args.fasta)
else:
    fhandle = open(args.fasta)

with fhandle:
    aln = AlignIO.read(fhandle, 'fasta')
    with open(args.gff, 'rb') as handle:
        rdr = csv.reader(handle, delimiter='\t')
        for row in rdr:
            if row[0].startswith('#'):
                continue
            if row[0] != replicon:
                continue
            if row[2] not in features:
                continue
            name = row[8]
            name = [x for x in name.split(';') if x.startswith('ID=')][0]
            name = re.sub('ID=', '', name)
            start = int(row[3]) - 1
            end = int(row[4])
            if row[6] == '+':
                AlignIO.write(aln[:,start:end], args.output + name + '.fasta',
                              'fasta')
            else:
                with open(args.output + name + '.fasta', 'wb') as out:
                    for rec in aln:
                        out.write('>' + rec[start:end].id + '\n')
                        out.write(str(rec[start:end].reverse_complement().seq) +
                                  '\n')
