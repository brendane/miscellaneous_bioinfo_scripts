#!/usr/bin/env python2.7
"""
    Calculate FPKM values using the output of htseq-count and a GFF3
    file.

    The total read count includes just the reads (or fragments) counted
    by htseq.

    fpkm_from_htseq_count.py --output <output prefix>
        <feature id> <feature type> <input file> <gff file>
"""

import argparse
import csv

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('feature_id')
parser.add_argument('feature_type')
parser.add_argument('htseq')
parser.add_argument('gff')
args = parser.parse_args()

outpre = args.output
fi = args.feature_id
ft = args.feature_type
infile = args.htseq
gff_file = args.gff

gene_lengths = {}
with open(outpre + '.lengths.tsv', 'w') as oh:
    oh.write('gene\tlength\n')
    with open(gff_file, 'rb') as ih:
        rdr = csv.reader(ih, delimiter='\t')
        for row in rdr:
            if row[2] != ft:
                continue
            info = {}
            for x in row[8].strip().split(';'):
                k = x.split('=')[0]
                v = '='.join(x.split('=')[1:])
                info[k] = v
            l = int(row[4]) - int(row[3]) + 1
            gene_lengths[info[fi]] = l
            oh.write(info[fi] + '\t' + str(l) + '\n')

count = {}
with open(infile, 'rb') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        if row[0].startswith('__'):
            continue
        count[row[0]] = int(row[1])
total = sum(count.values())

with open(outpre + '.fpkm.tsv', 'w') as oh:
    for g, c in count.iteritems():
        gl = gene_lengths[g]
        fpkm = (10E9 * c) / (total * gl)
        oh.write(g + '\t' + str(fpkm) + '\n')
