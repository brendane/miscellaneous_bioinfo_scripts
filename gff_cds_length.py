#!/usr/bin/env python3
"""
    Calculate FPKM values using the output of htseq-count and a GFF3
    file. Sums the length of the CDS elements in the GFF file for each
    gene.

    The total read count includes just the reads (or fragments) counted
    by htseq.

    fpkm_from_htseq_count.py --output <output prefix> [--header ]
        <feature id> <input file> <gff file>

    --header means the file has a header with column names. The first
    column is assumed to be genes and the other columns samples. Without
    --header, assumption is a file with one sample in the second column
    and gene names in the first.

    .[0-9] is removed from the end of gene IDs.
"""

import argparse
import collections
import csv
import re
import sys

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--header', default=False, action='store_true')
parser.add_argument('feature_id')
parser.add_argument('htseq')
parser.add_argument('gff')
args = parser.parse_args()

outpre = args.output
fi = args.feature_id
infile = args.htseq
gff_file = args.gff

gene_lengths = collections.defaultdict(int)
with open(gff_file, 'r') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        if row[2] != 'CDS':
            continue
        info = {}
        for x in row[8].strip().split(';'):
            k = x.split('=')[0]
            v = '='.join(x.split('=')[1:])
            info[k] = v
        l = int(row[4]) - int(row[3]) + 1
        g = re.sub('\\.[0-9]$', '', info[fi])
        gene_lengths[g] += l

with open(outpre + '.lengths.tsv', 'w') as oh:
    oh.write('gene\tlength\n')
    for g in gene_lengths:
        oh.write(g + '\t' + str(gene_lengths[g]) + '\n')

counts = collections.defaultdict(lambda: collections.defaultdict(int))
fields = None
gene_order = []
with open(infile, 'r') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    if args.header:
        fields = next(rdr)
    for row in rdr:
        if row[0].startswith('__'):
            continue
        gene_order.append(row[0])
        if args.header:
            for i, e in enumerate(row):
                if i == 0:
                    continue
                counts[fields[i]][row[0]] = int(row[i])
        else:
            counts['_'][row[0]] = int(row[1])
totals = {k:sum(counts[k].values()) for k in counts}

if args.header:
    with open(outpre + '.fpkm.tsv', 'w') as oh:
        oh.write('gene\t' + '\t'.join(fields[1:]) + '\n')
        for g in gene_order:
            if g not in gene_lengths:
                sys.stderr.write('%s not in gff, skipping\n' % g)
                continue
            gl = gene_lengths[g]
            oh.write(g)
            for sample in fields[1:]:
                c = counts[sample][g]
                fpkm = (1E9 * c) / (totals[sample] * gl)
                oh.write('\t' + str(fpkm))
            oh.write('\n')
else:
    with open(outpre + '.fpkm.tsv', 'w') as oh:
        t = totals['_']
        for g in gene_order:
            if g not in gene_lengths:
                sys.stderr.write('%s not in gff, skipping\n' % g)
                continue
            gl = gene_lengths[g]
            c = counts['_'][g]
            fpkm = (1E9 * c) / (t * gl)
            oh.write(g + '\t' + str(fpkm) + '\n')
