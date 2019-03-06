#!/usr/bin/env python3
"""
    Given a gff file and a txt annotation file from ENSEMBL, produce
    input files for inrich.

    gff_and_ensembl_to_inrich.py --output <prefix> <gff> <txt>
"""

import argparse
import collections
import csv
import re
import sys

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('gff')
parser.add_argument('annot')
args = parser.parse_args()


gff_info = {}
gene_order = []
gid_map = {}
with open(args.gff, 'r') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        if row[0].startswith('#'):
            continue
        if row[2] == 'gene':
            chrom = row[0]
            start = int(row[3])
            end = int(row[4])
            annot = row[8]
            gid = re.sub('gene:', '',
                         {x.split('=')[0]:x.split('=')[1] for x in annot.split(';')}['ID'])
            gff_info[gid] = (chrom, start, end)
            gene_order.append(gid)
            gid_map[gid] = len(gid_map) + 1
        else:
            continue


go_annot = {}
gene_go = collections.defaultdict(set)
gene_names = {}
with open(args.annot, 'r') as ih:
    rdr = csv.DictReader(ih, delimiter='\t')
    for row in rdr:
        gid = row['Gene stable ID']
        if gid not in gff_info:
            sys.stderr.write('%s is only in ensembl file\n' % gid)
        desc = row['Gene description']
        go = row['GO term accession']
        go_desc = row['GO term name']
        gene_go[gid].add(go)
        go_annot[go] = go_desc
        gene_names[gid] = desc

with open(args.output + '.genelocs.tsv', 'w') as oh:
    for gid in gene_order:
        chrom, start, stop = gff_info[gid]
        try:
            n = gene_names[gid]
        except KeyError:
            n = ''
        oh.write(chrom + '\t' + str(start) + '\t' + str(stop) + '\t' +
                 str(gid_map[gid]) + '\t' + n + '\n')

with open(args.output + '.goterms.tsv', 'w') as oh:
    for gid in gene_order:
        if gid not in gene_go:
            continue
        for go in gene_go[gid]:
            oh.write(str(gid_map[gid]) + '\t' + go + '\t' + go_annot[go] + '\n')
