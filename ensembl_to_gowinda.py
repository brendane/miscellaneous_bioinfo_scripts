#!/usr/bin/env python3
"""
    Given a gff file and a txt annotation file from ENSEMBL, produce
    input files for gowinda.

    ensembl_to_gowinda.py --output <prefix> <gff>
"""

import argparse
import collections
import csv
import re
import sys

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--go-domain')
parser.add_argument('--go-term-column', default='GO term accession')
parser.add_argument('--go-term-description', default='GO term name')
parser.add_argument('annot')
args = parser.parse_args()

domain = args.go_domain
gtc = args.go_term_column
gtd = args.go_term_description

go_annot = {}
gene_go = collections.defaultdict(set)
gene_names = {}
with open(args.annot, 'r') as ih:
    rdr = csv.DictReader(ih, delimiter='\t')
    for row in rdr:
        if domain is not None and row['GO domain'] != domain:
            continue
        gid = row['Gene stable ID']
        go = row[gtc]
        gene_go[go].add(gid)
        go_annot[go] = row[gtd]
        gene_names[gid] = row['Gene description']

with open(args.output, 'w') as oh:
    for go in go_annot:
        oh.write(go + '\t' + go_annot[go] + '\t' + ' '.join(gene_go[go]) + '\n')
