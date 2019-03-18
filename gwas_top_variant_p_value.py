#!/usr/bin/env python3
"""
    Get the top variant p-value for each gene in USDA1106

    Output:
        gene
        most significant p-value
        number of variants in gene (or LD groups?)
        mean of -log(p) values of all variants (or LD groups?) in gene
"""

import argparse
import collections
import csv
import gzip
import math
import re

def get_tags(gffstring):
    ret = {}
    for x in gffstring.split(';'):
        k, v = x.split('=')
        ret[k] = v
    return ret

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('gff')
parser.add_argument('gwas')
parser.add_argument('closest')
parser.add_argument('cdhit')
args = parser.parse_args()


locus_tags = {}
ofun = open
if args.gff.endswith('.gz'):
    ofun = gzip.open
with ofun(args.gff, 'r') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        ann = get_tags(row[8])
        if 'locus_tag' in ann and 'ID' in ann:
            locus_tags[ann['locus_tag']] = ann['ID']


annot = collections.defaultdict(set)
ofun = open
if args.closest.endswith('.gz'):
    ofun = gzip.open
with ofun(args.closest, 'r') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        if row[3] != '':
            annot[row[3]].add(re.sub('ID=(.+) :.+', '\\1', row[13]))

ofun = open
if args.cdhit.endswith('.gz'):
    ofun = gzip.open
cl = None
usda1106_gene = set()
denovo_pos = collections.defaultdict(set)
with ofun(args.cdhit, 'rt') as ih:
    for line in ih:
        l = line.strip()
        if l[0] == '>':
            if len(usda1106_gene) > 0:
                denovo_pos[cl] = usda1106_gene
            cl = l.split()[1]
            usda1106_gene = set()
            continue
        gs = l.split()[2]
        if gs.split('.')[0] == '>1106PB':
            usda1106_gene.add(locus_tags[gs.split('.')[1]])


top_vars = collections.defaultdict(lambda: 1)
all_vars = collections.defaultdict(lambda: 1)
n_vars = collections.defaultdict(lambda: 1)
ofun = open
if args.gwas.endswith('.gz'):
    ofun = gzip.open
with ofun(args.gwas, 'rt') as ih:
    rdr = csv.DictReader(ih, delimiter='\t')
    for row in rdr:
        p = float(row['p_lrt'])
        if len(annot[row['rs']]) > 0:
            genes = annot[row['rs']]
        else:
            continue
        for g in genes:
            if 'cluster' in g:
                gs = denovo_pos[re.sub('cluster_', '', g)]
                for gg in gs:
                    if p < top_vars[gg]:
                        top_vars[gg] = p
                        n_vars[gg] += 1
                        all_vars[gg].append(-math.log(p))
            else:
                if p < top_vars[g]:
                    top_vars[g] = p
                    n_vars[g] += 1
                    all_vars[g].append(-math.log(p))
with open(args.output, 'w') as oh:
    for g in top_vars:
        if g != '.':
            oh.write(g + '\t' + str(top_vars[g]) + '\t' + str(n_vars[g]) +
                     '\t' + str(sum(all_vars[g])/n_vars[g]) + '\n')
