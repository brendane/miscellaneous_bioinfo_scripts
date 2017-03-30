#!/usr/bin/env python2
"""
    Identify LD clusters associated with top SNPs.

    top_snps_ld_clusters.py <input file> <output prefix> <N>
"""

#==============================================================================#

import argparse
import collections
import csv

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('input')
parser.add_argument('outpre')
parser.add_argument('N', type=int)
args = parser.parse_args()

N = args.N

data = {}
with open(args.input, 'rb') as handle:
    rdr = csv.DictReader(handle, delimiter='\t')
    for row in rdr:
        data[int(row['rank'])] = row

genes = collections.defaultdict(set)
groups = collections.defaultdict(set)
rank = 1
while rank <= N:
    try:
        gene = data[rank]['gene']
    except KeyError:
        print "WARNING: gene for rank %i is missing" % rank
        rank += 1
        continue
    group = int(data[rank]['group'])
    if gene != '.':
        genes[gene].add(group)
    groups[group].add(gene)
    #if len(genes) == N:
    #    break
    rank += 1

with open(args.outpre + '.genes.tsv', 'wb') as out:
    out.write('gene\tgroups\n')
    for gene, grps in genes.iteritems():
        out.write(gene + '\t' + ';;'.join(str(g) for g in grps) + '\n')

with open(args.outpre + '.groups.tsv', 'wb') as out:
    out.write('group\tgenes\n')
    for group, genes in groups.iteritems():
        out.write(str(group) + '\t' + ';;'.join(genes) + '\n')

