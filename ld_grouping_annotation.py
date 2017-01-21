#!/usr/bin/env python2
"""
    Given a file with SNPs grouped by LD, extract information on the
    genes and p-values for each group.

    ld_grouping_annotation.py <ld file> <annotation file>

    Assumes that the annotation file is sorted by p-value, with lowest
    first.
"""

#==============================================================================#

import collections
import csv
import re
import sys

csv.field_size_limit(sys.maxsize)

#==============================================================================#

ld_file = sys.argv[1]
ann_file = sys.argv[2]

ld_groups = collections.defaultdict(set)
with open(ld_file, 'rb') as handle:
    rdr = csv.reader(handle, delimiter='\t')
    for row in rdr:
        grp = int(row[4])
        for snp in row[5].split(','):
            chrom = re.sub('-.+', '', snp)
            pos = int(re.sub('.+-', '', snp))
            ld_groups[grp].add((chrom, pos))

annot = {}
with open(ann_file, 'rb') as handle:
    rdr = csv.reader(handle, delimiter='\t')
    for rank, row in enumerate(rdr):
        chrom = row[0]
        pos = int(row[2])
        gene = row[12]
        pval = float(row[3])
        annot[(chrom, pos)] = (pval, gene, rank+1)

sys.stdout.write('group\tchrom\tpos\trank\tp\tgene\n')
for grp in sorted(ld_groups):
    info = []
    genes = set()
    pvals = []
    ranks = set()
    for snp in ld_groups[grp]:
        if snp not in annot:
            continue
        pval, gene, rank = annot[snp]
        info.append((rank, pval, gene, snp[0], snp[1]))
    info.sort()
    for i in info:
        sys.stdout.write(str(grp) + '\t' + i[3] + '\t' +
                         str(i[4]) + '\t' + str(i[0]) +  '\t' +
                         str(i[1]) + '\t' + str(i[2]) + '\n')
