#!/usr/bin/env python3
"""
    Combine *.assoc.txt.gz files from GEMMA and annotations from bedtools
    and custom python scripts (*.closest.genes.tsv). Produce output with
    one row per gene.

    combine_gemma_annotation_ld_groups.py <assoc.txt.gz> <closest.genes.tsv>

    Assumes annotation looks like XXX=geneID : annotation, where "XXX"
    is usually "ID". Also, extract the "p_lrt" column from the GEMMA
    output.

    Also, can translate between cluster_XXXX from CD-Hit by looking for a
    ref_id=... component of the annotation. This is specific to a
    particular dataset. If no ref_id is found, the gene will be skipped.
    Thus, this script is meant for working with a single reference genome.
"""

import argparse
import collections
import csv
import gzip
import math
import re
import sys

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('assoc')
parser.add_argument('closest')
args = parser.parse_args()

genes = collections.defaultdict(set)
if args.closest.endswith('.gz'):
    op = gzip.open
else:
    op = open
with op(args.closest, 'rt') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        grp = row[3]
        if row[-1] == '' or row[-1] == '.':
            continue
        if row[-1].startswith('ID=cluster_'):
            gene = re.sub('.+ref_id=', '', row[-1])[:-1]
            if gene.startswith('ID=cluster'):
                # No reference ID, skip this one
                continue
        else:
            gene = row[-1].strip().split(' : ')[0].split('=')[1]
        genes[grp].add(gene)

pvals = {}
if args.assoc.endswith('.gz'):
    op = gzip.open
else:
    op = open
with op(args.assoc, 'rt') as ih:
    rdr = csv.DictReader(ih, delimiter='\t')
    for row in rdr:
        pvals[row['rs']] = float(row['p_lrt'])

oh = sys.stdout
oh.write('LD_group\tgene\tp\tneg_log_p\n')
for grp, p in pvals.items():
    if len(genes[grp]) > 0:
        for g in genes[grp]:
            oh.writelines([grp, '\t', g, '\t', str(p), '\t',
                           str(-math.log(p)), '\n'])
