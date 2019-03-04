#!/usr/bin/env python3
"""
    Combine ENSEMBL annotation file (tab-delimited) with NCBI gff file.
    Made specifically for dolphin project.
"""

import argparse
import collections
import csv
import re
import sys

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('ensembl')
parser.add_argument('gff')
args = parser.parse_args()

gff_info = {}
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
            gid = re.sub('GeneID:', '',
                         {x.split('=')[0]:x.split('=')[1] for x in annot.split(';')}['Dbxref'])
            gff_info[gid] = (chrom, start, end)
        else:
            continue

nogid = set()
go_annot = {}
gene_go = collections.defaultdict(set)
errors = set()
with open(args.ensembl, 'r') as ih:
    rdr = csv.DictReader(ih, delimiter='\t')
    for row in rdr:
        egid = row['Gene stable ID']
        estart = int(row['Gene start (bp)'])
        eend = int(row['Gene end (bp)'])
        desc = row['Gene description']
        gid = row['NCBI gene ID']
        echrom = row['Chromosome/scaffold name']
        go = row['GO term accession']
        go_desc = row['GO term name']
        if gid not in gff_info:
            if gid == '': nogid.add(egid)
            errors.add(('%s in ensembl but not ncbi' % gid))
            continue
        if gff_info[gid][1] != estart:
            errors.add(('%s positions do not match (%i, %i) in ncbi and (%i, %i) in ens.'
                        %(gid, gff_info[gid][1], gff_info[gid][2], estart, eend)))
        gene_go[gid].add(go)
        go_annot[go] = go_desc

for error in errors:
    sys.stdout.write(error + '\n')
print('No NCBI gene ID: %i' % len(nogid))

print('NCBI genes annotated: %i' % len(gene_go))
print('Total NCBI genes in gff %i' %len(gff_info))
