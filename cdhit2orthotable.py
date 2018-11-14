#!/usr/bin/env python2.7
"""
    Turn a clstr file from CD-Hit into a tsv file with gene names
    suitable for use in a database.

    cdhit2orthotable.py --output <output file> [--strain-map <strain map file]
        <clstr file>

    Assumes genes are named <strain>.<gene name>
"""

import argparse
import collections
import csv
import gzip
import re
import sys

def ddlist():
    return collections.defaultdict(list)

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--strain-map')
parser.add_argument('clstr')
args = parser.parse_args()

infile = args.clstr
outfile = args.output
strain_map_file = args.strain_map

strains = set()
clusters = collections.defaultdict(ddlist)

cluster_number = -1
ofun = open
if infile.endswith('.gz'):
    ofun = gzip.open
with ofun(infile, 'rb') as ih:
    for line in ih:
        line = line.strip()
        if line.startswith('>'):
            cluster_number = int(re.sub('>Cluster ', '', line))
        else:
            gene = line.split(' ')[1]
            strain = gene[1:].split('.')[0]
            if 'cluster-' in strain:
                continue
            strains.add(strain)
            clusters[cluster_number][strain].append(gene.split('.')[1])

strs = sorted(strains)

strain_map = {}
if strain_map_file is not None:
    with open(strain_map_file, 'rb') as ih:
        rdr = csv.reader(ih, delimiter='\t')
        for row in rdr:
            strain_map[row[0]] = row[1]
else:
    strain_map = {s:s for s in strs}

cls = sorted(clusters)
with open(outfile, 'wb') as oh:
    oh.write('cdhit_ortholog\tstrain\tgene\n')
    for c in cls:
        for strain, genes in clusters[c].iteritems():
            s = strain_map[strain]
            for gene in genes:
                oh.write('Cluster' + str(c) + '\t' + s + '\t' + gene + '\n')
