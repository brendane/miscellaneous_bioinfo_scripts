#!/usr/bin/env python2.7
"""
    Get seeds from CD-HIT output.

    cdhit_seeds_fasta.py <cluster file> <fasta database>

    The output is a fasta file with records named by cluster (e.g.
    cluster-0, cluster-1)
"""

import gzip
import re
import sys

from Bio import SeqIO

cluster = None
start_of_cluster = False
cluster_genes = {}
ofun = open
if sys.argv[1].endswith('.gz'):
    ofun = gzip.open
with ofun(sys.argv[1], 'rb') as ih:
    for line in ih:
        if line.startswith('>'):
            start_of_cluster = True
            cluster = line.split()[1]
        elif start_of_cluster:
            start_of_cluster = False
            gene = re.sub('\.\.\.', '',
                          re.sub('^>', '', line.split()[2]))
            cluster_genes[gene] = cluster

for rec in SeqIO.parse(sys.argv[2], 'fasta'):
    g = rec.id
    if g in cluster_genes:
        c = cluster_genes[g]
        sys.stdout.write('>cluster-' + c + '\n')
        sys.stdout.write(str(rec.seq) + '\n')
