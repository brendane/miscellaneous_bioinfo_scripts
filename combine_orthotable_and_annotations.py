#!/usr/bin/env python3
"""
    Given orthosets extracted from OrthoFinder and a table of annotation
    information, provide a summary of gene names and descriptions.

    combine_orthotable_and_annotations.py <orthotable> <annotation>
"""

import csv
import sys

annot = {}
with open(sys.argv[2], 'rt') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        annot[row[0]] = (row[1], row[2])

with open(sys.argv[1], 'rt') as ih:
    rdr = csv.DictReader(ih, delimiter='\t')
    for row in rdr:
        genes = row['genes'].split(',')
        names = set()
        descriptions = []
        for gene in genes:
            a = annot[gene]
            names.add(a[0])
            descriptions.append(a[1])
        consensus_description = ???
        sys.stdout.writelines([row['subset'], '\t', ','.join(names), '\t',
                               consensus_description, '\n'])
