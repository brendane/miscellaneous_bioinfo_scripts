#!/usr/bin/env python3
"""
    Given the output of extract_orthosets_from_orthofinder.py,
    extract unaligned files of sequences.

    extract_orthosets_seqs.py <orthosets> <group> <fasta> <output directory>
"""
 
import argparse
import csv

from Bio import SeqIO

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('orthosets')
parser.add_argument('group')
parser.add_argument('fasta')
parser.add_argument('output')
args = parser.parse_args()

seqs = SeqIO.index(args.fasta, 'fasta')

with open(args.orthosets) as ih:
    rdr = csv.DictReader(ih, delimiter='\t')
    for row in rdr:
        if row['taxon'] != args.group:
            continue
        with open(args.output + '/' + row['subset'] + '.fasta', 'w') as oh:
            for g in row['genes'].split(','):
                try:
                    SeqIO.write(seqs[g], oh, 'fasta')
                except KeyError:
                    raise Exception('%s not found in %s' % (g, args.fasta))
