#!/usr/bin/env python3
"""
    Given the Orthogroups.tsv file from OrthoFinder 2.3, extract
    unaligned files of sequences.

    extract_orthogroups_seqs.py [--ignore-missing] [--gene-list]
        <Orthogroups.tsv> <fasta> <output directory>
"""
 
import argparse
import csv

from Bio import SeqIO

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--ignore-missing', default=False, action='store_true')
parser.add_argument('--gene-lists')
parser.add_argument('orthogroups')
parser.add_argument('fasta')
parser.add_argument('output')
args = parser.parse_args()

seqs = SeqIO.index(args.fasta, 'fasta')

keep = set()
if args.gene_lists is not None:
    with open(args.gene_lists, 'rt') as ih:
        for line in ih:
            keep.add(line.strip())

with open(args.orthogroups) as ih:
    rdr = csv.DictReader(ih, delimiter='\t')
    for row in rdr:
        if len(keep) > 0 and row['Orthogroup'] not in keep:
            continue
        genes = []
        with open(args.output + '/' + row['Orthogroup'] + '.fasta', 'w') as oh:
            for field in rdr.fieldnames:
                if field == 'Orthogroup': continue
                for gene in row[field].split(', '):
                    g = field + '_' + gene
                    try:
                        SeqIO.write(seqs[g], oh, 'fasta')
                    except KeyError:
                        if args.ignore_missing:
                            pass
                        else:
                            raise Exception('%s not found in %s' % (g, args.fasta))
