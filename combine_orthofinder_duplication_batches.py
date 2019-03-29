#!/usr/bin/env python3
"""
    Combine output files from my modified OrthoFinder tree reconciliation
    script (trees2ologs_modified.py). Also, add identical duplicates that
    were removed before running blast, and add in single-gene orthogroups.

    combine_orthofinder_duplication_batches.py --output <output file>
        <duplications.csv dir> <dups dir> <Orthogroups.csv file>

    Naming assumptions:
    - *.dups files are named <strain>.dups
    - Genes in the *.dups files are just the gene name, no strain
    - Duplications.*.csv files are named thus, with * = batch number
    - Genes in the "Genes 1" and "Genes 2" are named <strain>_<gene name>
    - Genes in Orthogroups.csv do not include the strain name.
"""

import argparse
import csv
import os
import re

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('ortho_dups_dir')
parser.add_argument('ident_dups_dir')
parser.add_argument('orthogroups')
args = parser.parse_args()

ident_dups = {}
all_genes = set() ## All genes listed in the first column of the *.dups files
for fname in os.listdir(args.ident_dups_dir):
    if fname.endswith('.dups'):
        with open(args.ident_dups_dir + '/' + fname, 'r') as ih:
            strain = re.sub('.dups', '', fname)
            for line in ih:
                genes = line.strip().split('\t')
                all_genes.add(strain + '_' + genes[0])
                if len(genes) > 1:
                    ident_dups[strain + '_' + genes[0]] = [strain + '_' + g for g in genes[1:]]


first = True
OGs_done = set()
genes_done = set()
with open(args.output, 'w') as oh:
    for fname in os.listdir(args.ortho_dups_dir):
        if fname.endswith('.csv') and fname.startswith('Duplications.'):
            with open(args.ortho_dups_dir + '/' + fname, 'r') as ih:
                rdr = csv.DictReader(ih, delimiter='\t')
                if first:
                    oh.write('\t'.join(rdr.fieldnames) + '\n')
                    first = False
                for row in rdr:
                    OGs_done.add(row['Orthogroup'])
                    oh.write('\t'.join([row['Orthogroup'], row['Species Tree Node'],
                                        row['Gene Tree Node'], row['Support'],
                                        row['Type']]))
                    for split in ['Genes 1', 'Genes 2']:
                        genes = row[split].strip().split(', ')
                        genes_with_dups = []
                        for gene in genes:
                            genes_with_dups.append(gene)
                            genes_done.add(gene)
                            if gene in ident_dups:
                                for g in ident_dups[gene]:
                                    genes_with_dups.append(g)
                                    genes_done.add(g)
                        oh.write('\t' + ', '.join(genes_with_dups))
                    oh.write('\n')
    with open(args.orthogroups, 'r') as ih:
        rdr = csv.DictReader(ih, delimiter='\t')
        strains_ortho = rdr.fieldnames[1:]
        for row in rdr:
            og = row['']
            if og not in OGs_done:
                OGs_done.add(og)
                genes_with_dups_and_strain = []
                for strain in strains_ortho:
                    genes = row[strain].strip().split(', ')
                    if len(genes) == 1 and genes[0] == '':
                        continue
                    for gene in genes:
                        gene_with_strain = strain + '_' + gene
                        genes_with_dups_and_strain.append(gene_with_strain)
                        genes_done.add(gene_with_strain)
                        if gene_with_strain in ident_dups:
                            for g in ident_dups[gene_with_strain]:
                                genes_with_dups_and_strain.append(g)
                                genes_done.add(g)
                oh.write('\t'.join([row[''], 'NoName',
                                    'NoName', '1.0', 'NoDup',
                                    ', '.join(genes_with_dups_and_strain), ''])
                         + '\n')
    ogi = len(OGs_done)
    for gene in all_genes:
        if gene in genes_done:
            continue
        genes_done.add(gene)
        gene_out = []
        gene_out.append(gene)
        if gene in ident_dups:
            for g in ident_dups[gene]:
                if g in genes_done:
                    raise Exception('Identical duplicate (%s) listed without parent (%s)'
                                    % (gene, g))
                gene_out.append(g)
                genes_done.add(g)
        oh.write('\t'.join(['OG' + str(ogi).zfill(7), 'NoName',
                            'NoName', '1.0', 'SingleGroup',
                            ', '.join(gene_out), ''])
                 + '\n')
        ogi += 1
