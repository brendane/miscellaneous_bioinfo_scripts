#!/usr/bin/env python3
"""
    Given a completed OrthoFinder run and a list of focal taxa, extract
    a list of genes for each ortholog.

    extract_orthosets_from_orthofinder.py --output <output file>
        [--min-support (0.0)]
        <duplications.csv> <orthogroups.csv> <taxa=focal strain list> ...

    Finds smallest groups of genes that contain all of the focal strains
    present in each orthogroup. Can find more than one set per orthogroup.

    Assumes that strain names have one underscore.
"""

import argparse
import csv
import itertools
import sys

def find_ortholog_sets(strains, strain_sets, gene_sets):
    ## Initial set up
    subs = []
    gene_subs = []
    strain_subs = []
    genes_accounted_for = set()
    strains_to_include = set.intersection(strains, strain_sets[0])
    genes_to_include = set()
    for g in gene_sets[0]:
        if '_'.join(g.split('_')[:2]) in strains:
            genes_to_include.add(g)
    ## Keep identifying subsets until all the genes are accounted for
    while len(genes_to_include - genes_accounted_for) > 0:
        min_size = len(genes_to_include) + 1
        candidate = 0
        ## Find the smallest set of genes that includes all the target
        ## strains. If there is more than one set, this will (I think!)
        ## choose the first because `n_focal_genes_in_set < min_size`
        ## instead of `n_focal_genes_in_set <= min_size`
        for i, (gs, ss) in enumerate(zip(gene_sets, strain_sets)):
            if len(set.intersection(gs, genes_accounted_for)) > 0:
                ## Overlaps existing subset, do not use
                continue
            if len(strains_to_include - ss) != 0:
                ## This set does not include all the focal strains
                continue
            n_focal_genes_in_set = len(set.intersection(gs, genes_to_include))
            if n_focal_genes_in_set < min_size:
                candidate = i
                min_size = n_focal_genes_in_set
        ## Add the subset just identified to the list, then subtract
        ## the genes from the list of genes to include in the next
        ## iteration. Update the list of strains using the updated
        ## gene list.
        subs.append(candidate)
        gene_subs.append(set.intersection(genes_to_include, gene_sets[candidate]))
        strain_subs.append(strains_to_include)
        genes_accounted_for.update(gene_sets[candidate])
        genes_to_include -= genes_accounted_for
        strains_to_include = {'_'.join(g.split('_')[:2]) for g in genes_to_include}
    return subs, gene_subs, strain_subs


parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--min-support', type=float, default=0)
parser.add_argument('dups')
parser.add_argument('orthos')
parser.add_argument('strainlists', nargs='+')
args = parser.parse_args()

## List of strains for each taxon
taxa = {}
for sl in args.strainlists:
    taxon, fname = sl.split('=')
    taxa[taxon] = set()
    with open(fname, 'r') as ih:
        for line in ih:
            taxa[taxon].add(line.strip())

## Genes and strains in each orthogroup
strains_in_ogs = {}
genes_in_ogs = {}
with open(args.orthos, 'r', newline='') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    strains = next(rdr)[1:]
    for row in rdr:
        og = row[0]
        strains_in_ogs[og] = set()
        genes_in_ogs[og] = set()
        for i, gs in enumerate(row[1:]):
            if gs != '':
                strains_in_ogs[og].add(strains[i])
                genes_in_ogs[og].update(strains[i] + '_' + g for g in gs.split(', '))

## Go through duplications file
oh = sys.stdout
ogs_in_dups = set()
oh.write('orthogroup\tsubset\ttaxon\tsingle_copy\tstrains\tgenes\n')
with open(args.dups, 'r', newline='') as ih:
    rdr = csv.DictReader(ih, delimiter='\t')

    for og, rows in itertools.groupby(rdr, lambda x: x['Orthogroup']):
        ogs_in_dups.add(og)
        gene_sets = [genes_in_ogs[og]]     # First item is just the whole orthogroup
        strain_sets = [strains_in_ogs[og]] # First item is just the whole orthogroup

        ## Get all sets of genes and strains for this orthogroup
        for row in rows:
            if float(row['Support']) < args.min_support:
                continue
            gs1 = set(row['Genes 1'].split(', '))
            gs2 = set(row['Genes 2'].split(', '))
            gene_sets.append(gs1)
            strain_sets.append({'_'.join(g.split('_')[:2]) for g in gs1})
            gene_sets.append(gs2)
            strain_sets.append({'_'.join(g.split('_')[:2]) for g in gs2})

        ## Identify gene subsets for each taxon
        for taxon in taxa:
            strains = taxa[taxon]
            for sub, g, s in zip(*find_ortholog_sets(strains, strain_sets, gene_sets)):
                oh.write(og + '\t' +
                         og + '.' + str(sub) + '\t' +
                         taxon + '\t' +
                         str(int(len(s) == len(g))) + '\t' +
                         ','.join(s) + '\t' +
                         ','.join(g) + '\n')


## Orthogroups not in duplications file report as they are
for og in genes_in_ogs:
    if og in ogs_in_dups:
        continue
    for taxon in taxa:
        s = set.intersection(taxa[taxon], strains_in_ogs[og])
        g = [g for g in genes_in_ogs[og] if '_'.join(g.split('_')[:2]) in taxa[taxon]]
        oh.write(og + '\t' +
                 og + '.0\t' +
                 taxon + '\t' +
                 str(int(len(s) == len(g))) + '\t' +
                 ','.join(s) + '\t' +
                 ','.join(g) + '\n')
