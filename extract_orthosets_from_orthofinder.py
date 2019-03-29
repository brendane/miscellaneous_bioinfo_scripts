#!/usr/bin/env python3
"""
    Given a completed OrthoFinder run and a list of focal taxa, extract
    a list of genes for each ortholog.

    extract_orthosets_from_orthofinder.py --output <output file>
        [--min-support (0.0)]
        <duplications.csv> <list of strain suffixes> <taxa=focal strain list> ...

    Finds smallest groups of genes that contain all of the focal strains
    present in each orthogroup. Can find more than one set per orthogroup.
"""

import argparse
import csv
import itertools
import sys

def genes2strains(genes, strains, suffixes):
    """ Because this function must work on all strains, not just those
    included in the analysis, it can return no match."""
    ret = []
    for gene in genes:
        found = False
        for suff in suffixes:
            if suff in gene:
                attempt = gene.split(suff)[0] + suff
                if attempt in strains:
                    found = True
                    ret.append(attempt)
                    break
        if not found:
            ret.append(None)
    return ret

def find_ortholog_sets(strains, strain_sets, gene_sets, all_taxa, suffs):
    ## Initial set up
    subs = []
    gene_subs = []
    strain_subs = []
    genes_accounted_for = set()
    strains_to_include = set.intersection(strains, strain_sets[0])
    genes_to_include = set()
    for strain, g in zip(genes2strains(gene_sets[0], all_taxa, suffs), gene_sets[0]):
        if strain in strains:
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
        strains_to_include = set(genes2strains(genes_to_include, all_taxa, suffs))
    return subs, gene_subs, strain_subs


parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--min-support', type=float, default=0.)
parser.add_argument('dups')
parser.add_argument('suffixes')
parser.add_argument('strainlists', nargs='+')
args = parser.parse_args()

suffs = []
with open(args.suffixes, 'r') as ih:
    for line in ih:
        suffs.append(line.strip())

## List of strains for each taxon
taxa = {}
all_taxa = set()
for sl in args.strainlists:
    taxon, fname = sl.split('=')
    taxa[taxon] = set()
    with open(fname, 'r') as ih:
        for line in ih:
            taxa[taxon].add(line.strip())
            all_taxa.add(line.strip())


## Go through duplications file
strains_observed = set()
with open(args.output, 'w') as oh:
    oh.write('orthogroup\tsubset\ttaxon\tsingle_copy\tn_strains\tcore\tstrains\tgenes\n')
    with open(args.dups, 'r', newline='') as ih:
        rdr = csv.DictReader(ih, delimiter='\t')

        for og, rows in itertools.groupby(rdr, lambda x: x['Orthogroup']):
            gene_sets = []
            strain_sets = []

            ## Get all sets of genes and strains for this orthogroup
            first = True
            for row in rows:
                if first:
                    gene_sets.append(set(row['Genes 1'].split(', ') + row['Genes 2'].split(', ')))
                    strain_sets.append(set(genes2strains(gene_sets[0], all_taxa, suffs)))
                    first = False
                if float(row['Support']) < args.min_support:
                    continue
                gs1 = set(row['Genes 1'].split(', '))
                gs2 = set(row['Genes 2'].split(', '))
                gene_sets.append(gs1)
                strain_sets.append(set(genes2strains(gs1, all_taxa, suffs)))
                gene_sets.append(gs2)
                strain_sets.append(set(genes2strains(gs2, all_taxa, suffs)))

            ## Identify gene subsets for each taxon
            for taxon in taxa:
                strains = taxa[taxon]
                for sub, g, s in zip(*find_ortholog_sets(strains, strain_sets, gene_sets, all_taxa, suffs)):
                    strains_observed.update(s)
                    oh.write(og + '\t' +
                             og + '.' + str(sub) + '\t' +
                             taxon + '\t' +
                             str(int(len(s) == len(g))) + '\t' +
                             str(len(s)) + '\t' +
                             str(int(len(s) == len(strains))) + '\t' +
                             ','.join(s) + '\t' +
                             ','.join(g) + '\n')

print('Strains not found:')
for taxon, strains in taxa.items():
    print(taxon + ':')
    print('\t' + ', '.join(strains - strains_observed))
