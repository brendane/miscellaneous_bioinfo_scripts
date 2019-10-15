#!/usr/bin/env python3
"""
    Combine a whole bunch of OrthoFinder-related output files along with
    other information.

    combine_orthofinder_dup_annotation_orthosets.py --output <output file>
        --orthosets <orthoset table from extract_orthosets_from_orthofinder_trees.py>
        --annotation <annotation file from combine_orthotable_and_annotations.py>
        --gene-annotation <annotation file with columns for strain, gene, product, partial, and pseudogene status>
        --dups <directory with .dups files (optional) containing identical duplicates>
        --species <species names file: one column for strain, one column for species>
        --no-pseudo <flag to exclude pseudogenes>
        --exclude <file with list of strain names to exclude>

    Assumes single underscore separates strain and gene names and that
    the strain name is included except in the dups files. Also, will
    combine annotation information across all orthosets in an orthogroup.
"""

import argparse
import csv
import os
import os.path as osp

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--orthosets')
parser.add_argument('--annotation')
parser.add_argument('--gene-annotation')
parser.add_argument('--dups')
parser.add_argument('--species')
parser.add_argument('--no-pseudo', default=False, action='store_true')
parser.add_argument('--exclude')
args = parser.parse_args()

dups = {}
if args.dups is not None:
    for fname in os.listdir(args.dups):
        if not fname.endswith('dups'):
            continue
        strain = osp.splitext(osp.basename(fname))[0]
        with open(osp.join(args.dups, fname), 'rt') as ih:
            for line in ih:
                genes = line.strip().split('\t')
                if len(genes) > 1:
                    dups[strain + '_' + genes[0]] = \
                            {strain + '_' + g for g in genes[1:]}


pseudo = set()
with open(args.gene_annotation, 'rt') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        if row[4] == '1':
            pseudo.add(row[0])

annot = {}
with open(args.annotation, 'rt') as ih:
    rdr = csv.DictReader(ih, delimiter='\t')
    for row in rdr:
        og = row['orthogroup']
        if row['gene'] == '':
            genes = set()
        else:
            genes = {g for g in row['gene'].split(',')}
        if row['descriptions'] == '':
            descs = set()
        else:
            descs = {d for d in row['descriptions'].split('; ')}
        if og in annot:
            annot[og][0].update(genes)
            annot[og][1].update(descs)
        else:
            annot[og] = [genes, descs]

spp = {}
gen = {}
with open(args.species, 'rt') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        spp[row[0]] = row[1]
        gen[row[0]] = row[1].split()[0]

exclude = set()
if args.exclude is not None:
    with open(args.exclude, 'rt') as ih:
        for line in ih:
            exclude.add(line.strip())

with open(args.output, 'wt') as oh:
    oh.write('orthogroup\torthoset\tn_genera\tn_species\tn_strains\tn_genes\t' +
             'gene_names\tgene_products\t' +
             'strains\tgenes\n')
    with open(args.orthosets, 'rt') as ih:
        rdr = csv.DictReader(ih, delimiter='\t')
        for row in rdr:
            og = row['orthogroup']
            os = row['subset']
            genes = []
            strains = set()
            species = set()
            genera = set()
            for g in row['genes'].split(','):
                if args.no_pseudo and g in pseudo:
                    continue
                strain = g.split('_')[0]
                if strain in exclude:
                    continue
                genes.append(g)
                if g in dups:
                    for d in dups[g]:
                        genes.append(d)
                strains.add(strain)
                species.add(spp[strain])
                genera.add(gen[strain])
            if len(genes) == 0:
                continue
            result =[og,
                     os,
                     len(genera),
                     len(species),
                     len(strains),
                     len(genes),
                     ','.join(n for n in annot[og][0] if n != ''),
                     '; '.join(d for d in annot[og][1] if d != ''),
                     ','.join(strains),
                     ','.join(genes)]
            oh.write('\t'.join(str(x) for x in result) + '\n')
