#!/usr/bin/env python
"""
    Convert bamutils RPKM output to VCF; meant for pan-genome.

    bamutils_rpkm2vcf.py
        --present <lower threshold for proportion of mean considered present>
        --absent <upper threshold for proportion of mean considered absent>
        <input directory>

    Output is to stdout.

    The input directory should have files named <strain>.<replicon>.tsv.
    
    Assumes input is in 0-based coordinates.

    --absent < --ambig < --greater
"""

#==============================================================================#

import argparse
import csv
import math
import os
import os.path as osp
import sys

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--absent', type=float)
parser.add_argument('--present', type=float)
parser.add_argument('indir')
args = parser.parse_args()

absent_threshold = args.absent
present_threshold = args.present
indir = args.indir

if absent_threshold >= present_threshold:
    raise Exception('--absent must be < --present')

strains_ = set()
replicons = set()
nfiles = 0
for fname in os.listdir(indir):
    if fname.endswith('.tsv'):
        strain, replicon, _ = fname.split('.')
        strains_.add(strain)
        replicons.add(replicon)
        nfiles += 1
if nfiles != len(strains_) * len(replicons):
    raise Exception('Some strain/replicons missing or extra files')

strains = sorted(strains_)
mean_depth = {}
strain_data = {}
gene_data = {}
genes_list = []

for strain in strains:
    strain_data[strain] = {}
    sum_depth = 0.
    n_genes = 0
    for replicon in replicons:
        fname = osp.join(indir, strain + '.' + replicon + '.tsv')
        with open(fname) as handle:
            lp = handle.tell()
            while handle.readline().startswith('#'):
                lp = handle.tell()
                continue
            handle.seek(lp)
            rdr = csv.DictReader(handle, delimiter='\t')
            for row in rdr:
                rpkm = float(row['RPKM'])
                gene = row['name']
                strain_data[strain][gene] = rpkm
                sum_depth += rpkm
                n_genes += 1
                if gene not in gene_data:
                    gene_data[gene] = (row['chrom'], row['start'])
                    genes_list.append(gene)
    mean_depth[strain] = sum_depth / n_genes

genotype_calls = {}
for strain in strains:
    genotype_calls[strain] = {}
    mean = mean_depth[strain]
    for gene in genes_list:
        ratio = strain_data[strain][gene] / mean
        if ratio >= present_threshold:
            call = 'T'
        elif ratio <= absent_threshold:
            call = 'A'
        else:
            call = 'N'
        genotype_calls[strain][gene] = call


sys.stdout.write('##fileformat=VCFv4.1\n')
sys.stdout.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT')
sys.stdout.write('\t' + '\t'.join(strains) + '\n')

for gene in genes_list:
    chrom = gene_data[gene][0]
    pos = gene_data[gene][1]
    var_id = gene
    gts = []
    for strain in strains:
        try:
            gts.append(genotype_calls[strain][gene])
        except:
            import pdb; pdb.set_trace()
    unique_gts = set(gts) - set('N')
    ref_gt = 'T'
    alt_gts = list(unique_gts - set(['T', 'N']))
    if len(alt_gts) == 0:
        alt_gts = ['.']

    sys.stdout.write(chrom + '\t' + str(pos) + '\t' + var_id + '\tT\t')
    sys.stdout.write(','.join(alt_gts) + '\t.\t.\t.\tGT')
    for gt in gts:
        if gt == 'N':
            sys.stdout.write('\t./.')
        elif gt == ref_gt:
            sys.stdout.write('\t0/0')
        else:
            g = str(alt_gts.index(gt) + 1)
            sys.stdout.write('\t' + g + '/' + g)

    sys.stdout.write('\n')
