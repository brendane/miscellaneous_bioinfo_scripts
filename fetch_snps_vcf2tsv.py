#!/usr/bin/env python2.7
"""
    Get variants from a VCF file using a bed file and format into tsv.
    
    fetch_snps_vcf2tsv.py --output <output file> [--strains <file with strains>]
        <vcf> <bed>

    Requires pyvcf.

    Changes heterozygotes to NaN.
"""

import argparse

import vcf

parser = argparse.ArgumentParser()
parser.add_argument('--output')
parser.add_argument('--strains')
parser.add_argument('vcf_file')
parser.add_argument('bed_file')
args = parser.parse_args()

outfile = args.output

strains = []
if args.strains is not None:
    with open(args.strains, 'rb') as ih:
        for line in ih:
            strains.append(line.strip())

regions = []
with open(args.bed_file, 'rb') as ih:
    for line in ih:
        fields = line.strip().split('\t')
        regions.append((fields[0], int(fields[1]), int(fields[2]), fields[3]))

with open(args.output, 'wb') as oh:
    with open(args.vcf_file, 'rb') as ih:
        v = vcf.VCFReader(ih)
        if len(strains) == 0:
            strains = v.samples
        oh.write('chrom\tpos\tregion\t' + '\t'.join(strains) + '\n')
        for r in regions:
            region_vars = list(v.fetch(r[0], r[1], r[2]))
            for var in region_vars:
                oh.write(var.CHROM + '\t' + str(var.POS) + '\t' + r[3])
                for strain in strains:
                    gt = var.genotype(strain).gt_alleles
                    if len(gt) > 1 and gt[0] != gt[1]:
                        gt = 'NaN'
                    else:
                        gt = str(gt[0])
                    oh.write('\t' + gt)
                oh.write('\n')
