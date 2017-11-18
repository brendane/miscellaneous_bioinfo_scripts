#!/usr/bin/env python2.7
"""
    Convert a VCF file into the format accepted by bugwas::lin_loc.

    vcf2bugwas.py --output <output prefix> <input file>

    Will treat heterozygous sites as missing data. Missing presence/
    absence data is treated as absent.

    Can perform imputation by choosing the homozygous allele with the
    greatest genotype likelihood (not if --pav is set); if there is no
    genotype likelihood information, it just sets the allele to
    the reference. Or can make missing data a third allele with --missallele.
"""

import argparse
import sys

import pysam

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--exclude')
parser.add_argument('--include')
parser.add_argument('--phenotype')
parser.add_argument('--pav', default=False, action='store_true')
parser.add_argument('--missallele', default=False, action='store_true')
parser.add_argument('--impute', default=False, action='store_true')
parser.add_argument('input')
args = parser.parse_args()

include = args.include
exclude = args.exclude
missing_extra_allele = args.missallele


vf = pysam.VariantFile(args.input)
samples = list(vf.header.samples)
if args.phenotype is not None:
    samples_ = []
    with open(args.phenotype, 'rb') as ih:
        for line in ih:
            s = line.strip().split('\t')[0]
            if s in samples:
                samples_.append(line.strip().split('\t')[0])
    samples = samples_

variants = []
i = 0
no_info = 0
with open(args.output + '.gen', 'wb') as oh:
    oh.write('ps\t' + '\t'.join(samples) + '\n')
    for rec in vf:
        if include is not None and rec.chrom != include:
            continue
        if exclude is not None and rec.chrom == exclude:
            continue
        oh.write(str(i))
        alleles = [rec.ref] + list(rec.alts)
        variants.append(rec.id)

        # Convert long alleles to one character - just make
        # them A, T, C, or G; doesn't matter whether they
        # the same as the actual allele
        if max(map(len, alleles)) > 1:
            alleles = ['A', 'T', 'C', 'G'][:len(alleles)]
        if missing_extra_allele:
            alleles = ['A', 'T', 'C', 'G'][:len(alleles)]
            missing_allele = ['A', 'T', 'C', 'G'][len(alleles)]
        else:
            missing_allele = 'N'

        for sample in samples:
            gt = rec.samples[sample]['GT']
            if gt[0] is None or (len(gt) > 1 and (gt[0] != gt[1])):
                if args.pav:
                    g = '0'
                else:
                    if args.impute:
                        gl = rec.samples[sample]['GL']
                        if gl[0] is None:
                            g = alleles[0] # Set to reference - no info
                            no_info += 1
                        elif gl[0] > gl[2]:
                            g = alleles[0]
                        else:
                            g = alleles[1]
                    else:
                        g = missing_allele
            else:
                if args.pav:
                    g = str(gt[0])
                else:
                    g = alleles[gt[0]]
            oh.write('\t' + g)
        oh.write('\n')
        i += 1

with open(args.output + '.variants', 'wb') as oh:
    for i, v in enumerate(variants):
        oh.write(str(i) + '\t' + v + '\n')

sys.stdout.write('Number of alleles set to reference with no info: %i\n' % no_info)
