#!/usr/bin/env python3
"""
    Set alleles to missing if they do not have a large enough ratio of
    reads supporting one allele over another.

    filter_vcf_genotype_allele_ratio.py --output <output file> <vcf> <ratio>

    Gets rid of sites at which no individuals are genotyped.
"""

import argparse
import csv

import pysam

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('vcf')
parser.add_argument('ratio', type=float)
args = parser.parse_args()

with open(args.output, 'wb') as oh:
    rdr = pysam.VariantFile(args.vcf)
    wtr = pysam.VariantFile(args.output, mode='w', header=rdr.header)
    for rec in rdr:
        ungt = 0
        ratios = []
        for sample_name in rec.samples:
            sample = rec.samples[sample_name]
            ad = sample['AD']
            if ad[0] is not None:
                try:
                    r = max(ad) / min(ad)
                    ratios.append(str(r))
                except ZeroDivisionError:
                    r = args.ratio + 1
                    ratios.append('Inf')
                if r < args.ratio:
                    sample['GT'] = None
                    ungt += 1
            else:
                ungt += 1
                ratios.append('NaN')
        print(str(rec.chrom) + '\t' + str(rec.pos) + '\t' + '\t'.join(ratios))
        if ungt < len(rec.samples):
            wtr.write(rec)
    wtr.close()
