#!/usr/bin/env python3
"""
    Set alleles to missing if they do not have a large enough ratio of
    reads supporting one allele over another.

    filter_vcf_genotype_allele_ratio.py --output <output file> <vcf> <ratio>
"""

import argparse
import csv

import vcf

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('vcf')
parser.add_argument('ratio', type=float)
args = parser.parse_args()

with open(args.output, 'wt') as oh:
    rdr = vcf.Reader(filename=args.vcf)
    wtr = vcf.Writer(oh, rdr)
    samples = rdr.samples
    for rec in rdr:
        for sample in rec.samples:
            ad = sample.data.AD
            if ad is not None:
                try:
                    r = max(ad) / min(ad)
                except ZeroDivisionError:
                    r = args.ratio + 1
                if r < args.ratio:
                    ## TODO: Figure out how to set genotype to missing
                    sample.gt_alleles = ['.']
        wtr.write_record(rec)
    wtr.flush()
    wtr.close()
