#!/usr/bin/env python2.7
"""
    Switch the alleles for two samples in a VCF file. Useful for testing
    a recombination detection program.

    swap_vcf_alleles.py --output <output file> <vcf> <swap table>

    The swap table is tab-delimited with the following columns:
        chrom
        start (1-based, inclusive)
        end   (1-based, inclusive)
        strain1
        strain2
"""

import argparse

import vcf

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('vcf')
parser.add_argument('swap')
args = parser.parse_args()

to_swap = []
with open(args.swap, 'rb') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for row in rdr:
        to_swap.append((row[0], int(row[1]), int(row[2]), row[3], row[4]))
to_swap.sort(key=lambda x: (x[0], x[1], x[2]))

with open(args.output, 'wb') as oh:
    rdr = vcf.Reader(filename=args.vcf)
    wtr = vcf.Writer(oh, rdr)
    samples = rdr.samples
    for rec in rdr:
        for ts in to_swap:
            if ts[0] == rec.CHROM and rec.POS >= ts[1] and rec.POS <= ts[2]:
                i = samples.index(ts[3])
                j = samples.index(ts[4])
                tmp = rec.samples[i]
                rec.samples[i] = rec.samples[j]
                rec.samples[j] = swap
        wtr.write_record(rec)
    wtr.flush()
    wtr.close()
