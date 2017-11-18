#!/usr/bin/env python2.7
"""
    Convert a diploid VCF file to haploid while preserving as much
    information as possible. Note that it's possible that some
    information should not be preserved under some circumstances.

    vcf_diploid2haploid.py <input file> <output file>

    Will treat heterozygotes as missing data.
"""

import argparse

import pysam

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('input')
parser.add_argument('output')
args = parser.parse_args()

with open(args.output, 'wb') as oh:
    vf = pysam.VariantFile(args.input)
    header = vf.header
    oh.write(str(header))
    for rec in vf:
        # Printing:
        # rec.chrom, rec.pos, rec.id, rec.ref, ','.join(rec.alts), rec.qual,
        # '.', '.', rec.format, samples...
        ks = rec.format.keys()
        recid = rec.id
        if recid is None:
            recid = rec.chrom + '-' + str(rec.pos)
        oh.write(rec.chrom + '\t' + str(rec.pos) + '\t' + recid + '\t' +
                 rec.ref + '\t' + ','.join(rec.alts) + '\t' + str(rec.qual) +
                 '\t.\t.\t' + ':'.join(ks))
        for sample in rec.samples:
            s_rec = rec.samples[sample]
            new_values = {}
            for k in ks:
                v = s_rec[k]
                if type(v) == int:
                    new_values[k] = v
                elif k == 'GT':
                    if len(v) > 1 and v[0] != v[1]:
                        new_values[k] = '.'
                    elif v[0] is None:
                        new_values[k] = '.'
                    else:
                        new_values[k] = v[0]
                elif v is None:
                    new_values[k] = '.'
                elif len(v) == 1:
                    new_values[k] = v[0]
                else:
                    if len(v) == 2:
                        new_values[k] = v[0]
                    if len(v) == 3:
                        new_values[k] = str(v[0]) + ',' + str(v[2])
            oh.write('\t' + ':'.join(str(new_values[k]) for k in ks))
        oh.write('\n')
