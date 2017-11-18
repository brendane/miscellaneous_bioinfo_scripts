#!/usr/bin/env python2.7
"""
    This script is meant to take a VCF file with just SNPs in which some
    SNPs may be listed as MNPs. It finds the actual position in the MNP
    that is variable and reports that, along with the variant ID.

    mnp2snp.py --output <output file> <input VCF>
"""

import argparse
import collections

import vcf

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--allele', action='store_true', default=False)
parser.add_argument('input')
args = parser.parse_args()

new_snps = collections.defaultdict(list)
rdr = vcf.Reader(filename=args.input)
for rec in rdr:
    chrom = rec.CHROM
    pos = rec.POS
    rs = rec.ID
    alleles = rec.alleles
    lens = map(len, alleles)
    all_same_length = reduce(lambda x, y: x and (y == lens[0]), lens[1:], True)
    if not all_same_length:
        print 'Indel at %s %i' % (chrom, pos)
        raise Exception('Indels in dataset')
    allele_length = lens[0]
    if allele_length == 1:
        for a in alleles[1:]:
            new_snps[rs].append((chrom, pos, rs, chrom, pos, str(a)))
    else:
        var_pos = []
        for i, b in enumerate(alleles[0]):
            variant = False
            for a in alleles[1:]:
                if str(a)[i] != b:
                    var_pos.append((pos + i, str(a)[i]))
                    break
        for vp in var_pos:
            new_snps[rs].append((chrom, vp[0], rs, chrom, pos, vp[1]))

with open(args.output, 'wb') as ohandle:
    if args.allele:
        ohandle.write('chrom\tpos\trs\torig_chrom\torig_pos\talt_allele\n')
    else:
        ohandle.write('chrom\tpos\trs\torig_chrom\torig_pos\n')
    for rs, info in sorted(new_snps.iteritems()):
        for nc, vp, rs, c, p, a in info:
            if args.allele:
                ohandle.write(nc + '\t' + str(vp) + '\t' + rs + '\t' +
                              c + '\t' + str(p) + '\t' + a + '\n')
            else:
                ohandle.write(nc + '\t' + str(vp) + '\t' + rs + '\t' +
                              c + '\t' + str(p) + '\n')
