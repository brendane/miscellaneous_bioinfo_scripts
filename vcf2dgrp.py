#!/usr/bin/env python2
"""
    Transform a VCF file into a DGRP file.

    vcf2dgrp.py --output <output file> <input file> <target replicon>

    OPTIONS

        --subset        A bed file with positions to include, default is
                        to include everything
        --min-gt        Minimum genotyping rate (0)
"""

#==============================================================================#

import argparse
import re

import vcf

#==============================================================================#

def ambig(g):
    return {('A','A'):'A', ('T','T'):'T', ('C','C'):'C', ('G','G'):'G',
            ('A','C'):'M', ('A','G'):'R', ('A','T'):'W', ('C','G'):'S',
            ('C','T'):'Y', ('G','T'):'K',
            ('C','A'):'M', ('G','A'):'R', ('T','A'):'W', ('G','C'):'S',
            ('T','C'):'Y', ('T','G'):'K'}[g]

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--subset')
parser.add_argument('--min-gt', default=0.0, type=float)
parser.add_argument('input')
parser.add_argument('replicon')
args = parser.parse_args()

if args.subset is None:
    positions = None
else:
    positions = set()
    with open(args.subset, 'rb') as handle:
        for line in handle:
            contig, _, pos = line.strip().split()[:3]
            if contig not in args.replicon:
                continue
            positions.add((contig, pos))

with open(args.output, 'wb') as out:
    rdr = vcf.Reader(filename=args.input)
    out.write(args.replicon + ',Ref,' + ','.join(rdr.samples) + '\n')
    for rec in rdr:
        if (positions is not None) and ((rec.CHROM, str(rec.POS)) not in positions):
            continue
        if not rec.is_snp:
            continue
        gts = []
        missed = 0
        for s in rec.samples:
            gt = re.split('\||/', s['GT'])
            if gt[0] == '.':
                g = 'N'
                missed += 1
            else:
                g = tuple(str(rec.alleles[int(x)]).upper() for x in gt)
                g = ambig(g)
            gts.append(g)
        if 1 - (missed / float(len(gts))) < args.min_gt:
            continue
        out.write(str(rec.POS))
        out.write(',' + rec.REF + ',' + ','.join(str(g) for g in gts) + '\n')
