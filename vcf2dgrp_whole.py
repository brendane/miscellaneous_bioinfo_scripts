#!/usr/bin/env python2
"""
    Transform a VCF file into a DGRP file. Targets entire genome.

    vcf2dgrp.py --reference <reference genome> --output <output file> <input file>

    OPTIONS

        --subset        A bed file with positions to include, default is
                        to include everything
        --min-gt        Minimum genotyping rate (0)
        --het-miss      Set heterozygous sites to missing
        --pad           Size of pad between contigs (1000)
        --reference-name Name of reference strain; include reference as
                        a sample

    UPDATE 02 Oct 2018: added option to include reference alleles as
    a sample.
"""

#==============================================================================#

import argparse
import re

from Bio import SeqIO
import vcf

#==============================================================================#

def ambig(g):
    return {('A','A'):'A', ('T','T'):'T', ('C','C'):'C', ('G','G'):'G',
            ('A','C'):'M', ('A','G'):'R', ('A','T'):'W', ('C','G'):'S',
            ('C','T'):'Y', ('G','T'):'K',
            ('C','A'):'M', ('G','A'):'R', ('T','A'):'W', ('G','C'):'S',
            ('T','C'):'Y', ('T','G'):'K'}[g]

def ambighm(g):
    return {('A','A'):'A', ('T','T'):'T', ('C','C'):'C', ('G','G'):'G',
            ('A','C'):'N', ('A','G'):'N', ('A','T'):'N', ('C','G'):'N',
            ('C','T'):'N', ('G','T'):'N',
            ('C','A'):'N', ('G','A'):'N', ('T','A'):'N', ('G','C'):'N',
            ('T','C'):'N', ('T','G'):'N'}[g]

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--reference')
parser.add_argument('--subset')
parser.add_argument('--min-gt', default=0.0, type=float)
parser.add_argument('--het-miss', default=False, action='store_true')
parser.add_argument('--pad', default=1000, type=int)
parser.add_argument('--reference-name', default=None)
parser.add_argument('input')
args = parser.parse_args()

if args.subset is None:
    positions = None
else:
    positions = set()
    with open(args.subset, 'rb') as handle:
        for line in handle:
            contig, _, pos = line.strip().split()[:3]
            positions.add((contig, pos))

ref_name = args.reference_name
use_ref_sample = ref_name is not None

# Get offsets
pad = 0
contigs = []
pads = {}
seq = SeqIO.index(args.reference, 'fasta')
for c in sorted(seq):
    contigs.append((c, pad +1))
    pads[c] = pad
    pad += len(seq[c]) + args.pad

# Write contigs file
with open(args.output + '.contigs', 'wb') as out:
    for (c, s) in contigs:
        out.write(c + '\t' + str(s) + '\n')

with open(args.output, 'wb') as out:
    rdr = vcf.Reader(filename=args.input)
    if use_ref_sample:
        out.write('genome,Ref,' + ref_name + ',' + ','.join(rdr.samples) + '\n')
    else:
        out.write('genome,Ref,' + ','.join(rdr.samples) + '\n')
    for i, rec in enumerate(rdr):
        if (positions is not None) and ((rec.CHROM, str(rec.POS)) not in positions):
            continue
        if not rec.is_snp:
            continue
        c, p = rec.CHROM, rec.POS
        pad = pads[c]
        gts = []
        missed = 0
        for s in rec.samples:
            gt = re.split('\||/', s['GT'])
            if gt[0] == '.':
                g = 'N'
                missed += 1
            else:
                g = tuple(str(rec.alleles[int(x)]).upper() for x in gt)
                if len(g) == 1:
                    g = (g[0], g[0])
                if args.het_miss:
                    g = ambighm(g)
                else:
                    g = ambig(g)
            gts.append(g)
        if 1 - (missed / float(len(gts))) < args.min_gt:
            continue
        out.write(str(rec.POS + pad))
        if use_ref_sample:
            out.write(',' + rec.REF + ',' + rec.REF + ',' + ','.join(str(g) for g in gts) + '\n')
        else:
            out.write(',' + rec.REF + ',' + ','.join(str(g) for g in gts) + '\n')
