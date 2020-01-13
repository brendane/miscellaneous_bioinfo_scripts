#!/usr/bin/env python3
"""
    Convert a tab-delimited table of LCB copy-number into a VCF file with
    presence and absence.

    lcbtsv2vcf.py [--keep <strains to keep>] <input tsv file>
"""

import argparse
import csv
import sys

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--keep')
parser.add_argument('tsv')
args = parser.parse_args()

keep = []
if args.keep is not None:
    with open(args.keep, 'rt') as ih:
        for line in ih:
            keep.append(line.strip())


with open(args.tsv, 'rt') as ih:
    rdr = csv.DictReader(ih, delimiter='\t')
    strains = rdr.fieldnames[6:]
    if args.keep is None:
        keep = strains
    sys.stdout.write('##fileFormat=VCFv4.2\n')
    sys.stdout.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' +
                     '\t'.join(keep) + '\n')
    for row in rdr:
        sys.stdout.write('lcb\t' + row['lcb'] + '\tlcb-' + row['lcb'] + '\t' +
                         'A\tP\t.\t.\t.\tGT')
        for strain in keep:
            if int(row[strain]) > 0:
                sys.stdout.write('\t1')
            else:
                sys.stdout.write('\t0')
        sys.stdout.write('\n')
