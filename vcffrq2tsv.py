#!/usr/bin/env python2
"""
    Convert the .frq file produced by vcftools --freq to a more parseable
    tsv format.

    vcffrq2tsv.py --output <output file> <input file>
"""

#==============================================================================#

import argparse
import csv
import re

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('input')
args = parser.parse_args()

with open(args.output, 'wb') as out:
    out.write('chrom\tpos\tN\tref\tref_frq\talt_frq\talleles\tfreqs\n')
    with open(args.input, 'rb') as handle:
        rdr = csv.reader(handle, delimiter='\t')
        rdr.next()
        for row in rdr:
            cont = row[0]
            pos = row[1]
            N = str(int(row[3])/2)
            alleles = [re.sub(':.+', '', x) for x in row[4:]]
            frqs = [re.sub('.+:', '', x) for x in row[4:]]
            ref_allele = alleles[0]
            out.write(cont + '\t' + pos + '\t'+ N + '\t' + ref_allele + '\t' +
                      frqs[0] + '\t' + str(1 - float(frqs[0])) + '\t' +
                      ','.join(alleles) + '\t' + ','.join(frqs) + '\n')
