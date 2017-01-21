#!/usr/bin/env python2
"""
    Produce some statistics from the output of plink --freq --within.

    strat_freq_counter.py --output <output file> <input file>
"""

#==============================================================================#

import argparse
import itertools

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('input')
args = parser.parse_args()

n_unique = 0
n_unique_fixed = 0
n_fixed = 0

with open(args.output, 'wb') as out:
    with open(args.input, 'rb') as handle:
        handle.readline()
        for snp, lines in itertools.groupby(handle, lambda x: x.strip().split()[1]):
            lines = [l.strip().split() for l in lines]
            # MAFs from groups that are represented by at least two
            # individuals
            mafs = [float(l[5]) for l in lines if int(l[7]) / 2 > 1]
            n_pops = len(mafs)

            if n_pops < 2:
                continue

            snp_type = 'none'

            maf_0 = sum(1 for m in mafs if m == 0)
            maf_1 = sum(1 for m in mafs if m == 1)

            # Allele fixed in a population
            #   MAF == 1 or == 0 in >= 1 population
            if (maf_0 + maf_1) > 0:
                snp_type = 'fixed'
                n_fixed += 1

            # Allele unique to one population
            #   MAF == 1 in all but one population or MAF == 0 in all but
            #   one population
            # Note that this comes after "fixed" b/c if an allele is
            # unique to one population, then the other allele must be
            # fixed in the other populations.
            if n_pops - maf_0 == 1 or n_pops - maf_1 == 0:
                snp_type = 'unique'
                n_unique += 1

            # Allele unique to one population and fixed in that population
            #   As above, but MAF == 0 or == 1 in the remaining population
            if (maf_0 + maf_1) == n_pops and (maf_1 == 1 or maf_0 == 1):
                snp_type = 'unique_fixed'
                n_unique_fixed += 1

            out.write(lines[0][0] + '\t' + lines[0][1] + '\t' + snp_type +
                      '\n')

print 'unique\tunique_fixed\tfixed'
print str(n_unique)  + '\t' + str(n_unique_fixed) + '\t' + \
        str(n_fixed)
