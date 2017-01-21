#!/usr/bin/env python2

#==============================================================================#

import argparse
import collections
import csv
import math
import re

#==============================================================================#

def calc_stats(data, snps):
    # pi, thetaW, D, thetaL, thetaH, nseg
    # All calculated using the formula from libsequence
    nseg = len(snps)
    if nseg > 0:

        # Pi and theta
        pi = 0.0
        theta_w = 0.0
        theta_l = 0.0
        theta_l_snps = 0
        max_t = 0
        for s in snps:
            d = data[s]
            al1, al2 = d[2], d[3]
            anc = d[7]
            n1, n2 = d[4], d[5]
            t = float(d[6])
            if t > max_t:
                max_t = t
            x = (n1 * (n1 - 1)) / (t * (t - 1))
            x += (n2 * (n2 - 1)) / (t * (t - 1))
            pi += 1 - x
            theta_w += 1. / sum(1 / float(i) for i in range(1, int(t)))
            if anc != 'N' and al1 == anc:
                nderiv = n2
                nanc = n1
            elif anc != 'N' and al2 == anc:
                nderiv = n1
                nanc = n2
            else:
                nderiv = None
                nanc = None
            if nderiv is not None:
                theta_l_snps += 1
                theta_l += nderiv / float(t - 1)


        # Tajima's D
        a1 = sum(1.0/i for i in range(1, int(max_t)))
        a2 = sum(1.0/(i ** 2) for i in range(1, int(max_t)))
        b1 = (max_t + 1.0) / (3.0 * (max_t - 1.0))
        b2 = (2.0 * (max_t ** 2 + max_t + 3.0)) / (9.0 * max_t * (max_t - 1.0))
        c1 = b1 - 1.0 / a1
        c2 = b2 - (max_t + 2.0) / (a1 * max_t) + a2 / (a1 ** 2)
        e1 = c1 / a1
        e2 = c2 / (a1 ** 2 + a2)
        d_denom = math.sqrt(e1 * nseg + e2 * nseg * (nseg - 1))
        D = (pi - theta_w) / d_denom

        # Fay and Wu's H
        if theta_l_snps > 0:
            b11 = sum(1.0/(i ** 2) for i in range(1, int(max_t) + 1))
            thetasq = nseg * (nseg - 1) / (a1 * a1 + b1)
            var_theta_l = (max_t * theta_w) / (2.0 * (max_t - 1)) + \
                    (2.0 * (max_t / (max_t - 1.0))**2 * (b11 - 1.0) - 1.0) * \
                    thetasq
            var_pi = (3.0 * max_t * (max_t + 1.0) * theta_w + \
                    2.0 * (max_t * max_t + max_t + 3.0) * thetasq) / \
                    (9 * max_t * (max_t - 1.0))
            cov = ((max_t + 1.0) / (3.0 * (max_t - 1.0))) * theta_w + \
                    ((7.0 * max_t * max_t + 3.0 * max_t - 2.0 - \
                      4.0 * max_t * (max_t + 1.0) * b11) / \
                   (2.0 * pow((max_t - 1.0), 2.0))) * thetasq
            h_denom = math.sqrt(var_theta_l + var_pi - 2.0 * cov)
            H = (pi - theta_l) / h_denom
        else:
            H = float('nan')

    else:
        pi = 0.0
        theta_w = 0.0
        theta_l = 0.0
        H = float('nan')
        D = float('nan')

    return {'nseg':nseg, 'pi':pi, 'theta_w':theta_w, 'theta_l':theta_l,
            'D':D, 'H':H}

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('outgroup')
parser.add_argument('genes')
parser.add_argument('vartypes')
parser.add_argument('frq')
args = parser.parse_args()

##
## Read all the data files
##

# Read outgroup states
outgroup = {} # key = (chr, pos), value = allele
with open(args.outgroup, 'rb') as handle:
    rdr = csv.reader(handle, delimiter='\t')
    for row in rdr:
        outgroup[(row[0], int(row[1]))] = row[2].upper()

# Genic bed file
genes = []
with open(args.genes, 'rb') as handle:
    rdr = csv.reader(handle, delimiter='\t')
    for row in rdr:
        genes.append((row[0], int(row[1]) + 1,  int(row[2]), row[3]))

# Variant type file
var_types = {}
with open(args.vartypes, 'rb') as handle:
    rdr = csv.reader(handle, delimiter='\t')
    for row in rdr:
        var_types[(row[0], int(row[2]))] = row[3]

# frq.count from vcftools --counts; should be sorted
data = []
gene_vars = collections.defaultdict(set)
with open(args.frq, 'rb') as handle:
    handle.readline()
    rdr = csv.reader(handle, delimiter='\t')
    i = 0
    for row in rdr:
        if len(row) != 6:
            # Not bi-allelic
            continue
        chrom = row[0]
        pos = int(row[1])
        allele1 = re.sub(':.+', '', row[4])
        allele2 = re.sub(':.+', '', row[5])
        n1 = int(re.sub('.:', '', row[4])) / 2
        n2 = int(re.sub('.:', '', row[5])) / 2
        try:
            out_allele = outgroup[(chrom, pos)]
        except KeyError:
            out_allele = 'N'
        total = int(row[3]) / 2
        var_type = var_types[(chrom, pos)]
        for gene in genes:
            if gene[0] == chrom and (gene[1] <= pos <= gene[2]):
                gene_vars[gene[3]].add(i)
        data.append((chrom, pos, allele1, allele2, n1, n2,
                     total, out_allele, var_type))
        i += 1

##
## Calculate statistics for each gene
##
for i, gene in enumerate(genes):
    gene_name = gene[3]
    snps = gene_vars[gene_name]
    stats = calc_stats(data, snps)
    if i == 0:
        print 'gene\tnseg\tpi\ttheta_w\ttheta_l\tD\tH'
    print gene_name + '\t' + str(stats['nseg']) + '\t' + \
            str(stats['pi']) + '\t' + str(stats['theta_w']) + '\t' + \
            str(stats['theta_l']) + '\t' + str(stats['D']) + '\t' + \
            str(stats['H'])

##
## Calculate statistics for the whole genome
##
stats = calc_stats(data, range(len(data)))
