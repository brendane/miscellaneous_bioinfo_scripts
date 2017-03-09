#!/usr/bin/env python
"""
    Assess how likely it is that two known variant sites in a read came
    from different original haplotypes.

    Will probably have to be re-written in a faster language.

    Assumes haploid genomes.

    VCF file must be sorted - and this script does not check.
"""

#==============================================================================#

import bisect
import itertools
import sys

import pysam

#==============================================================================#

def mlt(x, y): return x*y

class Haplotypes:
    """ Class to keep track of haplotype data. Currently does not hold
    any information on genotype likelihoods. Positions are 1-based."""

    def __init__(self):
        self.positions = {}
        self.haplotypes = {}
        self.replicons = []
        self.samples = []

    def populate(self, vf):
        for replicon, records in itertools.groupby(vf, lambda x: x.chrom):
            self.samples = vf.header.samples
            self.replicons.append(replicon)
            self.haplotypes[replicon] = {s:[] for s in self.samples}
            h = self.haplotypes[replicon]
            self.positions[replicon] = []
            p = self.positions[replicon]
            for rec in records:
                p.append(rec.pos)
                for s in self.samples:
                    # Take just the first allele b/c we're assuming haploid
                    h[s].append(rec.samples[s].alleles[0])

    def variant_pos(self, replicon, ref_start, ref_end):
        """ Adapted from the bisect documentation. Get variant indices
        given reference coordinates. """
        p = self.positions[replicon]
        l = bisect.bisect_left(p, ref_start)
        if p[l] > ref_end:
            return []
        r = bisect.bisect_right(p, ref_end)
        if p[r] < ref_start:
            return []
        return range(l, r)

    def get_haplotypes(self, replicon, ref_start, ref_end):
        """ Get the SNP haplotypes (really strain base calls) for a range
        of positions. Missing data is encoded as None. Return value is
        a list of lists. Each sublist contains the base calls for a
        strain ordered by position."""
        pos_idx = self.variant_pos(replicon, ref_start, ref_end)
        if len(pos_idx) == 0:
            return([], [])
        h = self.haplotypes[replicon]
        haps = []
        pos  = [self.positions[replicon][i] for i in pos_idx]
        for s in self.samples:
            haps.append([h[s][i] for i in pos_idx])
        return (pos, haps)

#==============================================================================#

infile = sys.argv[1]
vcffile = sys.argv[2]

haps = Haplotypes()
v = pysam.VariantFile(vcffile)
haps.populate(v)

import pdb; pdb.set_trace()

b = pysam.Samfile(infile, 'rb')

for read in b:

    if read.is_unmapped:
        continue

    # Get read information
    # Positions are 0-based
    replicon = read.reference_name
    ref_positions = read.get_reference_positions()
    start = read.reference_start
    end = read.reference_end
    qual = read.query_qualities
    seqn = read.seq

    # Get the haplotypes for this read (just SNPs)
    read_haps = haps.get_haplotypes(replicon, start, end)

    # Iterate through adjacent pairs of SNPs:
    #   1. Determine existing haplotypes
    #   2. Calculate probability of being from the same strain
    #      and probability of being from different strains
    #   3. Calculate and record some summary statistic that says
    #      how P(different) compares to P(same)
    #   4. Record based on the position halfway between SNP1 and SNP2
    #
    # Really do this using all the SNPs to the left and right of the
    # split.
    #
    # Might want to filter on base quality if it isn't incorporated
    # into calculations
    #
    # Watch out for 1-based and 0-based issues; also figure out what
    # an unaligned part of the read looks like.

    #
    # snp locations in the read = x
    #
    # Read:    ---------x----------x-------------------x----
    #                   1          2                   3
    #               
    # First iteration: Test 1 vs (2,3)
    #
    # Second iteration: Test (1,2) vs 3
    #

    snp_pos_in_read = read_haps[0]
    haps_in_read = read_haps[1]
    if len(snp_pos_in_read) > 1:
        for i in xrange(len(snp_pos)-1):

            # Target SNP and all to the left
            p0 = snp_pos_in_read[:(i+1)]
            p0_idx = [ref_positions.index(j) for j in p0]
            read_hap0 = [seqn[j] for j in p0_idx]

            # All to the right of target SNP
            p1 = snp_pos_in_read[(i+1):]
            p1_idx = [ref_positions.index(j) for j in p1]
            read_hap1 = [seqn[j] for j in p1_idx]

            # I don't want 0 and 1
            max_match = 0.999
            min_match = 0.001
            probs0 = []
            probs1 = []
            for h in haps_in_read:
                # Probability that read_hap0 is from this strain
                # Simple match / no match method:
                for j in xrange(read_hap0):
                    pr = 1
                    if h[j] is None:
                        pr *= (1 - 0.25)
                    else if h[j] == read_hap0[j]:
                        pr *= (1 - max_match)
                    else:
                        pr *= (1 - min_match)
                probs0.append(pr)

                # Probability that read_hap1 is from this strain
                for j in xrange(read_hap1):
                    pr = 1
                    if h[j] is None:
                        pr *= (1 - 0.25)
                    else if h[j] == read_hap1[j]:
                        pr *= (1 - max_match)
                    else:
                        pr *= (1 - min_match)
                probs1.append(pr)

            # Now calculate overall probability of being from the same
            # strain and overall probability of being from different
            # strains
            p_same = 1
            for pr0, pr1 in itertools.izip(probs0, probs1):
                p_same *= (1 - (pr0 * pr1))
            p_same = 1 - p_same
            p_diff = 1
            for j in xrange(len(pr0)):
                for k in xrange(len(pr1)):
                    if j == k:
                        continue
                    p_diff *= (1 - (pr0 * pr1))
            p_diff = 1 - p_diff

    else:
        # Not enough SNPs in this read to check anything
        pass

