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
import collections
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
        if replicon in self.positions:
            p = self.positions[replicon]
        else:
            return []
        l = bisect.bisect_left(p, ref_start)
        if l >= len(p):
            return []
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
sys.stderr.write('Read VCF file\n')

b = pysam.Samfile(infile)

recombination_count = collections.defaultdict(list)

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
    read_haps = haps.get_haplotypes(replicon, start+1, end)

    snp_pos_in_read = read_haps[0]
    haps_in_read = read_haps[1]
    if len(snp_pos_in_read) > 1:
        for i in xrange(len(snp_pos_in_read)-1):

            # Target SNP and all to the left
            p0 = snp_pos_in_read[:(i+1)]
            p0_idx = []
            for j in p0:
                try:
                    p0_idx.append(ref_positions.index(j-1))
                except ValueError:
                    # This can happen if there are indels
                    pass
            if len(p0_idx) == 0:
                continue
            read_hap0 = [seqn[j] for j in p0_idx]

            # All to the right of target SNP
            p1 = snp_pos_in_read[(i+1):]
            p1_idx = []
            for j in p1:
                try:
                    p1_idx.append(ref_positions.index(j-1))
                except ValueError:
                    pass
            if len(p1_idx) == 0:
                continue
            read_hap1 = [seqn[j] for j in p1_idx]

            breakpoint = (p0[-1] + p1[0]) / 2.

            probs0 = []
            probs1 = []
            for h in haps_in_read:
                # Probability that read_hap0 is from this strain
                # Simple match / no match method:
                for j in xrange(len(read_hap0)):
                    pr = 1
                    err = 10**(-qual[p0_idx[j]])/10.
                    if h[j] is None:
                        pr *= 0.25
                    elif h[j] == read_hap0[j]:
                        pr *= 1 - err
                    else:
                        pr *= err / 3.
                probs0.append(pr)

                # Probability that read_hap1 is from this strain
                for j in xrange(len(read_hap1)):
                    pr = 1
                    err = 10**(-qual[p1_idx[j]])/10.
                    if h[j] is None:
                        pr *= 0.25
                    elif h[j] == read_hap1[j]:
                        pr *= 1 - err
                    else:
                        pr *= err / 3.
                probs1.append(pr)

            max_same = 0.
            max_diff = 0.
            for j in xrange(len(probs0)):
                for k in xrange(len(probs1)):
                    pr = probs0[j] * probs1[k]
                    if j == k and pr > max_same:
                        max_same = pr
                    if j != k and pr > max_diff:
                        max_diff = pr

            recombination_count[breakpoint].append(max_diff > max_same)

            # NOTE: Adding read quality didn't seem to change the results,
            # although it did make the program run a bit longer.

    else:
        # Not enough SNPs in this read to check anything
        pass

print 'pos\trecombine\tno_recombine'
for bp in sorted(recombination_count):
    rc = recombination_count[bp]
    print str(bp) + '\t' + str(sum(1 for x in rc if x)) + '\t' + \
            str(sum(1 for x in rc if not x))

