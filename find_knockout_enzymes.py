#!/usr/bin/env python2.7
"""
    Find restriction enzymes for inserting a cassette into a gene and
    inserting the cassette + gene into a vector.

    find_knockout_enzymes.py --output <output file>
        <gene sequence> <cassette sequence>
        <enzymes> <vector sequence, start, end>

    start and end of target region should be 1-based, inclusive.

    Assumes vector is circular and restriction enzymes are palindromes.

    Options:
        --pad (100)     Amount of the gene sequence that should be
                        off-limits to cutting by any chosen enzyme.

   TODO:
        - include some bases on either side of each sequence (except
          the cassette, full gene, and gene pad)
        - report the actual cut position (position just before cut or
          just after in rc) rather than the recognition sequence range
        - check for PstI and SbfI sites in MCS -- they should be there
"""

import argparse
import itertools
import re

from Bio import SeqIO, Seq


iupac = {'ACGT':'N',
         'AGT':'D', 'ACG':'V', 'ACT':'H', 'CGT':'B',
         'CT':'Y', 'AG':'R', 'AT':'W', 'CG':'S', 'GT':'K', 'AC':'M',
         'A':'A', 'C':'C', 'G':'G', 'T':'T'}
r_iupac = {v:k for k, v in iupac.iteritems()}

def iupac_match(s1, s2):
    bases1 = set(r_iupac[s1.upper()])
    bases2 = set(r_iupac[s2.upper()])
    if len(bases1.intersection(bases2)):
        return True
    else:
        return False

def find_cut_sites(seq, rs):
    seq = seq.upper()
    rs = rs.upper()
    rc_rs = Seq.reverse_complement(rs)
    matches = set()
    for i in xrange(len(seq) - len(rs) + 1):
        match = True
        rc_match = True
        for j in xrange(len(rs)):
            if not iupac_match(seq[i+j], rs[j]):
                match = False
                break
        for j in xrange(len(rc_rs)):
            if not iupac_match(seq[i+j], rc_rs[j]):
                rc_match = False
                break
        if match:
            matches.add((i+1, i+len(rs)))
        if rc_match:
            matches.add((i+1, i+len(rc_rs)))
    return matches


## In progress...
def find_cut_sites2(seq, rs, trailing):
    seq = seq.upper()
    cs = rs.index(' ') ## 0-based
    rc_cs = len(rs) - cs - 1
    rs = re.sub(' ', '', rs.upper())
    rc_rs = Seq.reverse_complement(rs)
    matches = set()
    rc_matches = set()
    for i in xrange(len(seq) - len(rs) + 1):
        match = True
        rc_match = True
        for j in xrange(len(rs)):
            if not iupac_match(seq[i+j], rs[j]):
                match = False
                break
        for j in xrange(len(rc_rs)):
            if not iupac_match(seq[i+j], rc_rs[j]):
                rc_match = False
                break
        if match:
            matches.add((i+1, i+len(rs)))
        if rc_match:
            rc_matches.add((i+1, i+len(rc_rs)))
    ## Go through the matches and pick out the actual cut site
    ## location. Only report if the cut site is in the target area.
    ## Then combine matches and rc_matches.
    cut_site_locs = set()
    start = trailing ## 1-based inclusive
    end = len(seq) - trailing + 1 ## 1-based inclusive
    for m1, m2 in itertools.izip(matches, rc_matches):
        s1, s2 = None, None
        if m1[0] + cs >= start and m1[0] + cs <= end:
            s1 = m1[0] + cs - trailing
        if m1[0] + cs == end + 1:
            s1 = end - trailing * 2
        if m2[0] + rc_cs >= start and m2[0] + rc_cs <= end:
            s2 = m2[0] + rc_cs - trailing
        if m2[0] + rc_cs == end + 1:
            s2 = end - trailing * 2
        if s1 is not None or s2 is not None:
            cut_site_locs.add((s1, s2))
    return cut_site_locs


parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--pad', type=int, default=100)
parser.add_argument('--output')
parser.add_argument('gene')
parser.add_argument('cassette')
parser.add_argument('enzymes')
parser.add_argument('vector', nargs=3)
args = parser.parse_args()

## Read sequences from fasta files
gene_seq = str(SeqIO.read(args.gene, 'fasta').seq)
cassette_seq = str(SeqIO.read(args.cassette, 'fasta').seq)
vector_seq = str(SeqIO.read(args.vector[0], 'fasta').seq)

## Split the vector sequence into the part that should be cut and the
## part that should not be cut. Note that we add region after
## the target region to the beginning of the region before the target --
## I think this accounts for the fact that the vector is circular and that
## the target region has been "chopped out" of the non-target.
target_vector_seq = vector_seq[(int(args.vector[1])-1):int(args.vector[2])]
non_target_vector_seq = vector_seq[int(args.vector[2]):] + \
        vector_seq[:(int(args.vector[1])-1)]
tvs_padded = vector_seq[(int(args.vector[1])-11):(int(args.vector[1])-1)] + \
        target_vector_seq + \
        vector_seq[int(args.vector[2]):(int(args.vector[2])+10)]
        
## Read enzyme cut sites
cut_sites = {}
with open(args.enzymes, 'rb') as ih:
    for line in ih:
        fields = line.strip().split('\t')
        if ' ' not in fields[1]:
            continue
        cut_sites[fields[0].strip()] = fields[1]

## Search for cut sites in gene
gene_enzymes = {}
for enzyme, rs in cut_sites.iteritems():
    match = find_cut_sites2(gene_seq, rs, 0)
    if len(match) > 0:
        gene_enzymes[enzyme] = match

## Search for cut sites in gene disregarding "padded" regions
gene_nopad_enzymes = {}
gene_nopad_seq = gene_seq[(args.pad-10):(-(args.pad-10))]
for enzyme, rs in cut_sites.iteritems():
    match = find_cut_sites2(gene_nopad_seq, rs, 10)
    if len(match) > 0:
        gene_nopad_enzymes[enzyme] = match

## Search for cut sites in "padded" regions
gene_pad_enzymes = {}
gene_pad_seq_1 = gene_seq[:args.pad]
gene_pad_seq_2 = gene_seq[(-(args.pad)):]
for enzyme, rs in cut_sites.iteritems():
    match = find_cut_sites2(gene_pad_seq_1, rs, 0)
    if len(match) > 0:
        gene_pad_enzymes[enzyme] = match
for enzyme, rs in cut_sites.iteritems():
    match = find_cut_sites2(gene_pad_seq_2, rs, 0)
    if len(match) > 0:
        gene_pad_enzymes[enzyme] = match

## Search for cut sites in cassette
cassette_enzymes = {}
for enzyme, rs in cut_sites.iteritems():
    match = find_cut_sites2(cassette_seq, rs, 0)
    if len(match) > 0:
        cassette_enzymes[enzyme] = match

## Search for cut sites in target region of vector
target_vector_enzymes = {}
for enzyme, rs in cut_sites.iteritems():
    match = find_cut_sites2(tvs_padded, rs, 10)
    if len(match) > 0:
        target_vector_enzymes[enzyme] = match

## Search for cut sites in non-target region of vector
non_target_vector_enzymes = {}
for enzyme, rs in cut_sites.iteritems():
    match = find_cut_sites2(non_target_vector_seq, rs, 0)
    if len(match) > 0:
        non_target_vector_enzymes[enzyme] = match

## Choose candidates
gene_cutter_candidates = set(e for e, s in gene_nopad_enzymes.iteritems() if 0 < len(s) < 3) - \
        set(gene_pad_enzymes) - \
        set(cassette_enzymes)
vector_cutter_candidates = set(target_vector_enzymes) - \
        set(non_target_vector_enzymes) - \
        set(gene_enzymes) - \
        set(cassette_enzymes)

print 'Step 1: %i candidates' % len(gene_cutter_candidates)
print 'Step 2: %i candidates' % len(vector_cutter_candidates)

## This should be added:
## Among candidates, make gene + cassette constructs and
## search for cut sites

## Write output
with open(args.output, 'wb') as oh:
    oh.write('Gene candidates\n')
    oh.write('Enzyme\tcut site\tnumber of gene cuts\n')
    for e in gene_cutter_candidates:
        oh.write(e + '\t' + cut_sites[e] + '\t' +
                 str(len(gene_nopad_enzymes[e])) + '\n')
    oh.write('\nPlasmid candidates\n')
    oh.write('Enzyme\tcut site\tnumber of plasmid cuts\n')
    for e in vector_cutter_candidates:
        oh.write(e + '\t' + cut_sites[e] + '\t' +
                 str(len(target_vector_enzymes[e])) + '\n')
