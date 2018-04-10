#!/usr/bin/env python2.7
"""
    Find restriction enzymes in records from a FASTA file and report the
    number of matches and their positions.

    find_re_sites.py --output <output file>
        <enzymes file> <fasta file>

    Assumes linear sequences.
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
parser.add_argument('--output')
parser.add_argument('enzymes')
parser.add_argument('fasta')
args = parser.parse_args()


## Read enzyme cut sites
cut_sites = {}
with open(args.enzymes, 'rb') as ih:
    for line in ih:
        fields = line.strip().split('\t')
        if ' ' not in fields[1]:
            continue
        cut_sites[fields[0].strip()] = fields[1]

## Search for cut sites in sequences
re_sites = {}
for rec in SeqIO.parse(args.fasta, 'fasta'):
    re_sites[rec.id] = {}
    seq = str(rec.seq)
    for enzyme, rs in cut_sites.iteritems():
        match = find_cut_sites2(seq, rs, 0)
        if len(match) > 0:
            re_sites[rec.id][enzyme] = match

## Write number of sites
with open(args.output + '.nsites.txt', 'wb') as oh:
    oh.write('enzyme\t' + '\t'.join(re_sites.keys()) + '\n')
    for enzyme in cut_sites.iterkeys():
        oh.write(enzyme)
        for seqid in re_sites.keys():
            if enzyme not in re_sites[seqid]:
                count = 0
            else:
                count = len(re_sites[seqid][enzyme])
            oh.write('\t' + str(count))
        oh.write('\n')

## Write site locations
with open(args.output + '.locations.txt', 'wb') as oh:
    oh.write('enzyme\t' + '\t'.join(re_sites.keys()) + '\n')
    for enzyme in cut_sites.iterkeys():
        oh.write(enzyme)
        for seqid in re_sites.keys():
            sites = []
            if enzyme in re_sites[seqid]:
                for site in re_sites[seqid][enzyme]:
                    sites.append(','.join(str(s) for s in site))
            oh.write('\t' + ';'.join(sites))
        oh.write('\n')
