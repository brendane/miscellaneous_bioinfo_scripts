#!/usr/bin/env python2.7
"""
    Extract SNPs from a reference genome + the results of show-snps +
    the results of show-coords.

    mummer_extract_snps.py coords snps

    Requires biopython.

    Does not pay attention to copy-number or duplication status.
"""

import argparse
import csv
import os.path as osp

from Bio import SeqIO

iupac = {'ACGT':'N',
         'AGT':'D', 'ACG':'V', 'ACT':'H', 'CGT':'B',
         'CT':'Y', 'AG':'R', 'AT':'W', 'CG':'S', 'GT':'K', 'AC':'M',
         'A':'A', 'C':'C', 'G':'G', 'T':'T'}
r_iupac = {v:k for k, v in iupac.iteritems()}

def iupac_combine(a, b):
    aa = r_iupac[a]
    bb = r_iupac[b]
    return iupac[''.join(sorted(set(aa + bb)))]


parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('coords')
parser.add_argument('snps')
parser.add_argument('name')
args = parser.parse_args()

# Read the coords file
ref_file = None
qry_file = None
ref_covered = {}
with open(args.coords, 'rb') as ih:
    ref_file, qry_file = ih.readline().strip().split()

    ih.readline(); ih.readline();
    header = ih.readline().strip()
    if header != '[S1]\t[E1]\t[S2]\t[E2]\t[LEN 1]\t[LEN 2]\t[TAGS]':
        raise Exception('Columns are not as expected in coords file')

    for line in ih:
        row = line.strip().split('\t')
        s, e, c = int(row[0]), int(row[1]), row[6]
        if c not in ref_covered:
            ref_covered[c] = []
        if s > e:
            ss = e
            e = s
            s = ss
        ref_covered[c].append((s, e))

# Get the part of the reference sequence that is covered
ref_idx = SeqIO.index(ref_file, 'fasta')
filled_in_seq = {}
for c in ref_idx:
    seq = [x.upper() for x in ref_idx[c].seq]
    for b in range(len(seq)):
        if c not in ref_covered:
            seq[b] = 'N'
        else:
            found = False
            for (s, e) in ref_covered[c]:
                if b >= s-1 and b <= e-1:
                    found = True
                    break
            if not found:
                seq[b] = 'N'
    filled_in_seq[c] = seq

# Modify with SNPs
with open(args.snps, 'rb') as ih:
    rf, qf = ih.readline().strip().split()
    if rf != ref_file or qf != qry_file:
        raise Exception('coords and snps files do not match')

    ih.readline(); ih.readline();
    header = ih.readline().strip()
    if header != '[P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[R]\t[Q]\t[FRM]\t[TAGS]':
        raise Exception('Columns are not as expected in snps file')

    prev = ('.', -1)
    for line in ih:
        row = line.strip().split('\t')
        p, c = row[0], row[10]
        r, q = row[1], row[2]
        if r == '.':
            continue
        if q == '.':
            q = 'N'
        rs = (c, p)
        dup = False
        if rs == prev:
            dup = True
        if not dup and r != filled_in_seq[c][int(p)-1]:
            raise Exception('Ref does not match at %s %s' % (c, p))
        prev = rs
        q = iupac_combine(q, filled_in_seq[c][int(p)-1])
        filled_in_seq[c][int(p)-1] = q

with open(args.output, 'wb') as oh:
    oh.write('##fileformat=VCFv4.1\n')
    oh.write('##source=mummer_extract_snps\n')
    oh.write('##reference=' + osp.basename(ref_file) + '\n')
    oh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    oh.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' +
             args.name + '\n')
    for c in sorted(filled_in_seq):
        r = str(ref_idx[c].seq)
        for i, b in enumerate(filled_in_seq[c]):
            # TODO: need to get the ref/alt for multiple alleles
            # and no variation working
            gt = str(1-int(b == r[i]))
            if b == 'N':
                gt = '.'
            if gt == '0':
                alt = '.'
            elif gt == '1':
                alt = b
            oh.write(c + '\t' + str(i+1) + '\t' + c + '-' + str(i+1) +
                     '\t' + r[i] + '\t' + b + '\t' + '0\t.\t.\tGT\t' +
                     str(1-int(b == r[i])) + '\n')