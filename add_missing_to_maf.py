#!/usr/bin/env python3
"""
    Add missing pieces to an MAF file. Especially meant for processing
    mugsy output files, which often do not quite cover all the input
    genomes. Expects mugsy naming conventions.

    add_missing_to_maf.py <fasta directory> [<input maf file>]

    Assumes no whitespace in the sequences.
"""

import argparse
import collections
import os
import sys
import subprocess
import tempfile

from Bio import SeqIO

def parse_maf(ih, print_line=False, handle=None):
    record = {'a':{}, 's':{}}
    for line in ih:
        if print_line:
            handle.write(line)
        if line.startswith('a'):
            if len(record['s']) > 0:
                yield record
            record = {'a':{}, 's':{}}
            record['a'] = {x.split('=')[0]:x.split('=')[1] for x in line.strip().split()[1:]}
            record['s'] = {}
        elif line.startswith('s'):
            fields = line.strip().split()[1:]
            record['s'][(fields[0], int(fields[1]), int(fields[2]), fields[3], int(fields[4]))] = fields[5]
        else:
            continue
    yield record


parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('fasta')
parser.add_argument('maf', nargs='*')
args = parser.parse_args()

lcb_bed = collections.defaultdict(list)

if args.maf is None:
    handle = sys.stdin
else:
    handle = open(args.maf[0])

with handle as ih:
    ## Loop through LCBs and record coordinates in each strain
    for rec in parse_maf(ih, True, sys.stdout):
        names = list(rec['s'].keys())
        for n in names:
            strain, replicon = n[0].split('.', 1)
            if n[3] == '+':
                ss = n[1]
                ee = n[1] + n[2]
            else:
                ee = n[4] - n[1]
                ss = n[4] - n[1] - n[2]
            lcb_bed[strain].append((replicon, ss, ee))

for strain, lcb_coords in lcb_bed.items():
    ## For each strain, used bedtools to find the missing pieces of the
    ## genome. Then write each piece as an LCB.
    lcb_coords.sort()
    tmp = tempfile.mkstemp(dir='.')[1]
    tmp_genome = tempfile.mkstemp(dir='.')[1]
    seq_idx = SeqIO.index(args.fasta + '/' + strain + '.fasta', 'fasta')

    ## Make a bed file with positions that are included
    with open(tmp, 'wt') as oh:
        for c, s, e in lcb_coords:
            oh.write(c + '\t' + str(s) + '\t' + str(e) + '\n')

    ## Make a "genome file" with contig/chromosome lengths
    lens = {}
    with open(tmp_genome, 'wt') as oh:
        for i in sorted(seq_idx.keys()):
            l = len(seq_idx[i])
            lens[i] = l
            oh.write(i + '\t' + str(l) + '\n')

    ## Run bedtools
    p = subprocess.Popen(' '.join(['bedtools', 'sort', '-i', tmp, '|',
                                   'bedtools', 'merge', '-i', '-', '|',
                                   'bedtools', 'complement', '-i', '-',
                                   '-g', tmp_genome]),
                         shell=True, stdout=subprocess.PIPE)
    x = p.communicate()[0].decode('utf-8').split('\n')
    os.unlink(tmp); os.unlink(tmp_genome)

    for i, rec in enumerate(x):
        if len(rec) == 0: continue
        chrom, s, e = rec.strip().split('\t')
        s = int(s); e = int(e)
        seqrec = seq_idx[chrom][s:e]
        sys.stdout.write('a label=%s_%i\n' % (strain, i))
        sys.stdout.write('\t'.join(['s', strain + '.' + chrom, str(s),
                                    str(e-s), '+', str(lens[chrom]),
                                    str(seqrec.seq)]) + '\n\n')
