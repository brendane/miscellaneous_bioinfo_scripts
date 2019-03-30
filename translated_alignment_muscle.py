#!/usr/bin/env python3
"""
    Align protein coding nucleotide sequences based on the alignments of
    their translations.

    translated_alignment_muscle.py <fasta input>
"""

import argparse
from io import StringIO
import re
import subprocess
import sys

from Bio import AlignIO, SeqIO
from Bio.Data.CodonTable import TranslationError

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--table', default='Bacterial')
parser.add_argument('input')
args = parser.parse_args()

nucs = {}
prots = []
for rec in SeqIO.parse(args.input, 'fasta'):
    ## Read input file and translate
    nucs[rec.description] = rec
    try:
        p = str(rec.seq.translate(table=args.table, cds=True))
    except TranslationError:
        p = str(rec.seq.translate(table=args.table))
    prots.append((rec.description, p))

## Align
muscle = subprocess.Popen(args=['muscle'],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE)
prot_input = ['>' + rec[0] + '\n' + re.sub('\\*', 'X', rec[1]) + '\n' for rec in prots]
prot_aln_string = muscle.communicate(bytearray('\n'.join(prot_input)+'\n', 'utf-8'))[0].decode('utf-8')
if prot_aln_string == '':
    raise Exception('Empty protein alignment')
prot_aln = AlignIO.read(StringIO(prot_aln_string), 'fasta') 

## Use protein alignment to align nucleotides
for p in prot_aln:
    n = nucs[p.description]
    aln_n = []
    codon = 0
    for aa in p:
        if aa == '-':
            aln_n.append('---')
        else:
            s = codon*3
            e = codon*3+3
            if e > len(n): e = len(n)
            aln_n.append(str(n[s:e].seq))
            codon += 1
    sys.stdout.writelines(['>', p.description, '\n', ''.join(aln_n), '\n'])
