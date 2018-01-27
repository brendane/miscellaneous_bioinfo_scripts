#!/usr/bin/env python2.7
"""
    Combine the output from multiple runs of mummer_extract_snps into
    a single file that is ready for use with HARP.

    combine_mummer_extracted_snps.py --output <output prefix>
        <file with list of input files>

    The input files should have replicons and positions in the same
    order and should have a genotype call for every position.

    The output file has several characteristics:
        - DGRP format for HARP
        - Ambiguous genotype calls given ambiguity codes
        - Only segregating sites
        - All contigs are combined into a single molecule with 1000
          bp in between; a second output file gives the order and start
          and end points of each real contig in the combined "contig".
"""

import argparse
import itertools

iupac = {'ACGT':'N',
         'AGT':'D', 'ACG':'V', 'ACT':'H', 'CGT':'B',
         'CT':'Y', 'AG':'R', 'AT':'W', 'CG':'S', 'GT':'K', 'AC':'M',
         'A':'A', 'C':'C', 'G':'G', 'T':'T'}
r_iupac = {v:k for k, v in iupac.iteritems()}

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('inputlist')
args = parser.parse_args()

# Open all the file handles
with open(args.inputlist, 'rb') as ih:
    file_handles = [open(f.strip(), 'rb') for f in ih]

# Get strain ids
strains = [h.readline().strip().split('\t')[-1] for h in file_handles]
print 'Found %i strains' % len(strains)

# Iterate over files
pad = 0
contigs = [] # (name, start)
gts = [] # [pos, ref, gt calls...]
prev_contig = None
prev_pos = 0
for lines in itertools.izip(*file_handles):
    gt = []
    for i, l in enumerate(lines):
        fields = l.strip().split('\t')
        if i == 0:
            c, p = fields[:2]
            if prev_contig is None:
                contigs.append((c, pad + 1))
            if prev_contig is not None and c != prev_contig:
                pad += prev_pos + 1000
                contigs.append((c, pad + 1))
            gt.append(int(p) + pad)
            gt.append(fields[2].upper())
            prev_contig = c
            prev_pos = int(p)
        else:
            if (c, p) != tuple(fields[:2]):
                raise Exception('Files do not match at file %i' %i)
        gt.append(fields[3])
    # Check if this site is segregating
    base_calls = set(r_iupac[g] for g in gt[2:] if g != 'N')
    if len(base_calls) > 1:
        gts.append(gt)


# Check that there aren't any lines left in files
# (rough check)
for i, h in enumerate(file_handles):
    if h.readline() != '':
        raise Exception('File %i has lines left' %i)

# Write the output
with open(args.output + '.snps.txt', 'wb') as oh:
    oh.write('genome,Ref,' + ','.join(strains) + '\n')
    for gt in gts:
        oh.write(','.join(str(x) for x in gt) + '\n')

# Write contig start positions in combined contig
with open(args.output + '.contigs.txt', 'wb') as oh:
    for (c, s) in contigs:
        oh.write(c + '\t' + str(s) + '\n')

# Close file handles
for h in file_handles:
    h.close()
