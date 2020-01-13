#!/usr/bin/env python3
"""
    Merge multi-stage whole genome alignments

    maf_merge_from_consensus.py --output <output file> <input maf>
        <input_name:input_file...>
"""

import argparse
import collections

from Bio import Seq

MafSeq = collections.namedtuple('MafSeq', 'seq_name start aligned_bases strand contig_length seq')

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
            record['s'] = []
        elif line.startswith('s'):
            fields = line.strip().split()[1:]
            record['s'].append(MafSeq(fields[0], int(fields[1]), int(fields[2]), fields[3], int(fields[4]), fields[5]))
        else:
            continue
    yield record

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('maf')
parser.add_argument('inputs', nargs='+')
args = parser.parse_args()

## Goal:
## MAF file with coordinates of individual genomes

## 1. Map between consensus coordinates and individual genome coordinates
##    args.inputs = batch_0:results/.../round1-batch.0/processed.maf
## Do this for each seq_name separately to parallelize
full_coords = {}
contig_lengths = {}
individual_contigs = collections.defaultdict(list)
for seq_name, maf_file in map(lambda x: x.split(':'), args.inputs):
    coords = {} # consensus position: genome_position, direction, base
    with open(maf_file, 'rt') as handle:
        for rec in parse_maf(handle):
            label = rec['a']['label']
            length = len(rec['s'][0].seq)
            print(seq_name, label)
            for i in range(length):
                # i = consensus position
                coords[(seq_name + '.' + label, i)] = {}
            for srec in rec['s']:
                direction = 1
                j = srec.start
                if srec.strand == '-':
                    direction = -1
                    j = srec.contig_length - srec.start - 1
                for i in range(length):
                    if i == 0:
                        contig_lengths[srec.seq_name] = srec.contig_length
                        individual_contigs[seq_name + '.' + label].append(srec.seq_name)
                    if srec.seq[i] == '-':
                        coords[(seq_name + '.' + label, i)][srec.seq_name] = (None, direction, srec.seq[i])
                    else:
                        coords[(seq_name + '.' + label, i)][srec.seq_name] = (j, direction, srec.seq[i])
                        j += direction
    full_coords[seq_name] = coords


## 2. Walk through alignment of consensuses, substituting in the
##    coordinates and bases from the individual genomes
## Do this for ind_seq_name separately to parallelize
with open(args.output, 'wt') as oh:
    with open(args.maf, 'rt') as handle:

        ## For each LCB in the alignment of consensus genome
        for rec in parse_maf(handle):
            label = rec['a']['label']
            length = len(rec['s'][0].seq)
            oh.write('a label=' + label + '\n')

            ## For each consensus sequence that contributes to this LCB
            for srec in rec['s']:
                lcb = srec.seq_name
                consensus = srec.seq_name.split('.')[0]
                coords = full_coords[consensus]

                ## For each individual sequence that is part of the consensus sequence
                for ind_seq_name in individual_contigs[lcb]:
                    direction = 1
                    j = srec.start
                    if srec.strand == '-':
                        direction = -1
                        j = srec.contig_length - srec.start - 1
                    ind_start = None
                    ind_end = None
                    ind_strand = None # Multiplied by direction of LCB in consensus alignment
                    ind_contig_length = contig_lengths[ind_seq_name]
                    aln_length = 0
                    sequence = ''

                    ## Loop through each base to figure out start, end,
                    ## strand, and aligned sequence
                    for i in range(length):
                        if srec.seq[i] == '-':
                            sequence += '-'
                        else:
                            ## Matching coordinates:
                            matching = coords[(lcb, j)][ind_seq_name]
                            pos = matching[0]
                            ind_strand = matching[1]
                            bp = matching[2]
                            if direction == -1: bp = Seq.reverse_complement(bp)
                            sequence += bp
                            if bp != '-': aln_length += 1
                            j += direction
                            if pos is not None and (ind_start is None or ind_start > pos):
                                ind_start = pos
                            if pos is not None and (ind_end is None or ind_end < pos):
                                ind_end = pos

                    ## Print out the record
                    actual_direction = ind_strand * direction
                    if len(set(sequence) - set('-')) == 0:
                        ## Covers entirely a gap region in the non-consensus alignment
                        continue
                    if actual_direction == 1:
                        oh.write('s\t' + ind_seq_name + '\t' + str(ind_start) + '\t' +
                                 str(aln_length) + '\t+\t' + str(ind_contig_length) +
                                 '\t' + sequence + '\n')
                    else:
                        oh.write('s\t' + ind_seq_name + '\t' + str(ind_contig_length - ind_end) + '\t' +
                                 str(aln_length) + '\t-\t' + str(ind_contig_length) +
                                 '\t' + sequence + '\n')

            oh.write('\n')

## 3. Merge all the individual components together if parallelized
