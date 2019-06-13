#!/usr/bin/env python3
"""
    Parse the MUMmer delta file format into a tab-delimited file with
    matching positions.

    parse_mummer_delta.py [--show-gaps] <delta file>

    --show-gaps     List gaps, using 0 as the position number
"""

import argparse
import sys

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--show-gaps', default=False, action='store_true')
parser.add_argument('delta')
args = parser.parse_args()

with open(args.delta, 'rt') as ih:
    f1, f2 = ih.readline().strip().split()
    ih.readline()
    for line in ih:
        if line.startswith('>'):
            fields = line.strip()[1:].split()
            ref_contig = fields[0]
            qry_contig = fields[1]
            ref_contig_len = int(fields[2])
            qry_contig_len = int(fields[3])
        else:
            fields = line.strip().split()
            if len(fields) == 1:
                dist_to_indel = int(fields[0])
                for i in range(abs(dist_to_indel)-1):
                    sys.stdout.write(qry_contig + '\t' + str(qry_move+qry_start) + '\t' +
                                     ref_contig + '\t' + str(ref_move+ref_start) + '\t' +
                                     strand + '\n')
                    ref_move += ref_dir
                    qry_move += qry_dir
                if dist_to_indel > 0:
                    if args.show_gaps:
                        sys.stdout.write(qry_contig + '\t' + '0' + '\t' +
                                         ref_contig + '\t' + str(ref_move+ref_start) + '\t' +
                                         strand + '\n')
                    ref_move += ref_dir
                elif dist_to_indel < 0:
                    if args.show_gaps:
                        sys.stdout.write(qry_contig + '\t' + str(qry_move+qry_start) + '\t' +
                                         ref_contig + '\t' + '0' + '\t' +
                                         strand + '\n')
                    qry_move += qry_dir
                else:
                    for i in range(abs(ref_end - (ref_move*ref_dir+ref_start))+1):
                        sys.stdout.write(qry_contig + '\t' + str(qry_move+qry_start) + '\t' +
                                         ref_contig + '\t' + str(ref_move+ref_start) + '\t' +
                                         strand + '\n')
                        ref_move += ref_dir
                        qry_move += qry_dir
            else:
                ref_start, ref_end, qry_start, qry_end, _, _, _ = \
                        (int(x) for x in line.split())
                ref_strand = 1; ref_dir = 1; qry_strand = 1; qry_dir = 1
                ref_move = 0; qry_move = 0
                if ref_start > ref_end:
                    ref_strand = -1
                    ref_dir = -1
                    #ref_move = -1
                if qry_start > qry_end:
                    qry_strand = -1
                    qry_dir = -1
                    #qry_move = -1
                if qry_strand == ref_strand:
                    strand = '+'
                else:
                    strand = '-'
