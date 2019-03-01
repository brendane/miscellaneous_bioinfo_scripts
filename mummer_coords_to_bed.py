#!/usr/bin/env python3
"""
    Convert a coords file to bed format by using the reference positions.

    mummer_coords_to_bed.py --output <output file> [--snp] [--no-header]
        <input file>

    --no-header: File not have header (not implemented yet)
    --target: report only (best guess at) target SNP postion. SNP
        positions are assumed to be in the name of the query like so:
        <contig>_<position>
"""

import argparse
import csv
import re

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--target', default=False, action='store_true')
parser.add_argument('--no-header', default=False, action='store_true')
parser.add_argument('input')
args = parser.parse_args()

with open(args.output, 'w') as oh:
    with open(args.input, 'r') as ih:
        ih.readline(); ih.readline; ih.readline(); ih.readline(); ih.readline()
        rdr = csv.reader(ih, delimiter='\t')
        for row in rdr:
            qry = row[7]
            ref = row[8]
            pid = float(row[6])
            reflen = int(row[5])
            qlen = int(row[4])
            start = int(row[2])
            end = int(row[3])
            
            if args.target:
                original_target_coord = int(re.sub('.+_([0-9]+)::.+', '\\1', qry))
                st, et = (int(x) for x in re.sub('.+::.+:', '', qry).split('-'))
                sq, eq = int(row[0]), int(row[1])
                target_in_qry = original_target_coord - st + 1
                if target_in_qry >= sq and target_in_qry <= eq:
                    ## Find the proportional distance between start and end
                    p = (target_in_qry - sq) / float(eq - sq + 1)
                    ref_target_coord = p * reflen + start
                    oh.write(ref + '\t' + str(int(ref_target_coord-1)) + '\t' +
                             str(int(ref_target_coord)) + '\t' + qry + '\t' +
                             str(qlen) + '\t' + str(pid) + '\n')
                else:
                    continue
            else:
                oh.write(ref + '\t' + str(start-1) + '\t' + str(end) +
                         '\t' + qry + '\t' +
                         str(qlen) + '\t' + str(pid) + '\n')

