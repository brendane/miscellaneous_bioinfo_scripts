#!/usr/bin/env python3
"""
    Convert a blast outfmt 6 file to bed format by using the reference positions.
    Assumes file is sorted with top hit first.

    blast_best_hit_to_bed.py --output <output file> [--target] <input file>

    --target: report only (best guess at) target SNP postion. SNP
        positions are assumed to be in the name of the query like so:
        <contig>_<position>
"""

import argparse
import csv
import itertools
import re

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--target', default=False, action='store_true')
parser.add_argument('input')
args = parser.parse_args()

with open(args.output, 'w') as oh:
    with open(args.input, 'r') as ih:
        rdr = csv.reader(ih, delimiter='\t')
        for qry, rows in itertools.groupby(rdr, lambda x:x[0]):
            for row in rows:
                ref = row[1]
                pid = float(row[2])
                start = int(row[8])
                end = int(row[9])
                reflen = end - start + 1
                qlen = int(row[7]) - int(row[6]) + 1
                
                if args.target:
                    original_target_coord = int(re.sub('.+_([0-9]+)::.+', '\\1', qry))
                    st, et = (int(x) for x in re.sub('.+::.+:', '', qry).split('-'))
                    sq, eq = int(row[6]), int(row[7])
                    target_in_qry = original_target_coord - st + 1
                    if target_in_qry >= sq and target_in_qry <= eq:
                        ## Find the proportional distance between start and end
                        p = (target_in_qry - sq) / float(eq - sq + 1)
                        ref_target_coord = p * reflen + start - 1
                        oh.write(ref + '\t' + str(int(ref_target_coord-1)) + '\t' +
                                 str(int(ref_target_coord)) + '\t' + qry + '\t' +
                                 str(qlen) + '\t' + str(pid) + '\n')
                        break
                    else:
                        continue
                else:
                    oh.write(ref + '\t' + str(start-1) + '\t' + str(end) +
                             '\t' + qry + '\t' +
                             str(qlen) + '\t' + str(pid) + '\n')
                    break

