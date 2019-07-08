#!/usr/bin/env python3
"""
    Given a tab-delimited file with NCBI BLAST hit information, determine
    the name of the gene closest to the hit.

    find_closest_gene_from_blast_hits.py --email <email> --output <output file>
        <input file>

    The input file should be tab-delimited and have these columns (no
    header):
    1. subject NCBI accession number
    2. % identity
    3. subject start
    4. subject end
    5. subject scientific name

    Assumes that the input file contains only the top hits.
"""

import argparse
import collections
import csv
import os
import os.path as osp
import tempfile

from Bio import Entrez
from Bio.Entrez.Parser import DataHandler
from Bio import SeqIO

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--email')
parser.add_argument('--output')
parser.add_argument('infile')
args = parser.parse_args()

## Without some of these lines, a bunch of stuff ends up in
## ~/.config, which can be a problem. I think only the last
## line is necessary.
DataHandler.global_dtd_dir = osp.join(os.getcwd(), 'DTDs')
DataHandler.global_xsd_dir = osp.join(os.getcwd(), 'XSDs')
DataHandler.local_dtd_dir = os.getcwd()
DataHandler.local_xsd_dir = os.getcwd()
DataHandler.directory = os.getcwd()

Entrez.email = args.email

hits = collections.defaultdict(list)
seq_files = {}
with open(args.infile, 'rt') as ih:
    rdr = csv.reader(ih, delimiter='\t')
    for i, row in enumerate(rdr):
        gid = row[0].split('|')[-2]
        pid = float(row[1])
        s = int(row[2])
        e = int(row[3])
        n = row[4]

        ## Check if we have already encountered this id, if not,
        ## download the sequence
        if gid not in hits:
            seq = Entrez.efetch(db='nucleotide', id=gid,
                                rettype='gb', retmode='text').read()
            tmpfile = tempfile.mkstemp(dir='.')
            with open(tmpfile[0], 'wt') as th:
                th.write(seq)
            seq_files[gid] = tmpfile[1]

        ## Store information
        hits[gid].append((s, e, pid, n, gid))

        ## (Temporary) quit after 10 sequences
        if i == 9:
            break

## For each hit, find the closest gene
with open(args.output, 'wt') as oh:
    for gid, hit_info in hits.items():
        features = [f for f in SeqIO.read(seq_files[gid], 'genbank').features if len(f.qualifiers) > 2 and f.type != 'source']
        for s, e, i, sn, g in hit_info:
            closest = None
            mind = 100E6 # Just a number that is bigger than any bacterial genome
            for f in features:
                fe = int(f.location.end)
                fs = int(f.location.start)
                if fs <= s and fe >= e:
                    d = 0
                else:
                    d = min(abs(fs - s), abs(fe - e))
                if d < mind:
                    mind = d
                    closest = f
            annot = ''
            if closest is not None:
                try:
                    annot = closest.qualifiers['product']
                except KeyError:
                    annot = ''
            oh.write('\t'.join(str(x) for x in [sn, g, i, ';'.join(annot)]) + '\n')
