#!/usr/bin/env python2
"""
    Check a fastq file that's been through some processing against the
    input. Checks that:
        - all reads in the processed file(s) (by read.id) are also in the
          raw file
        - processed read length <= raw read length
        - each processed read file has no more than one copy of each read
          (this won't work with interleaved files)

    check_processed_fastq_by_raw <raw> <processed...>
"""

import gzip
import sys

from Bio import SeqIO

open_fun = open
if sys.argv[1].endswith('.gz'):
    open_fun = gzip.open
raw_reads = {}
with open_fun(sys.argv[1], 'rb') as ihandle:
    for rec in SeqIO.parse(ihandle, 'fastq'):
        raw_reads[rec.id] = len(rec)

for pfile in sys.argv[2:]:
    observed = set()     # Only check within a file - might be paired-end reads
    open_fun = open
    if pfile.endswith('.gz'):
        open_fun = gzip.open
    with open_fun(pfile, 'rb') as ihandle:
        for rec in SeqIO.parse(ihandle, 'fastq'):
            if rec.id not in raw_reads:
                import pdb; pdb.set_trace()
                print '%s not found - fail on %s' %(rec.id, pfile)
                exit(1)
            if len(rec) > raw_reads[rec.id]:
                print '%s too big - fail on %s' %(rec.id, pfile)
                exit(1)
            if rec.id in observed:
                print '%s duplicated - fail on %s' %(rec.id, pfile)
                exit(1)
            observed.add(rec.id)

print 'pass'
exit(0)
