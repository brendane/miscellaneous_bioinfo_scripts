#!/usr/bin/env python2.7
"""
    Count the number of reads aligned in two bam files.

    compare_aligned_reads.py <bam> <bam>

    Writes to stdout.

    Uses the "AS" tag, which stands for "alignment score" to determine
    which file has the better mapping for a read. This tag might not be
    produced by all aligners.
"""

#==============================================================================#

import sys

import pysam as ps

#==============================================================================#

bam_file_1 = sys.argv[1]
bam_file_2 = sys.argv[2]

reads_1 = {}
reads_2 = {}
for fname, rd_dict in [(bam_file_1, reads_1), (bam_file_2, reads_2)]:
    for read in ps.Samfile(fname, 'rb'):
        if read.is_unmapped:
            continue
        a = [y for x, y in read.tags if x == 'AS'][0]
        name = read.qname + '_' + ('1' if read.is_read1 else '2')
        if name in rd_dict:
            old_a = rd_dict[name]
            if a > old_a:
                rd_dict[name] = a
        else:
            rd_dict[name] = a

n_mapped_both = 0
n_mapped_better_1 = 0
n_mapped_better_2 = 0
n_mapped_equal = 0
for read, m1 in reads_1.iteritems():
    if read in reads_2:
        n_mapped_both += 1
        m2 = reads_2[read]
        if m2 > m1:
            n_mapped_better_2 += 1
        elif m1 > m2:
            n_mapped_better_1 += 1
        else:
            n_mapped_equal += 1

print 'File 1\t' + bam_file_1
print 'File 2\t' + bam_file_2
print 'Mapped in 1\t' + str(len(reads_1))
print 'Mapped in 2\t' + str(len(reads_2))
print 'Mapped in both\t' + str(n_mapped_both) + '\t' + \
        str(n_mapped_both / float(len(reads_1)) * 100) + '%\t' + \
        str(n_mapped_both / float(len(reads_2)) * 100) + '%'
print 'Mapped equally in both\t' + str(n_mapped_equal) + '\t' + \
        str(n_mapped_equal / float(len(reads_1)) * 100) + '%\t' + \
        str(n_mapped_equal / float(len(reads_2)) * 100) + '%'
print 'Mapped better in 1\t' + str(n_mapped_better_1) + '\t' + \
        str(n_mapped_better_1 / float(len(reads_1)) * 100) + '%\t' + \
        str(n_mapped_better_1 / float(len(reads_2)) * 100) + '%'
print 'Mapped better in 2\t' + str(n_mapped_better_2) + '\t' + \
        str(n_mapped_better_2 / float(len(reads_1)) * 100) + '%\t' + \
        str(n_mapped_better_2 / float(len(reads_2)) * 100) + '%'
print 'Mapped as well or better in 2\t' + str(n_mapped_better_2 + n_mapped_equal) + '\t' + \
        str((n_mapped_better_2 + n_mapped_equal) / float(len(reads_1)) * 100) + '%\t' + \
        str((n_mapped_better_2 + n_mapped_equal) / float(len(reads_2)) * 100) + '%'
