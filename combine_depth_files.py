#!/usr/bin/env python2.7
"""
    Combine tsv files with read depth. Files should have all the same
    lines and should have a column for chromosome/contig, a column for
    position, and a column for depth.

    script/bin/combine_depth_files.py <file(s)...>

    The directory name of each file is used as the strain name.
"""

import gzip
import itertools
import os.path as osp
import sys

handles = []
strains = []
for fname in sys.argv[1:]:
    strains.append(osp.basename(osp.dirname(fname)))
    if fname.endswith('.gz'):
        handles.append(gzip.open(fname))
    else:
        handles.append(open(fname))

sys.stdout.write('replicon\tpos\t' + '\t'.join(strains) + '\n')
for lines in itertools.izip(*handles):
    r = None
    p = None
    for i, line in enumerate(lines):
        fields = line.strip().split('\t')
        if i == 0:
            sys.stdout.write(fields[0] + '\t' + fields[1])
            r = fields[0]
            p = fields[1]
        else:
            if fields[0] != r or fields[1] != p:
                raise Exception('Unmatched position at %s in %s' % (line, strain))
        sys.stdout.write('\t' + fields[2])
    sys.stdout.write('\n')
