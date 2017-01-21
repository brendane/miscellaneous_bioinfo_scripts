#!/usr/bin/env python2
"""
    Adjust gene IDs for PopGenome to play nicely with my other code.
"""

import re
import sys

with open(sys.argv[1], 'rb') as handle:
    for line in handle:
        if line.startswith('#'):
            sys.stdout.write(line)
            continue
        fields = line.strip().split('\t')
        tags = fields[8]
        gid = re.sub(';.+', '', re.sub('ID=', '', tags))
        rest = re.sub('ID='+gid+';', '', tags)
        n = int(re.sub('.+\.', '', gid))
        rid = re.sub('\..+', '', gid)
        n += 1
        #n -= 1
        gid = rid + '.' + str(n)
        sys.stdout.write('\t'.join(fields[:8]) + '\tID=' + gid + ';' +
                         rest + '\n')
