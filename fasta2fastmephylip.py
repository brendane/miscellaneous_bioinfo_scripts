#!/usr/bin/env python3
"""
    Apparently FastME needs phylip format, even though it was written
    in 2016. So, here is a format converter.

    fasta2fastmephylip.py <input>
"""

import re
import sys

data = []
with open(sys.argv[1], 'rt') as ih:
    name = None
    seq = []
    for line in ih:
        if line.startswith('>'):
            if name is not None:
                data.append((name, ''.join(seq)))
            name = line.strip()[1:]
            seq = []
            if len(name) > 64:
                raise Exception('%s is too long' % name)
        else:
            seq.append(re.sub('\s', '', line.strip()))
    if name is not None:
        data.append((name, ''.join(seq)))

sys.stdout.write(str(len(data)) + ' ' + str(len(data[0][1])) + '\n')
for n, s in data:
    sys.stdout.write(n + ' ' + s + '\n')
