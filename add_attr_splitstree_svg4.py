#!/usr/bin/env python2
"""
    Script to modify output from SplitsTree because I don't think any
    software exists to read the network format.

    This version uses a diverging color scheme and is not hard-coded for
    a particular range of data values. Also, removes the text label for
    all branches, even if there isn't a phenotype value.

    Requirements:
        - Choose "simple" method of laying out the labels
        - Keep hitting Ctrl-minus until the font size is as small as it
          can go - this means there won't be much of a gap b/n the circles
          and the lines
        - Save image from SplitsTree as an SVG; this was tested with
          v4.14
        - Make a tab-delimited file with strain names in the first
          column and the attribute to use in the other column; attribute
          should be a number
        - BeautifulSoup v4 must be installed
    
    Output is to stdout.

    add_attr_splitstree_svg.py <input file> <attributes file>

    Alternatives:
        - Export network to GML format and plot with igraph - but this
          doesn't display right
    
    Right now the script writes the file as HTML, so there is some manual
    editing required to turn it into a true svg file or piping through
    grep to get rid of the html and body tags.
"""

import math
import sys

from bs4 import BeautifulSoup

def orange(x):
    r = int(255 - x * 120 + 50)
    g = int(150 - x * 150 + 25)
    b = int(60 - x * 60 + 10)
    return ((r, g, b), x)

def purple(x):
    x = 1 - x
    r = int(200 - x * 200 + 25)
    g = int(100 - x * 100 + 10)
    b = int(255 - x * 120 + 50)
    return ((r, g, b), x)

def color(x):
    if math.isnan(x):
        return ((255, 255, 255), 0.)
    elif x > 0.51:
        return orange(x)
    elif x < 0.49:
        return purple(x)
    else:
        return ((50, 50, 50), 100/255.)


svgfile = sys.argv[1]
attrfile = sys.argv[2]

strain_info = {}
min_attr = None
max_attr = None
with open(attrfile, 'rb') as ihandle:
    for line in ihandle:
        fields = line.split('\t')
        try:
            x = float(fields[1])
        except ValueError:
            strain_info[fields[0]] = float('nan')
            continue
        if min_attr is None or x < min_attr:
            min_attr = x
        if max_attr is None or x > max_attr:
            max_attr = x
        strain_info[fields[0]] = float(fields[1])

for strain, attr in strain_info.iteritems():
    x = strain_info[strain]
    strain_info[strain] = (x - min_attr) / (max_attr - min_attr)

with open(svgfile, 'rb') as ihandle:
    soup = BeautifulSoup(ihandle, 'lxml')
    text_nodes = soup.find_all('text')
    for tn in text_nodes:
        label = tn.text
        if label in strain_info:
            attr = strain_info[label]
            rgb, opacity = color(attr)
            col_string = '#' + ''.join(map(chr, rgb)).encode('hex')
            x = int(tn.attrs['x'])
            y = int(tn.attrs['y'])
            new_element = soup.new_tag('circle', cx=x, cy=y, r=3,
                                       style='fill:%s; stroke:%s; fill-opacity:%f; stroke-opacity:%f;' %
                                       (col_string, col_string, opacity, opacity))
            index = tn.parent.contents.index(tn)
            tn.parent.insert(index + 1, new_element)
        tn.extract()
    print soup.prettify()
