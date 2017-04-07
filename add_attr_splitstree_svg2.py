#!/usr/bin/env python2
"""
    Script to modify output from SplitsTree because I don't think any
    software exists to read the network format properly.

    This version classifies strains into increase, decrease, and same.

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

import sys

from bs4 import BeautifulSoup

increase_color = '#edf8b1'
neutral_color = '#7fcdbb'
decrease_color = '#2c7fb8'

svgfile = sys.argv[1]
attrfile = sys.argv[2]
cutoff = float(sys.argv[3])

strain_info = {}
with open(attrfile, 'rb') as ihandle:
    for line in ihandle:
        fields = line.split('\t')
        strain_info[fields[0]] = float(fields[1])

with open(svgfile, 'rb') as ihandle:
    soup = BeautifulSoup(ihandle, 'lxml')
    text_nodes = soup.find_all('text')
    for tn in text_nodes:
        label = tn.text
        if label in strain_info:
            attr = strain_info[label]
            if attr > cutoff:
                col_string = increase_color
            elif attr < (-1 * cutoff):
                col_string = decrease_color
            else:
                col_string = neutral_color
            x = int(tn.attrs['x'])
            y = int(tn.attrs['y'])
            new_element = soup.new_tag('circle', cx=x, cy=y, r=2.5,
                    style='fill: %s; stroke: %s' % (col_string, col_string))
            index = tn.parent.contents.index(tn)
            tn.parent.insert(index + 1, new_element)
            tn.extract()
    print soup.prettify()
