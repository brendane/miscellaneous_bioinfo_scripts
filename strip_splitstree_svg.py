#!/usr/bin/env python2
"""
    Script to modify output from SplitsTree because I don't think any
    software exists to read the network format.

    This script takes input from the add_attr_splitstree_svg scripts
    and gets rid of the big rectangle around the figure and also the
    scale bar, which makes it more suitable for pasting into a file
    with other figures.

    strip_splitstree_svg.py <input file>
    
    Right now the script writes the file as HTML, so there is some manual
    editing required to turn it into a true svg file or piping through
    grep to get rid of the html and body tags.
"""

import math
import sys

from bs4 import BeautifulSoup

svgfile = sys.argv[1]

max_x = None
max_y = None
min_x = None
min_y = None

with open(svgfile, 'rb') as ihandle:
    soup = BeautifulSoup(ihandle, 'lxml')
    
    # Make the rectangles transparent; Getting rid of them
    # messes up some of the coordinates.
    rect_nodes = soup.find_all('rect')
    for rn in rect_nodes:
        rn.attrs['fill'] = 'none'

    # Find the bounding box of the leaf nodes
    circle_nodes = soup.find_all('circle')
    for cn in circle_nodes:
        x, y = float(cn.attrs['cx']), float(cn.attrs['cy'])
        if max_x is None or max_x < x:
            max_x = x
        if max_y is None or max_y < y:
            max_y = y
        if min_x is None or min_x > x:
            min_x = x
        if min_y is None or min_y > y:
            min_y = y

    # Get rid of any lines that are (nearly) outside of the bounding box
    line_nodes = soup.find_all('line')
    for ln in line_nodes:
        x1, x2 = float(ln.attrs['x1']), float(ln.attrs['x2'])
        y1, y2 = float(ln.attrs['y1']), float(ln.attrs['y2'])
        if max(x1, x2) > max_x or min(x1, x2) < min_x or \
            max(y1, y2) > max_y or min(y1, y2) < min_y:
            ln.extract()

    print soup.prettify()
