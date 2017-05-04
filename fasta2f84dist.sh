#!/bin/bash
#
# Use phylip dnadist to create a distance matrix from a fasta file. Will
# overwrite anything called "infile" or "outfile" in the working directory.
#


# fasta -> phylip
convert="""
import sys
from Bio import AlignIO
aln = AlignIO.read(sys.argv[1], 'fasta')
AlignIO.write(aln, sys.stdout, 'phylip')
"""

fasta2f84dist() {

    # Convert to phylip format
    python -c "$convert" "$1" > "infile" \
        || { echo "converting to phylip format failed"; exit 1; }

    # Make phylip control file
    # Just use defaults; watch out for the interleaved formats
    # option.
    pc="phylip_control.txt"
    rm -f "outfile"
    rm -f $pc
    echo "Y" >> $pc
    dnadist < $pc \
        || { echo "running dnadist failed"; exit 1; }
    mv "outfile" "$2"
}

fasta2f84dist "$1" "$2"
