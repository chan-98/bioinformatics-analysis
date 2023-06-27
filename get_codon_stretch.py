#!/usr/bin/env python3

"""
This code counts the number of stretches of any codon present in a multi-fasta file.
Usage: get_codon_stretch.py sequence.fasta <pattern>
It prints the sequence ID, start position, end position and num of repeats of the pattern to std_out.
"""

import sys
from Bio import SeqIO

def getPatternStretches(seq, pattern, index):
        startIndex = 0
        stretch = ""
        stretchesList = []

        k = len(pattern)
        n = len(seq)

        for i in range(index, n, k):
                if seq[i:i+k] == pattern:
                        if startIndex == 0:
                                startIndex = i
                        stretch = stretch + seq[i:i+k]
                else:
                        if stretch:
                                stretchesList.append((startIndex, stretch))
                        stretch = ""
                        startIndex = 0
        return stretchesList

if len(sys.argv) != 3:
        print("USAGE:\n\tpython3 sys.argv[0] [fasta file] [codon] > [output].tsv")
else:
        file = sys.argv[1]
        pattern = sys.argv[2]
        k = len(pattern)
        print("Seq-ID\tStart-index\tEnd-index\tNumber-of-repeats")
        for record in SeqIO.parse(file, "fasta"):
                seqID = str(record.id)
                seq = str(record.seq)

                stretchAt0 = getPatternStretches(seq, pattern, 0)
                stretchAt1 = getPatternStretches(seq, pattern, 1)
                stretchAt2 = getPatternStretches(seq, pattern, 2)
                stretch = stretchAt0 + stretchAt1 + stretchAt2

                for (start, stretch) in stretch:
                        end = start + len(stretch)
                        repeats = int(len(stretch)/3)
                        print(seqID + "\t" + str(start) + "\t" + str(end) + "\t" + str(repeats))
