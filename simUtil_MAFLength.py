#!/usr/bin/env python
# simUtil_MAFLength.py
# dent earl, dearl (a) soe ucsc edu
# 16 dec 2009
#
# scan through a MAF, summing up the largest
# sequence in each block e.g.:
#
# a blah:
# s aoue aoeu 85 + ATGC...
# s aoue aoeu 87 + ATGC...
#
# a bleh:
# s aoeu aoeu 4 - ATGC 
# s aoeu aoeu 5 - ATGC
#
# yields a total MAF length of 87+5 = 92
#
##############################
import os
import re
import sys

def usage():
    sys.stderr.write('USAGE: %s file.maf\n' %(sys.argv[0]))
    sys.exit(2)

def main():
    if len(sys.argv) != 2:
        usage()
    if not os.path.isfile( sys.argv[1] ):
        usage()
    infile = open( sys.argv[1], 'r' )
    patLen = re.compile('^s\s+\S+\s+\d+\s+(\d+)\s+[-+]')
    blockLen = 0
    alnLen = 0
    while infile:
        line = infile.readline()
        if not line:
            break
        p=patLen.match(line)
        if p:
             blockLen = max( blockLen, int(p.group(1)) )
        else:
            if blockLen != 0:
                alnLen  += blockLen
                blockLen = 0
    if(blockLen != 0):
        alnLen += blockLen
    print alnLen

if __name__ == "__main__":
    main()
