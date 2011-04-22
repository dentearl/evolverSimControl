#!/usr/bin/env python
# simUtil_maskedFAStats.py
# dent earl, dearl (a) soe ucsc edu
# 8 dec 2011
#
# scan through a FASTA, counting up the 
# masked and unmasked sequence characters
#
##############################
import os
import re
import sys

def usage():
    sys.stderr.write('USAGE: %s file.fa\n' %(sys.argv[0]))
    sys.exit( 2 )

def main():
    if len(sys.argv) != 2:
        usage()
    if not os.path.isfile( sys.argv[1] ):
        usage()
    infile = open( sys.argv[1], 'r' )
    DNA = re.compile('[ACTG]')
    dna = re.compile('[actg]')
    FASTALen = 0
    fastaLen = 0
    for line in infile:
        if line[0] == '>':
            continue
        FASTALen = FASTALen + len( DNA.findall( line ) )
        fastaLen = fastaLen + len( dna.findall( line ) )

    total = FASTALen + fastaLen
    print 'masked:   %9d (%6.2f%%)' % ( fastaLen, 100 * float(fastaLen) / total )
    print 'unmasked: %9d (%6.2f%%)' % ( FASTALen, 100 * float(FASTALen) / total )
    print 'total:    %9d' % total

if __name__ == "__main__":
    main()
