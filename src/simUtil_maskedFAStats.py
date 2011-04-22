#!/usr/bin/env python
# simUtil_maskedFAStats.py
# dent earl, dearl (a) soe ucsc edu
# 8 dec 2011
#
# scan through a FASTA, counting up the 
# masked and unmasked sequence characters
#
##################################################
# Copyright (C) 2009-2011 by
# Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
# Benedict Paten (benedict@soe.ucsc.edu, benedictpaten@gmail.com)
# ... and other members of the Reconstruction Team of David Haussler's
# lab (BME Dept. UCSC)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
##################################################
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
