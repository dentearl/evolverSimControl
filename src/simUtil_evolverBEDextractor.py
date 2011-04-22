#!/usr/bin/env python
"""
eval_evolverBEDextractor.py
dent earl, dearl (a) soe ucsc edu
2 nov 2010
A script that will take the path to an evolver
generated annots.gff annotation file and extracts
the annotations annotation-specific BED files.
Six files are produced:
annots.CDS.bed
annots.NGE.bed
annots.NXE.bed
annots.UTR.bed
annots.island.bed
annots.tandem.bed
"""
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
import os
import re
import sys
from optparse import OptionParser

def usage():
    sys.stderr.write('USAGE: %s  < annots.gff\n' % (sys.argv[0]))
    sys.exit(2)

def initOptions(parser):
    pass

def checkOptions(options):
    pass

def openFiles():
    fileNameMap = {'CDS':'annots.CDS.bed',
                   'NGE':'annots.NGE.bed',
                   'NXE':'annots.NXE.bed',
                   'UTR':'annots.UTR.bed',
                   'island':'annots.island.bed',
                   'tandem':'annots.tandem.bed'}
    fileHandleMap = {}
    for f in fileNameMap:
        fileHandleMap[f] = open( os.path.abspath( fileNameMap[f] ), 'w' )
    return fileHandleMap

def closeFiles( fileHandleMap ):
    for f in fileHandleMap:
        fileHandleMap[f].close()

def convertGFFtoBED( line ):
    #      0       1             2                3        4           5      6      7
    # GFF: seqname source        feature          start(1) end(inclu.) score  strand frame
    # BED: chrom   chromStart(0) chromEnd(exclu.) name     score       strand 
    t = line.split('\t')
    bed=[]
    bed.append( t[0] ) # chr name
    bed.append( str( int(t[3])-1 ) ) # start pos
    bed.append( str( int(t[4]) ) ) # end pos
    bed.append( t[2] ) # feature name
    bed.append( '0' )
    bed.append( t[6] )
    bedLine = '\t'.join( bed )
    return '%s\n' % bedLine

def main():
    parser=OptionParser()
    initOptions(parser)
    (options, args) = parser.parse_args()
    checkOptions(options)
    
    pat = re.compile('^.*?\t.*?\t(.*?)\t\d+')

    fileHandleMap = openFiles()
    
    for line in sys.stdin:
        r = re.match(pat, line)
        if r:
            bedLine = convertGFFtoBED( line )
            fileHandleMap[ r.group(1) ].write( bedLine )

    closeFiles( fileHandleMap )

if __name__ == "__main__":
    main()
