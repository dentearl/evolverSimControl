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
##############################
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
