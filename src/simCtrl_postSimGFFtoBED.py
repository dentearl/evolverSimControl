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
import sys
from optparse import OptionParser

def initOptions(parser):
    parser.add_option('--gff', dest = 'gff', 
                      help = 'gff to expand into beds')
    parser.add_option('--speciesName', dest = 'speciesName', type = 'string',
                      help = 'species name to append at the start of field one, i.e. "hg19.chr0 ..."')
    parser.add_option('--outDir', dest = 'outDir', default = os.curdir,
                      help = 'location to write out beds')
    parser.add_option('--outPrefix', dest = 'outPrefix', default = '',
                      help = 'prefix to put at front of file names, i.e. "hg19.annots.CDS.bed"')

def checkOptions(options, args, parser):
    if options.gff is None:
        parser.error('specify --gff')
    if not os.path.exists(options.gff):
        parser.error('--gff %s does not exist' % options.gff)

def closeFiles(fileHandleMap):
    for f in fileHandleMap:
        fileHandleMap[f].close()

def convertGFFtoBED(line, options):
    #      0       1             2                3        4           5      6      7
    # GFF: seqname source        feature          start(1) end(inclu.) score  strand frame
    # BED: chrom   chromStart(0) chromEnd(exclu.) name     score       strand 
    t = line.split('\t')
    bed=[]
    bed.append(t[0]) # chr name
    bed.append(str(int(t[3])-1)) # start pos
    bed.append(str(int(t[4]))) # end pos
    bed.append(t[2]) # feature name
    bed.append('0')
    bed.append(t[6])
    if options.speciesName:
        bed[0] = '%s.%s' % (options.speciesName, bed[0])
    bedLine = '\t'.join(bed)
    return '%s\n' % bedLine

def recordLineToFile(fileHandleMap, name, bedLine, options):
    try:
        fileHandleMap[name].write(bedLine)
    except KeyError:
        if options.outPrefix:
            s = '%s.annots.%s.bed' % (options.outPrefix, name)
        else:
            s = 'annots.%s.bed' % name
        fileHandleMap[name] = open(os.path.join(options.outDir, s), 'w')
        fileHandleMap[name].write(bedLine)

def main():
    parser = OptionParser()
    initOptions(parser)
    options, args = parser.parse_args()
    checkOptions(options, args, parser)
    
    fileHandleMap = {}
    
    f = open(options.gff, 'r')
    for line in f:
        line = line.strip()
        d = line.split('\t')
        if len(d) == 9:
            bedLine = convertGFFtoBED(line, options)
            recordLineToFile(fileHandleMap, d[2], bedLine, options)

    closeFiles(fileHandleMap)

if __name__ == "__main__":
    main()
