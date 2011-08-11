#!/usr/bin/env python
"""
gtfStopCodonMerger.py is an ugly, ugly kludge, and it bothers me.
dent earl, dearl (a) soe ucsc edu

Takes in a gtf file. Merges all stop_codon annotations into the
appropriate CDS and expands the CDS's region.

"""
from optparse import OptionParser
import os
import re
import sys
import signal # deal with broken pipes                        

signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

class Stop:
    """Used to keep track of stop_codon annotations
    """
    def __init__(self):
        self.begin = 0
        self.end = 0

def main():
    usage = ('usage: %prog < file.gtf\n\n'
             '%prog takes in a gtf from stdin and merges all stop_codon annotations into\n'
             'the appropriate CDS and exands the CDS\'s region.')
    parser = OptionParser(usage = usage)
    options, args = parser.parse_args()
    
    annots = {}
    # read in all data
    for s in sys.stdin:
        if not s:
            break
        s = s.rstrip()
        m = re.search('gene_id "(.*?)"', s)
        if m is not None:
            if m.group(1) in annots:
                annots[m.group(1)].append(s)
            else:
                annots[m.group(1)]=[s]
        else:
            print 'There is a line that doesn\'t have a gene_id'
            print s

    # process gene_id sets together
    newAnnots = {}
    for a in annots:
        stopCodon = Stop()
        # get the stop_codon
        for s in annots[a]:
            m = re.search('\tstop_codon\t(\d+)\t(\d+)\t', s)
            if m is not None:                
                stopCodon.begin = int(m.group(1))
                stopCodon.end   = int(m.group(2))
        if not stopCodon.begin:
            # this set of annotations did not have a stop_codon
            if a in newAnnots:
                newAnnots[a].append(s)
            else:
                newAnnots[a] = [s]
        else:
            tempRecord = []
            for s in annots[a]:
                m = re.search('\tCDS\t(\d+)\t(\d+)\t', s)
                if m:
                    if int(m.group(2)) == stopCodon.begin - 1:
                        tempRecord.append(re.sub('CDS\t(\d+)\t\d+\t', 
                                                 r'CDS\t\1\t' + str(stopCodon.end) + '\t', s))
                    else:
                        tempRecord.append(s)
                else:
                    m = re.search('\tstop_codon\t', s)
                    if not m:
                        tempRecord.append(s)
            newAnnots[a] = tempRecord
    for na in newAnnots:
        for a in newAnnots[na]:
            print a
                
if __name__ == "__main__":
    main()
