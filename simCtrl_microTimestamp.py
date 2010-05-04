#!/usr/bin/env python
"""
microTimestamp.py
 dent earl, dearl (a) soe ucsc edu

microTimestamp.py is a small helper script that takes in
a cycle directory that has completed some small part of a cycle.
By way of example, and the primary motivation for the writing of
the script, cycleMain_2 requires each chromosome to go through a
serial run of evolver_evo, evolver_cvt, and trf. This step takes
about 80% of the time of a cycle. 
Timestamps are placed in the cycleInfo.xml, under
info
 -timestamps
  -micro
   -yourTimeStamp

23 April 2010
"""
######################################
import xml.etree.ElementTree as ET
import os, sys, time
from optparse import OptionParser
import eval.lib.libSimControl as LSC
#import eval.lib.libSimCycle   as LSY

def usage():
    print "USAGE: %s --cycleDir [cycleDirectory/] --name [cycleStep_2_cycleMain_2_evolver_evo_start]" %(sys.argv[0])
    print __doc__
    sys.exit(2)

def performTimestamp(file, stamp, chrom=''):
    """performTimestamp() takes a micro.CHR.info.xml file and a string called stamp
    and places the standard timestamp (human readable + epoch) under the
    hierarchy info:timestamps:micro:stamp
    """
    if os.path.exists(file):
        infoTree=ET.parse(file)
        root=infoTree.getroot()
        tObj=root.find('timestamps')
        mObj=tObj.find('micro')
        if mObj == None:
            # create new tag
            mObj=ET.SubElement(tObj, 'micro')
    else:
        root=ET.Element('info')
        tObj=ET.SubElement(root,'timestamps')
        mObj=ET.SubElement(tObj,'micro')
    
    if chrom:
        cObj=mObj.find(chrom)
        if cObj == None:
            cObj=ET.SubElement(mObj, chrom)
        nObj=ET.SubElement(cObj, stamp)
    else:
        nObj=ET.SubElement(mObj, stamp)
    timeHuman=ET.SubElement(nObj, 'humanUTC')
    timeHuman.text=str(time.strftime("%a, %d %b %Y %H:%M:%S (UTC)", time.gmtime()))
    timeEpoch=ET.SubElement(nObj, 'epochUTC')
    timeEpoch.text=str(time.time())
    info=ET.ElementTree(root)
    info.write(file)

def main(argv):
    parser=OptionParser()
    parser.add_option('-c', '--cycleInfo',dest='cycleInfo', type='string',
                      help='cycleInfo.xml file.')
    parser.add_option('--chr',dest='chrom', type='string',
                      help='Chromosome name.')
    parser.add_option('-n', '--name',dest='timeName', type='string',
                      help='Timestamp name.')
    (options, args) = parser.parse_args()
    if (options.cycleInfo == None):
        sys.stderr.write('%s: Error, specify --cycleInfo .\n' % sys.argv[0])
        usage()
    if (options.timeName == None):
        sys.stderr.write('%s: Error, specify name of timestamp.\n' % sys.argv[0])
        usage()
    performTimestamp(os.path.abspath(options.cycleInfo), options.timeName, options.chrom)


if __name__ == "__main__":
    main(sys.argv)
