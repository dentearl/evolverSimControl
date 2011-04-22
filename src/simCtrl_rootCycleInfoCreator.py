#!/usr/bin/env python
"""
 dummyCycleInfoCreator.py
 dent earl, dearl (a) soe ucsc edu
 16 april 2010
 When the user supplies simTree.py with a newick tree that starts
 at a branch point, the simulation copies in the simulation root into
 a directory named 'root' and extends the newick tree to include this
 root as the branch point. For everything to work soomthly the root dir
 needs to have a cycleInfo.xml file. This script creates that file.
"""
##############################
import xml.etree.ElementTree as ET
from optparse import OptionParser
import os, sys

def usage():
    print 'USAGE: %s --dir <dir>' %(sys.argv[0])
    print __doc__
    sys.exit(2)

def initOptions(parser):
    parser.add_option('-d', '--dir',dest='dir',
                      help='Destination directory.')

def checkOptions(options):
    if not os.path.isdir(options.dir):
        sys.stderr.write('%s: Error, directory "%s" is not a directory!\n' % (sys.argv[0], options.dir))
        usage()
    options.dir=os.path.abspath(options.dir)
    

def main():
    parser=OptionParser()
    initOptions(parser)
    (options, args) = parser.parse_args()
    checkOptions(options)
    root=ET.Element('info')
    e=ET.SubElement(root, 'cycleIsRoot')
    e.text=str(True)
    info=ET.ElementTree(root)
    info.write(os.path.join(options.dir,'cycleInfo.xml'))
    

if __name__ == "__main__":
    main()
