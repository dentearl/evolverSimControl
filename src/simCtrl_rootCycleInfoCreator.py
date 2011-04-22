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
