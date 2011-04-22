#!/usr/bin/env python
"""
completeTimestamp.py
 dent earl, dearl (a) soe ucsc edu

completeTimestamp.py is a small helper script that takes in
a cycle directory that has just completed either the cycle
steps {begin, mid_1, mid_2, end} or has just completed the
stats steps {stats_1, stats_2, stats_3, stats_4} and then
updates the cycleInfo.xml file with the appropriate
timestamp.

16 March 2010
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
import os, sys
from optparse import OptionParser
import simulation.lib.libSimControl as LSC

def usage():
    print "USAGE: %s --cycleDir [cycleDirectory/] --type [main,stats]" %(sys.argv[0])
    print __doc__
    sys.exit(2)

def main(argv):
    parser=OptionParser()
    parser.add_option('-c', '--cycleDir',dest='cycleDir', type='string',
                      help='Cycle directory.')
    parser.add_option('-t', '--timeType',dest='timeType', type='string',
                      help='Timestamp type, either "main" or "stats".')
    (options, args) = parser.parse_args()
    if (options.cycleDir == None):
        sys.stderr.write('%s: Error, specify cycle dir.\n' % sys.argv[0])
        usage()
    if (not os.path.isdir(options.cycleDir)):
        sys.stderr.write('%s: Error, cycle dir "%s" not a directory!\n' % (sys.argv[0], options.cycleDir))
        usage()
    if (options.timeType == None):
        sys.stderr.write('%s: Error, specify timeType.\n' % sys.argv[0])
        usage()
    if (options.timeType != 'main') and (options.timeType != 'stats'):
        sys.stderr.write('%s: Error, timeType "%s" not in set {main, stats}. Choose one.\n' % (sys.argv[0], options.timeType))
        usage()

    typeDict={'main':'cycleStep_4_cycleMain_4_end', 'stats':'cycleStep_8_cycleStats_4_end'}
    # the cycleStep needs to be stamped...
    LSC.subTypeTimestamp(os.path.join(options.cycleDir,'cycleInfo.xml'),
                             options.timeType, typeDict[options.timeType])
    # and the entire cycle type needs to be stamped
    LSC.typeTimestamp(os.path.join(options.cycleDir,'cycleInfo.xml'),
                      options.timeType, 'End')
    

if __name__ == "__main__":
    main(sys.argv)
