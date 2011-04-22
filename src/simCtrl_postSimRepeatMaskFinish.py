#!/usr/bin/env python
"""
postSimRepeatMaskWrapper.py
dent earl, dearl (a) soe ucsc edu
9 september 2010
a script to be run following the completion of a
simulation. The script will check the simulation dir
and for each cycle with a cycle/repeatMask/ dir it
will extract the repeat masked fasta file and stick
it at cycle/seq.masked.fa
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
import glob
import os
import subprocess
import sys
from sonLib.bioio import newickTreeParser
from sonLib.bioio import printBinaryTree
from optparse import OptionParser
import xml.etree.ElementTree as ET
import simulation.lib.libSimControl as LSC
import simulation.lib.libSimTree as LST

programs = ['twoBitToFa', 'parasol']
LSC.verifyPrograms(programs)
( TWOBIT2FA_BIN, PARASOL ) = programs

def usage():
    print 'USAGE: '+sys.argv[0]+' --simDir <dir>'
    print __doc__
    sys.exit(2)

def initOptions(parser):
    parser.add_option('-p', '--simDir',dest='simDir',
                      help='Simulation directory.')

def checkOptions(options):
    if (options.simDir == None):
        sys.stderr.write('%s: Error, specify simulation dir.\n' % sys.argv[0])
        usage()
    options.simDir  = os.path.abspath(options.simDir)
    if not os.path.exists(os.path.join(options.simDir, 'simulationInfo.xml')):
        sys.stderr.write('%s: Error, unable to find simulationInfo.xml.\n' % sys.argv[0])
        usage()
    infoTree = ET.parse(os.path.join(options.simDir, 'simulationInfo.xml'))
    treeObj = infoTree.find('tree')
    options.inputNewick=treeObj.text
    treeObj = infoTree.find('rootDir')
    options.rootDir=treeObj.text

def findDirectories( simDir ):
    allItems  = glob.glob(os.path.join(simDir, '*'))
    allDirs   = []
    readyDirs = []
    for a in allItems:
        if os.path.isdir(a):
            allDirs.append(a)
    for a in allDirs:
        if os.path.exists(os.path.join(a,'seq.name.fa')):
            readyDirs.append(a)
    return readyDirs

def main():
    parser=OptionParser()
    initOptions(parser)
    (options, args) = parser.parse_args()
    checkOptions(options)
    dirsToProcess = findDirectories( options.simDir )
    jobsArray = [] 
    cmdArray  = []
    for thisDir in dirsToProcess:
        if not os.path.exists( os.path.join( thisDir, 'repeatMask', 'seq.rmsk.2bit')  ):
            sys.stderr.write('%s: Error, file %s does not exist!\n' %(sys.argv[0], os.path.join(thisDir, 'repeatMask', 'seq.rmsk.2bit')))
            sys.exit( 1 )
        logPath = os.path.abspath('/dev/null')
        logFILE = open(logPath, 'w')
        convertCMD =  [ TWOBIT2FA_BIN ]
        convertCMD.append(os.path.join(thisDir, 'repeatMask', 'seq.rmsk.2bit') )
        convertCMD.append(os.path.join(thisDir, 'seq.masked.fa') )
        cmdArray.append( ' '.join( convertCMD ) )
        job = subprocess.Popen( convertCMD, stderr=logFILE, stdout=logFILE )
        jobsArray.append( job )
    i=-1
    for job in jobsArray:
        i=i+1
        job.wait()
        if(job.returncode):
            sys.stderr.write('%s: Error in a sub-processes "%s", returncode: %s.\n' %(sys.argv[0], cmdArray[i], job.returncode))
    

if __name__ == "__main__":
    main()


