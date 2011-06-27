#!/usr/bin/env python
"""
postSimMAFextractor.py
dent earl, dearl (a) soe ucsc edu
15 april 2010
a script to be run following the completion of a
simulation. The script will check the simulationInfo.xml
file and then figure out the order to both
extract the relvant node's MAFs but then combine
the MAFs into larger alignments in the correct order.
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
import evolverSimControl.lib.libSimControl as lsc
import evolverSimControl.lib.libSimControlClasses as lscc
from jobTree.scriptTree.stack import Stack
from optparse import OptionParser
import os
import re
from sonLib.bioio import newickTreeParser
from sonLib.bioio import printBinaryTree
import sys
import xml.etree.ElementTree as ET

programs = ['evolver_cvt', 'evolver_transalign', 'mafJoin']
            
lsc.verifyPrograms(programs)

def initOptions(parser):
    parser.add_option('--simDir', dest = 'simDir',
                      help = 'Simulation directory.')
    parser.add_option('--maxBlkWidth',dest = 'maxBlkWidth',
                      default = 10000, type = 'int',
                      help = ('Maximum mafJoin maf block output size. May be reduced '
                              'towards 250 for complicated phylogenies. '
                              'default=%default'))
    parser.add_option('--maxInputBlkWidth', dest = 'maxInputBlkWidth',
                      default = 1000, type = 'int',
                      help = ('Maximum mafJoin maf block input size. mafJoin will cut '
                              'inputs to size, may result in long runs for very simple '
                              'joins. May be reduced towards 250 for complicated '
                              'phylogenies. default=%default'))
    parser.add_option('--noBurninMerge', dest = 'noBurninMerge',
                      action='store_true', default = False,
                      help = ('Will not perform a final merge of simulation '
                              'to the burnin. default=%default'))

def checkOptions(options, parser):
    if options.simDir is None:
        parser.error('specify --simDir.\n')
    if not os.path.exists(options.simDir):
        parser.error('Directory does not exist: --simDir %s\n' % options.simDir)
    if not os.path.isdir(options.simDir):
        parser.error('Option --simDir %s is not a directory.\n' % options.simDir)
    options.simDir  = os.path.abspath(options.simDir)
    if options.jobTree is None:
        parser.error('specify --jobTree.\n')
    if os.path.exists(options.jobTree):
        parser.error('--jobTree %s already exists! If your simulation crashed, '
                    'relaunch it with "jobTreeRun --jobTree %s/" \n' % 
                    (options.jobTree, options.jobTree))
    if not os.path.exists(os.path.join(options.simDir, 'simulationInfo.xml')):
        parser.error('unable to find simulationInfo.xml in --simDir %s\n' % options.simDir)
    infoTree = ET.parse(os.path.join(options.simDir, 'simulationInfo.xml'))
    treeObj = infoTree.find('tree')
    options.inputNewick=treeObj.text
    rootNameObj = infoTree.find('rootDir')
    options.rootName = os.path.basename(rootNameObj.text)
    options.rootDir  = os.path.abspath(rootNameObj.text)

def launchJobTree(options):
    jobResult = Stack(lscc.MasterMafGenerator(options)).startJobTree(options)
    if jobResult:
        raise RuntimeError('The jobTree contained %d failed jobs!\n' % jobResult)

def main():
    usage=('usage: %prog --simDir path/to/dir [options]\n\n'
           '%prog requires mafJoin which is part of mafTools and is available\n'
           'at https://github.com/dentearl/mafTools/ . ')
    parser = OptionParser(usage = usage)
    initOptions(parser)
    Stack.addJobTreeOptions(parser)
    options, args = parser.parse_args()
    
    checkOptions(options, parser)
    
    launchJobTree(options)
    
if __name__ == "__main__":
    main()
