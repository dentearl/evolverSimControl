#!/usr/bin/env python
"""
eval_evolverpairwiseMAFextractor.py
dent earl, dearl (a) soe ucsc edu
16 nov 2009
A script that will take the path to a simulation
out directory and makes repeated calls to cvt
(in fabulous PARALLEL-vision!) in order
to extract FASTA files from either all cycles, or
only the leaves.
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
import re
import subprocess
import sys
from optparse import OptionParser
from sonLib.bioio import newickTreeParser
import xml.etree.ElementTree as ET
import evolverSimControl.lib.libSimControl as lsc

def initOptions(parser):
    parser.add_option('--simDir',dest='simDir',
                      help='Out directory from simTree.py.')
    parser.add_option('--allCycles', action='store_true', dest='allCycles',
                      default=False, help='Will extract fastas from all cycles, not just leafs. '
                      'default=%default')

def checkOptions(options, parser):
    if options.simDir is None:
        parser.error('specify --simDir.\n')
    if not os.path.isdir(options.simDir):
        parser.error('simulation directory "%s" not a directory!\n' % options.simDir)
    options.simDir = os.path.abspath(options.simDir)
    if not os.path.exists(os.path.join(options.simDir, 'simulationInfo.xml')):
        parser.error('unable to find %s.\n' % os.path.join(options.simDir, 'simulationInfo.xml'))
    infoTree = ET.parse(os.path.join(options.simDir, 'simulationInfo.xml'))
    treeObj = infoTree.find('tree')
    options.inputNewick=treeObj.text
    treeObj = infoTree.find('rootDir')
    options.rootDir=treeObj.text

def directoriesOnly(aList):
    """directoriesOnly() takes a list of items from a directory
    and purges anything from the list that is not a directory.
    """
    bList = []
    for i in aList:
        if os.path.isdir(i):
            bList.append(i)
    return bList

def extractLeaves(nt, leafDict):
    """Given a newick tree object, it returns a dict of
    leaf objects. Operates recursively.
    """
    if nt is None:
        return None
    nt.distance = 0
    if nt.right is None and nt.left is None:
        leafDict[nt.iD] = True
    else:
        extractLeaves(nt.right, leafDict = leafDict)
        extractLeaves(nt.left , leafDict = leafDict)    
    
def main():
    usage = ('usage: %prog --simDir path/to/dir [options]\n\n'
             '%prog takes in a simulation directory and then extracts the sequences\n'
             'of some of the cycles to their directories in fasta format.')
    parser = OptionParser(usage = usage)
    initOptions(parser)
    options, args = parser.parse_args()
    checkOptions(options, parser)
    
    cycles = glob.glob(os.path.join(options.simDir, '*'))
    cycles = directoriesOnly(cycles)
    leaves = {}
    nt = newickTreeParser(options.inputNewick, 0.0)
    extractLeaves(nt, leaves)
    for d in cycles:
        if not options.allCycles and not os.path.basename(d) in leaves:
            continue

        cmds = []
        outPipes = []
        inPipes = []
        nameA     = os.path.basename(d)
        nameA     = nameA.replace('[','')
        nameA     = nameA.replace(']','')
        cleanName = nameA.replace('\'','')
        
        cmd = [lsc.which('evolver_cvt')]
        cmd.append('-fromrev')
        cmd.append(os.path.join(d,'seq.rev'))
        cmd.append('-tofasta')
        cmd.append(os.path.join(d, 'seq.fa.tmp'))
        inPipes.append(None)
        outPipes.append(None)
        cmds.append(cmd)
        
        cmd = [lsc.which('mv')]
        cmd.append(os.path.join(d, 'seq.fa.tmp'))
        cmd.append(os.path.join(d, 'seq.fa'))
        inPipes.append(None)
        outPipes.append(None)
        cmds.append(cmd)
        
        cmd = [lsc.which('sed')]
        cmd.append(r"s/^>/>%s./;" % cleanName)
        inPipes.append(os.path.join(d, 'seq.fa'))
        outPipes.append(os.path.join(d, 'seq.name.fa.tmp'))
        cmds.append(cmd)
        
        cmd = [lsc.which('mv')]
        cmd.append(os.path.join(d, 'seq.name.fa.tmp'))
        cmd.append(os.path.join(d, 'seq.name.fa'))
        inPipes.append(None)
        outPipes.append(None)
        cmds.append(cmd)
        
        lsc.runCommands(cmds, os.curdir, outPipes = outPipes, inPipes = inPipes, mode = 's')
    
if __name__ == "__main__":
    main()
